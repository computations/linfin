#include <CLI/App.hpp>
#include <CLI/Config.hpp>
#include <CLI/Formatter.hpp>
#include <algorithm>
#include <cmath>
#include <corax/corax.hpp>
#include <cstdint>
#include <cstring>
#include <expected>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <logger.hpp>
#include <string>
#include <string_view>
#include <utility>
#include <yaml-cpp/yaml.h>

using Path = std::filesystem::path;
using Tree = corax_utree_t;
using std::expected;
using std::unexpected;

inline uint32_t bextr(uint32_t a, size_t i) { return (a >> i) & 1ULL; }

inline uint32_t popcount(uint32_t a) { return __builtin_popcount(a); }

inline uint32_t find_first_set(uint32_t a) { return __builtin_ctz(a); }

class AccumulationTable {
public:
  using AccumulationType      = uint32_t;
  static constexpr size_t max = std::numeric_limits<uint32_t>::max();

  AccumulationTable(size_t lineages, size_t queries)
      : _table{(AccumulationType *)calloc(lineages * queries,
                                          sizeof(AccumulationType))},
        _lineage_count{lineages},
        _query_count{queries} {}

  ~AccumulationTable() {
    if (_table) { free(_table); }
  }

  AccumulationTable(const AccumulationTable &)            = delete;
  AccumulationTable &operator=(const AccumulationTable &) = delete;

  AccumulationTable(AccumulationTable &&other)
      : _table{other._table},
        _lineage_count{other._lineage_count},
        _query_count{other._query_count} {
    other._table = nullptr;
  }

  AccumulationTable &operator=(AccumulationTable &&other) {
    std::swap(*this, other);

    if (other._table) { free(_table); }
    return *this;
  }

  AccumulationType &get(size_t lineage_index, size_t query_index) {
    if (lineage_index >= _lineage_count || query_index >= _query_count) {
      LOG_ERROR("Index is too large");
    }
    return _table[lineage_index * _query_count + query_index];
  }

private:
  AccumulationType *_table;
  size_t            _lineage_count;
  size_t            _query_count;
};

class TaxaMask {
public:
  TaxaMask(size_t total_size)
      : _mask{(uint32_t *)calloc(compute_element_count(total_size),
                                 sizeof(uint32_t))},
        _size_in_bits(total_size) {}
  ~TaxaMask() { release_memory(); }

  explicit TaxaMask(const TaxaMask &other)
      : _mask{(uint32_t *)malloc(compute_element_count(other._size_in_bits))},
        _size_in_bits(other._size_in_bits) {
    memcpy(_mask, other._mask, size_in_bytes());
  }

  TaxaMask &operator=(const TaxaMask &) = delete;

  TaxaMask(TaxaMask &&rhs)
      : _mask{rhs._mask}, _size_in_bits{rhs._size_in_bits} {
    rhs._mask = nullptr;
  }

  TaxaMask &operator=(TaxaMask &&rhs) {
    std::swap(*this, rhs);
    rhs.release_memory();
    return *this;
  }

  uint32_t const *begin() const { return _mask; }
  uint32_t const *end() const { return _mask + size_in_elements(); }

  inline uint32_t operator[](size_t index) const { return _mask[index]; }

  uint32_t extract_tip_state(size_t index) const {
    size_t mask_index  = index / _bits_per_element;
    size_t mask_offset = index % _bits_per_element;

    return bextr(_mask[mask_index], mask_offset);
  }

  void set_bits(size_t start) {
    for (size_t i = start; i < _size_in_bits; ++i) {
      size_t mask_index  = i / _bits_per_element;
      size_t mask_offset = i % _bits_per_element;

      _mask[mask_index] |= 1u << i;
    }
  }

  constexpr size_t size_in_elements() const {
    return compute_element_count(_size_in_bits);
  }

  constexpr size_t size_in_bytes() const {
    return compute_element_count(_size_in_bits) * sizeof(uint32_t);
  }

  constexpr size_t size_in_bits() const { return _size_in_bits; }

private:
  void release_memory() {
    if (_mask) { free(_mask); }
    _mask = nullptr;
  }

  constexpr size_t compute_element_count(size_t total_size) const {
    size_t ele_count = total_size / _bits_per_element;
    if ((total_size % _bits_per_element) > 0) { ele_count += 1; }
    return ele_count;
  }

  static constexpr size_t _bits_per_element = sizeof(uint32_t) * 8;

  uint32_t *_mask;
  size_t    _size_in_bits;

  static_assert(sizeof(size_t) == sizeof(uint32_t *));
};

class TaxaList {
public:
  TaxaList(const std::vector<std::string> &labels) : _labels{labels} {}

  TaxaList(const YAML::Node &yaml) {
    for (auto n : yaml) { _labels.push_back(n.as<std::string>()); }
  }

  enum class find_error {
    not_found,
  };

  expected<size_t, find_error>
  find_label_index(const std::string_view &label) const {
    auto res = std::lower_bound(_labels.begin(), _labels.end(), label);
    if (res == _labels.end() || *res != label) {
      return unexpected{find_error::not_found};
    }
    return std::distance(_labels.begin(), res);
  }

  void sort() { std::sort(_labels.begin(), _labels.end()); }

  size_t size() const { return _labels.size(); }

  auto begin() { return _labels.begin(); }
  auto end() { return _labels.end(); }

  TaxaMask make_mask(size_t offset) const {
    TaxaMask lm(offset + size());

    lm.set_bits(offset);

    return lm;
  }

  std::string_view operator[](size_t index) const { return _labels[index]; }

private:
  std::vector<std::string> _labels;
};

class TreeList {
public:
  static TreeList parse_tree_file(const Path &treeset_file) {
    TreeList      tree_list;
    std::string   line;
    std::ifstream infile(treeset_file);
    while (std::getline(infile, line)) {
      auto t = corax_utree_parse_newick_string_unroot(line.c_str());
      if (t == nullptr) {
        LOG_ERROR("Failed to parse tree: '{}'", line);
        LOG_ERROR("msg: {}", corax_errmsg);
      }
      tree_list._trees.push_back(t);
    }

    // tree_list.make_tree_node_ids_consistent();

    return tree_list;
  }

  auto begin() const { return _trees.begin(); }
  auto end() const { return _trees.end(); }

  ~TreeList() {
    for (auto t : _trees) {
      if (t) { corax_utree_destroy(t, nullptr); }
    }
  }

  enum class find_error {
    not_found,
  };

  expected<uint32_t, find_error> get_node_id(const std::string_view &label) {
    auto &tree = _trees.front();
    for (size_t i = 0; i < tree->tip_count; ++i) {
      auto &node = tree->nodes[i];
      if (label == node->label) { return {node->node_index}; }
    }
    return unexpected{find_error::not_found};
  }

  void normalize_node_ids(const TaxaList &lineages, const TaxaList &queries) {
    size_t cur_idx = 0;

    size_t queries_offset = lineages.size();
    size_t all_offset     = lineages.size() + queries.size();

    auto &tree = _trees.front();
    for (size_t i = 0; i < tree->tip_count; ++i) {
      auto &node = tree->nodes[i];

      auto idx = lineages.find_label_index(node->label);
      if (idx.has_value()) {
        node->node_index = idx.value();
        continue;
      }

      idx = queries.find_label_index(node->label);
      if (idx.has_value()) {
        node->node_index = idx.value() + queries_offset;
        continue;
      }

      node->node_index  = cur_idx + all_offset;
      cur_idx          += 1;
    }
    make_tree_node_ids_consistent();
  }

private:
  void make_tree_node_ids_consistent() {
    auto &standard_tree = _trees.front();
    for (size_t i = 1; i < _trees.size(); ++i) {
      auto res = corax_utree_consistency_set(standard_tree, _trees[i]);
      if (res != CORAX_SUCCESS) {
        throw std::runtime_error{"Failed to set consistency"};
      }
    }
  }
  std::vector<Tree *> _trees;
};

struct Split {
  static constexpr size_t _bits_per = sizeof(uint32_t) * 8;

  uint32_t extract_tip_state(size_t index) const {
    size_t mask_index  = index / _bits_per;
    size_t mask_offset = index % _bits_per;

    return bextr(_split[mask_index], mask_offset);
  }

  inline uint32_t mask_and_popcount(const TaxaMask &mask) {
    uint32_t popcnt = 0;
    for (size_t i = 0; i < mask.size_in_elements(); ++i) {
      popcnt += popcount(_split[i] & mask[i]);
    }
    return popcnt;
  }

  inline void score(AccumulationTable &table,
                    const TaxaMask    &lineage_mask,
                    const TaxaMask    &query_mask) {
    auto count = mask_and_popcount(lineage_mask);

    if (!(count == 1 || count == lineage_mask.size_in_bits() - 1)) { return; }

    bool invert = count == 1 ? false : true;

    /* find the first element with a set bit */
    size_t set_element_index = 0;
    while (true) {
      auto res = _split[set_element_index] & lineage_mask[set_element_index];
      res      = invert ? ~res : res;
      if (res != 0) { break; }
      set_element_index++;
    }

    uint32_t labv = _split[set_element_index] & lineage_mask[set_element_index];
    labv          = invert ? ~labv & lineage_mask[set_element_index] : labv;

    size_t set_bit_offset = find_first_set(labv);

    size_t lineage_index = set_element_index * _bits_per + set_bit_offset;

    for (size_t i = lineage_mask.size_in_bits(); i < query_mask.size_in_bits();
         ++i) {
      auto tip_state = extract_tip_state(i);
      tip_state      = invert ? !tip_state : tip_state;
      table.get(lineage_index, i - lineage_mask.size_in_bits()) += tip_state;
    }
  }

  /* This is a _non owning_ pointer */
  uint32_t *_split;
};

class SplitSet {
public:
  ~SplitSet() { release_splits(); }

  SplitSet(Tree *const tree) {
    static_assert(sizeof(Split) == sizeof(corax_split_base_t *));
    _splits = reinterpret_cast<Split *>(
        corax_utree_split_create(tree->vroot, tree->tip_count, nullptr));

    _split_count = tree->edge_count - tree->tip_count;
    _split_len   = compute_element_count(tree->tip_count);
  }

  SplitSet(const SplitSet &)            = delete;
  SplitSet &operator=(const SplitSet &) = delete;

  SplitSet(SplitSet &&other)
      : _splits{other._splits},
        _split_count{other._split_count},
        _split_len{other._split_len} {
    other._splits      = nullptr;
    other._split_len   = 0;
    other._split_count = 0;
  }

  SplitSet &operator=(SplitSet &&other) {
    std::swap(*this, other);
    if (other._splits) { release_splits(); }
    other._splits      = nullptr;
    other._split_len   = 0;
    other._split_count = 0;
    return *this;
  }

  size_t split_count() const { return _split_count; }

  Split &operator[](size_t index) { return _splits[index]; }

  Split *begin() { return _splits; }
  Split *end() { return _splits + _split_count; }

  void accumulate(AccumulationTable &accumulation_table,
                  const TaxaMask    &lineage_mask,
                  const TaxaMask    &query_mask) {
    for (size_t i = 0; i < _split_count; ++i) {
      _splits->score(accumulation_table, lineage_mask, query_mask);
    }
  }

private:
  constexpr size_t compute_element_count(size_t total_size) {
    size_t ele_count = total_size / Split::_bits_per;
    if ((total_size % Split::_bits_per) > 0) { ele_count += 1; }
    return ele_count;
  }

  void release_splits() {
    if (_splits) {
      free(_splits[0]._split);
      free(_splits);
    }
  }

  Split *_splits;
  size_t _split_count;
  size_t _split_len;
};

class SplitSetList {
public:
  SplitSetList(const TreeList &trees) : _splits{trees.begin(), trees.end()} {}

  SplitSet &operator[](size_t index) { return _splits[index]; }
  size_t    size() const { return _splits.size(); }

  AccumulationTable accumulate(const TaxaList &lineage_list,
                               const TaxaList &query_list) {
    AccumulationTable accumulation_table(lineage_list.size(),
                                         query_list.size());

    auto lineage_mask = lineage_list.make_mask(0);
    auto query_mask   = query_list.make_mask(lineage_list.size());

    for (auto &split : _splits) {
      split.accumulate(accumulation_table, lineage_mask, query_mask);
    }

    return accumulation_table;
  }

  size_t total_splits() const { return size() * _splits.front().split_count(); }

private:
  std::vector<SplitSet> _splits;
};

struct ProgramOptions {
  Path treeset_file;
  Path output_prefix;
  Path yaml_config;
};

int main(int argc, char **argv) {
  logger::get_log_states().add_stream(
      stdout,
      logger::log_level::info | logger::log_level::warning
          | logger::log_level::important | logger::log_level::error
          | logger::log_level::progress);
  ProgramOptions options;
  CLI::App       app{"A project for chase :)"};

  app.add_option("--treeset", options.treeset_file);
  app.add_option("--config", options.yaml_config);
  // app.add_option("--prefix", options.output_prefix);

  CLI11_PARSE(app, argc, argv);

  LOG_INFO("Parsing trees");
  auto tree_list = TreeList::parse_tree_file(options.treeset_file);

  auto yaml = YAML::LoadFile(options.yaml_config);

  constexpr auto LINEAGE_KEY = "lineages";
  constexpr auto QUERIES_KEY = "queries";

  TaxaList lineage_list(yaml[LINEAGE_KEY]);
  TaxaList query_list(yaml[QUERIES_KEY]);

  LOG_INFO("Sorting taxa lists");

  lineage_list.sort();
  query_list.sort();

  LOG_INFO("Normalizing trees");
  tree_list.normalize_node_ids(lineage_list, query_list);

  LOG_INFO("Making splits");
  SplitSetList split_set_list{tree_list};
  auto        &ss = split_set_list[0];

  LOG_INFO("Computing results");
  auto table = split_set_list.accumulate(lineage_list, query_list);

  LOG_INFO("Total Splits: {}", split_set_list.total_splits());
  for (size_t i = 0; i < lineage_list.size(); ++i) {
    for (size_t j = 0; j < query_list.size(); ++j) {
      LOG_INFO("lineage: {}, query: {}, count: {}",
               lineage_list[i],
               query_list[j],
               table.get(i, j));
    }
  }
}
