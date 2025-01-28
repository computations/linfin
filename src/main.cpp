#include <CLI/App.hpp>
#include <CLI/Config.hpp>
#include <CLI/Formatter.hpp>
#include <algorithm>
#include <cmath>
#include <corax/corax.hpp>
#include <cstdint>
#include <expected>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <logger.hpp>
#include <string>
#include <string_view>
#include <vector>

using Path = std::filesystem::path;
using Tree = corax_utree_t;
using std::expected;
using std::unexpected;

inline uint32_t bextr(uint32_t a, size_t i) { return (a >> i) & 1ULL; }

class LinageList {
public:
  LinageList(const std::vector<std::string> &labels) : _labels{labels} {
    _counts.resize(_labels.size());
  }

  enum class find_error {
    not_found,
  };

  expected<size_t, find_error>
  find_label_index(const std::string_view &label) const {
    auto res = std::lower_bound(_labels.begin(), _labels.end(), label);
    if (res == _labels.end()) { return unexpected{find_error::not_found}; }
    return std::distance(_labels.begin(), res);
  }

  void sort() { std::sort(_labels.begin(), _labels.end()); }

  size_t size() const { return _labels.size(); }

  auto begin() { return _labels.begin(); }
  auto end() { return _labels.end(); }

private:
  std::vector<std::string> _labels;
  std::vector<size_t>      _counts;
};

class LineageMask {
public:
  ~LineageMask() { release_memory(); }

  LineageMask(const LineageMask &)            = delete;
  LineageMask &operator=(const LineageMask &) = delete;

  LineageMask(LineageMask &&rhs) : _mask{rhs._mask}, _size{rhs._size} {
    rhs._mask = nullptr;
  }

  LineageMask &operator=(LineageMask &&rhs) {
    std::swap(*this, rhs);
    rhs.release_memory();
    return *this;
  }

  uint32_t const *begin() const { return _mask; }
  uint32_t const *end() const { return _mask + _size; }

private:
  void release_memory() {
    if (_mask) { free(_mask); }
    _mask = nullptr;
  }

  static constexpr size_t _bits_per = sizeof(uint32_t) * 8;
  uint32_t               *_mask;
  size_t                  _size;
  static_assert(sizeof(size_t) == sizeof(uint32_t *));
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

    tree_list.make_tree_node_ids_consistent();

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

  void normalize_node_ids(const LinageList &lineages) {
    size_t cur_idx = 0;
    auto  &tree    = _trees.front();
    for (size_t i = 0; i < tree->tip_count; ++i) {
      auto &node = tree->nodes[i];
      auto  idx  = lineages.find_label_index(node->label);
      if (idx.has_value()) {
        node->node_index = idx.value();
      } else {
        node->node_index  = lineages.size() + cur_idx;
        cur_idx          += 1;
      }
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
    _split_len   = tree->tip_count / Split::_bits_per;
    if ((tree->tip_count % Split::_bits_per) > 0) { _split_len += 1; }
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

private:
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

private:
  std::vector<SplitSet> _splits;
};

struct ProgramOptions {
  Path treeset_file;
  Path output_prefix;
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
  // app.add_option("--prefix", options.output_prefix);

  CLI11_PARSE(app, argc, argv);

  auto       tree_list = TreeList::parse_tree_file(options.treeset_file);
  LinageList lineage_list{{"a", "b", "c"}};
  tree_list.normalize_node_ids(lineage_list);
  SplitSetList split_set_list{tree_list};

  auto &ss = split_set_list[0];

  for (auto &split : ss) {
    LOG_INFO("split: {:010b}, a: {:b}, b: {:b}",
             *split._split,
             split.extract_tip_state(0),
             split.extract_tip_state(1));
  }
}
