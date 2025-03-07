#pragma once
#include "Accumulation.hpp"
#include "Taxa.hpp"
#include "Tree.hpp"
#include "Util.hpp"

#include <sstream>

struct Split {
  static constexpr size_t _bits_per = sizeof(uint32_t) * 8;

  uint32_t extract_tip_state(size_t index) const {
    size_t mask_index  = index / _bits_per;
    size_t mask_offset = index % _bits_per;

    return bextr(_split[mask_index], mask_offset);
  }

  std::string to_string(const TaxaList &linages,
                        const TaxaList &queries) const {
    std::ostringstream left_oss;
    std::ostringstream right_oss;

    for (size_t i = 0; i < linages.size(); ++i) {
      if (extract_tip_state(i)) {
        left_oss << linages[i];
      } else {
        right_oss << linages[i];
      }
    }

    for (size_t i = 0; i < queries.size(); ++i) {
      if (extract_tip_state(linages.size() + i)) {
        left_oss << queries[i];
      } else {
        right_oss << queries[i];
      }
    }

    left_oss << "|" << right_oss.str();
    return left_oss.str();
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

    auto lineage_index = find_first_set(labv);
    for (size_t j = lineage_mask.size_in_bits(); j < query_mask.size_in_bits();
         ++j) {
      auto tip_state = extract_tip_state(j);
      tip_state      = invert ? !tip_state : tip_state;
      table.get(lineage_index, j - lineage_mask.size_in_bits()) += tip_state;
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
      _splits[i].score(accumulation_table, lineage_mask, query_mask);
    }
  }

  void print(const TaxaList &lineages, const TaxaList &queries) const {
    for (size_t i = 0; i < _split_count; ++i) {
      LOG_INFO("{}", _splits[i].to_string(lineages, queries));
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
