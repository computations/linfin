#pragma once

#include "Util.hpp"

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdlib.h>
#include <string>
#include <utility>
#include <vector>
#include <yaml-cpp/yaml.h>

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
