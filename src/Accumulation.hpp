#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <logger.hpp>
#include <stdlib.h>
#include <utility>

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
