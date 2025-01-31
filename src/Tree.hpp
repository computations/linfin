#pragma once
#include "Taxa.hpp"
#include "Util.hpp"

#include <corax/corax.hpp>
#include <fstream>
#include <filesystem>
#include <logger.hpp>

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
