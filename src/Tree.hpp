#pragma once
#include "Taxa.hpp"
#include "Util.hpp"

#include <algorithm>
#include <corax/corax.hpp>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <logger.hpp>
#include <stdexcept>

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

      node->node_index = cur_idx + all_offset;
      cur_idx += 1;
    }
    make_tree_node_ids_consistent();
  }

private:
  struct label_index_pair {
    char  *label;
    size_t index;
  };

  std::vector<label_index_pair> fill_label_table(Tree *tree) {
    std::vector<label_index_pair> label_table;
    label_table.reserve(tree->tip_count);
    // char **label_table = (char **)malloc(tree->tip_count * sizeof(char *));
    for (size_t i = 0; i < tree->tip_count; ++i) {
      auto node = tree->nodes[i];
      label_table.push_back({node->label, node->node_index});
    }

    struct {
      bool operator()(const label_index_pair &a, const label_index_pair &b) {
        return strcmp(a.label, b.label) < 0;
      }
    } comp;

    std::sort(label_table.begin(), label_table.end(), comp);

    return label_table;
  }

  void set_node_ids_by_label(Tree                                *tree,
                             const std::vector<label_index_pair> &label_table) {
    for (size_t i = 0; i < tree->tip_count; ++i) {
      auto node = tree->nodes[i];

      auto res = std::lower_bound(label_table.begin(),
                                  label_table.end(),
                                  node->label,
                                  [](const label_index_pair &a, char *key) {
                                    return strcmp(a.label, key) < 0;
                                  });

      if (res == label_table.end()) {
        throw std::runtime_error{"Failed to set node id by label"};
      }

      node->node_index = res->index;
    }
  }

  void make_tree_node_ids_consistent() {
    auto &standard_tree = _trees.front();

    auto label_table = fill_label_table(standard_tree);

    for (size_t i = 1; i < _trees.size(); ++i) {
      set_node_ids_by_label(_trees[i], label_table);
    }
  }
  std::vector<Tree *> _trees;
};
