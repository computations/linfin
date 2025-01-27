#include <CLI/App.hpp>
#include <CLI/Config.hpp>
#include <CLI/Formatter.hpp>
#include <corax/corax.hpp>
#include <filesystem>
#include <fstream>
#include <logger.hpp>
#include <stdexcept>
#include <string>
#include <vector>

using Path         = std::filesystem::path;
using Tree         = corax_utree_t;
using TreeList     = std::vector<Tree *>;
using Split        = corax_split_t;
using SplitList    = std::vector<Split>;
using SplitSet     = corax_split_set_t;
using SplitSetList = std::vector<corax_split_set_t>;

struct ProgramOptions {
  Path treeset_file;
  Path output_prefix;
};

TreeList parse_tree_file(const Path &treeset_file) {
  std::vector<Tree *> trees;
  std::string         line;
  std::ifstream       infile(treeset_file);
  while (std::getline(infile, line)) {
    auto t = corax_utree_parse_newick_string_unroot(line.c_str());
    if (t == nullptr) {
      LOG_ERROR("Failed to parse tree: '{}'", line);
      LOG_ERROR("msg: {}", corax_errmsg);
    }
    trees.push_back(t);
  }
  return trees;
}

void make_trees_consistent(const TreeList &trees) {
  auto &standard_tree = trees.front();
  for (size_t i = 1; i < trees.size(); ++i) {
    auto res = corax_utree_consistency_set(standard_tree, trees[i]);
    if (res != CORAX_SUCCESS) {
      throw std::runtime_error{"Failed to set consistency"};
    }
  }
}

SplitSet *parse_splits(Tree *tree) {
  auto ss = corax_utree_splitset_create(tree);

  return ss;
}

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

  auto trees = parse_tree_file(options.treeset_file);
  LOG_INFO("Found {} trees", trees.size());

  make_trees_consistent(trees);

  for (auto tree : trees) {
    auto ss = parse_splits(tree);
    for (size_t i = 0; i < ss->split_count; ++i) {
      LOG_INFO("Split {}: {:010b}", i, ss->splits[i][0]);
    }
    corax_utree_splitset_destroy(ss);
  }

  for (auto t : trees) { corax_utree_destroy(t, nullptr); }
}
