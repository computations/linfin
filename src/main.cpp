#include "Split.hpp"
#include "Taxa.hpp"
#include "Tree.hpp"

#include <CLI/App.hpp>
#include <CLI/Config.hpp>
#include <CLI/Formatter.hpp>
#include <cmath>
#include <corax/corax.hpp>
#include <cstring>
#include <expected>
#include <logger.hpp>
#include <string>
#include <yaml-cpp/yaml.h>

struct ProgramOptions {
  Path treeset_file;
  Path output_prefix;
  Path yaml_config;
};

int main(int argc, char **argv) {
  logger::get_log_states().add_stream(
      stdout,
      logger::log_level::info | logger::log_level::warning |
          logger::log_level::important | logger::log_level::error |
          logger::log_level::progress);
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
