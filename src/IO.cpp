#include "IO.hpp"

#include <array>
#include <fstream>
#include <numeric>
#include <string_view>

using namespace std::string_view_literals; // for the 'sv' suffix

template <typename T, size_t N>
inline std::string make_csv_row(const std::array<T, N> &entries) {
  return std::accumulate(std::next(entries.begin()),
                         entries.end(),
                         std::string(*entries.begin()),
                         [](const T &acc, const T &entry) -> std::string {
                           return acc + ", " + entry;
                         })
       + "\n";
}

template <size_t N>
inline std::string make_csv_row(const std::array<std::string_view, N> &fields) {
  return std::accumulate(
             std::next(fields.begin()),
             fields.end(),
             std::string(*fields.begin()),
             [](std::string acc, const std::string_view &entry) -> std::string {
               acc += ", ";
               acc += entry;
               return acc;
             })
       + "\n";
}

void write_results_to_csv(const AccumulationTable     &table,
                          const TaxaList              &lineage_list,
                          const TaxaList              &query_list,
                          const std::filesystem::path &output_filename) {
  std::ofstream        csv_file{output_filename};
  constexpr std::array fields{"lineage"sv, "query"sv, "matches"sv};
  csv_file << make_csv_row(fields);
  for (size_t i = 0; i < lineage_list.size(); ++i) {
    for (size_t j = 0; j < query_list.size(); ++j) {
      auto matches_string = std::to_string(table.get(i, j));
      csv_file << make_csv_row(std::array{
          lineage_list[i], query_list[j], std::string_view(matches_string)});
    }
  }
}
