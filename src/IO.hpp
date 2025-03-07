#pragma once

#include "Accumulation.hpp"
#include "Taxa.hpp"

#include <filesystem>

void write_results_to_csv(const AccumulationTable     &table,
                          const TaxaList              &lineage_list,
                          const TaxaList              &query_list,
                          const std::filesystem::path &output_filename);
