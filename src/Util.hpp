#pragma once
#include <corax/corax.hpp>
#include <cstdint>
#include <expected>
#include <filesystem>

using Path = std::filesystem::path;
using Tree = corax_utree_t;
using std::expected;
using std::unexpected;

inline uint32_t bextr(uint32_t a, std::size_t i) { return (a >> i) & 1ULL; }

inline uint32_t popcount(uint32_t a) { return __builtin_popcount(a); }

inline uint32_t find_first_set(uint32_t a) { return __builtin_ctz(a); }
