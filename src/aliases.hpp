#pragma once

#include <map>
#include <vector>

#include "boost/multiprecision/cpp_int.hpp"

namespace lpn
{
using long_int = boost::multiprecision::cpp_int;  // NOLINT
using FactorSet = std::unordered_map<long_int, size_t>;
using FactorSets = std::vector<FactorSet>;
};  // namespace lpn
