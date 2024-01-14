#pragma once

#include <map>
#include "basic.h"

namespace lpn
{

using boost::multiprecision::cpp_int;

std::unordered_map<cpp_int, size_t> RhoFactorization(cpp_int n, cpp_int c, cpp_int max);

};