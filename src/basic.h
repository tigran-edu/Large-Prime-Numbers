#pragma once

#include <map>
#include "boost/multiprecision/cpp_int.hpp"

namespace lpn
{
using boost::multiprecision::cpp_int;

cpp_int gcd(cpp_int a, cpp_int b);
cpp_int FastExponentiation(cpp_int a, cpp_int b);
cpp_int FastExponentiationWithMod(cpp_int a, cpp_int b, cpp_int m);
bool PseudoprimeTest(cpp_int p);
std::unordered_map<cpp_int, size_t> BasicFactorization(cpp_int a);

};
