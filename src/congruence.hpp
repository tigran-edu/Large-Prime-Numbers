#pragma once

#include "basic.hpp"

namespace lpn
{

int ComputeLegendreSymbol(long_int n, long_int p);

class QuadraticCongruences
{
   public:
    static long_int Solve(const long_int & n, const long_int & p);

   private:
    static long_int SolveCongruence(const long_int & n, const long_int & h, const long_int & p);
    static long_int FindStartValue(const long_int & n, const long_int & p);
    static std::vector<bool> ToBinaryFormat(long_int val);
};

};  // namespace lpn
