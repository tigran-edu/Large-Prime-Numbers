#pragma once

#include "basic.hpp"

namespace lpn
{

int ComputeLegendreSymbol(long_int n, long_int p);
bool IsQuadraticResidue(const long_int & n, const long_int & p);
std::vector<size_t> FindQuadraticResiduePrimes(const long_int & n, size_t factor_size);

class QuadraticCongruences
{
   public:
    static long_int Solve(const long_int & n, const long_int & p);

   private:
    static long_int SolveCongruence(const long_int & n, const long_int & h, const long_int & p);
    static long_int FindStartValue(const long_int & n, const long_int & p);
};

};  // namespace lpn
