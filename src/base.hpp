#pragma once

#include "gaussian.hpp"
#include <boost/multiprecision/cpp_int.hpp>

namespace lpn
{

struct Solution
{
    FactorSets factors;
    std::vector<long_int> values;
    std::vector<size_t> primes;
};

class FactorizationBase
{
   protected:
    using Matrix = GaussianBasic::Matrix;
    using Line = GaussianBasic::Line;

    static bool IsPerfectSquare(const Line & line);
    static std::vector<size_t> GetParticipantsPositions(const Line & line);
    static FactorSet FindFactor(const Solution & solution, const Matrix & matrix, const long_int & n);
    static long_int ComputeY(const Solution & solution, const std::vector<size_t> & positions, const long_int & n);
    static long_int ComputeX(const Solution & solution, const std::vector<size_t> & positions, const long_int & n);
};
};  // namespace lpn
