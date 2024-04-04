#pragma once

#include <vector>
#include <boost/dynamic_bitset.hpp>
#include "basic.hpp"

namespace lpn
{
class GaussianBasic
{
   public:
    class Line
    {
       public:
        Line(size_t n, size_t m) : participants(n), mask(m) {}

        Line & operator^=(Line & line)
        {
            mask ^= line.mask;
            participants ^= line.participants;
            return *this;
        }

        boost::dynamic_bitset<> participants;
        boost::dynamic_bitset<> mask;
    };

    using Matrix = std::vector<Line>;

    static Matrix CreateMatrix(const FactorSets & factors, const std::vector<size_t> & primes);

    GaussianBasic(const FactorSets & factors, const std::vector<size_t> & primes);

    Matrix Solve();

    size_t FindFirstNonZeroInLine(size_t line_pos);

    void Add(size_t col, size_t line);

   private:
    size_t m_;
    size_t n_;
    Matrix matrix_;
    const FactorSets & factors_;
};
};  // namespace lpn
