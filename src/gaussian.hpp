#pragma once

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
        Line() = default;
        Line(size_t n, size_t m);

        Line & operator^=(Line & line);

        boost::dynamic_bitset<> participants;
        boost::dynamic_bitset<> mask;
    };

    using Matrix = std::vector<Line>;

    static Matrix CreateMatrix(const FactorSets & factors, const std::vector<size_t> & primes);

    GaussianBasic(const FactorSets & factors, const std::vector<size_t> & primes);

    Matrix Solve();

    size_t FindFirstNonZeroInLine(size_t line_pos) const;

    void Add(size_t col, size_t line);

   private:
    size_t m_;
    size_t n_;
    Matrix matrix_;
};
};  // namespace lpn
