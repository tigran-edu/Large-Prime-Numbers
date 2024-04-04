#pragma once

#include <vector>
#include <boost/dynamic_bitset.hpp>
#include "basic.hpp"

namespace lpn
{
class GaussianBasic
{
   public:
    struct Bitset
    {
        Bitset(size_t n, size_t m) : participants(n), mask(m) {}

        void operator^=(Bitset & bitset)
        {
            mask ^= bitset.mask;
            participants ^= bitset.participants;
        }

        boost::dynamic_bitset<> participants;
        boost::dynamic_bitset<> mask;
    };

    GaussianBasic(const std::vector<FactorSet> & factors, const std::vector<size_t> & primes);

    std::vector<Bitset> Solve();

    void Add(size_t pos, size_t line);

    void Print();

   private:
    size_t m_;
    size_t n_;
    std::vector<FactorSet> factors_;
    std::vector<Bitset> matrix_;
};
};  // namespace lpn
