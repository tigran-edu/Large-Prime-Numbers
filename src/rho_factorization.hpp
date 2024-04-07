#pragma once

#include "basic.hpp"
#include <optional>

namespace lpn
{

class RhoFactorization
{
   private:
    struct BasicConfig
    {
        static const inline long_int kDefCValue = 3;
        static constexpr size_t kDefMaxIterValue = 1000000;
        static constexpr size_t kDefFreqValue = 10;
        static constexpr size_t kDefMaxAttempValue = 100;
    };

   public:
    RhoFactorization();
    RhoFactorization(const long_int & c, size_t max_iter, size_t frequency, size_t max_attemp);
    FactorSet Factorize(const long_int & n, size_t starting_point);

   private:
    size_t ComputeNextRange(size_t range);
    void UpdateX1(const long_int & x2, long_int & x1);
    void UpdateX2(const long_int & n, size_t range, long_int & x2);
    std::optional<long_int> TryToFindDivisor(const long_int & n, size_t terms, size_t range, long_int & x1,
                                             long_int & x2);
    long_int FindDivisor(const long_int & n, long_int & x1);
    long_int Next(const long_int & x2, const long_int & n);
    const long_int c_;
    const size_t max_iter_;
    const size_t frequency_;
    const size_t max_attemp_;
};

};  // namespace lpn
