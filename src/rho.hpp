#pragma once

#include "aliases.hpp"
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
        static const inline std::array<long_int, 7> kPrimes = {2, 3, 5, 7, 11, 13};
        static constexpr size_t kDefMaxIterValue = 1000000;
        static constexpr size_t kDefFreqValue = 10;
        static constexpr size_t kDefMaxAttempValue = 100;
    };

   public:
    static FactorSet Factorize(const long_int & n, size_t starting_point);
    static FactorSet Factorize(const long_int & n, const long_int & c, size_t max_iter, size_t frequency,
                               size_t max_attemp, size_t starting_point);

   private:
    RhoFactorization();
    RhoFactorization(const long_int & c, size_t max_iter, size_t frequency, size_t max_attemp);

   private:
    FactorSet FindFactor(const long_int & n, size_t starting_point);
    size_t ComputeNextRange(size_t range) const;
    void UpdateX1(const long_int & x2, long_int & x1);
    void UpdateX2(const long_int & n, size_t range, long_int & x2);
    std::optional<long_int> TryToFindDivisor(const long_int & n, size_t terms, size_t range, long_int & x1,
                                             long_int & x2);

    long_int FindDivisor(const long_int & n, long_int & x1);
    long_int Next(const long_int & x2, const long_int & n) const;

    const long_int c_;
    const size_t max_iter_;
    const size_t frequency_;
    const size_t max_attemp_;
};

};  // namespace lpn
