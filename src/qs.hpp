#pragma once

#include <atomic>
#include <optional>

#include "gaussian.hpp"
#include "base.hpp"

namespace lpn
{

class Sieve
{
   private:
    struct BasicConfig
    {
        static constexpr size_t kDefSegmentSize = 50'000'000;
        static constexpr size_t kDefFactorSize = 2000;
        static constexpr float kExpansionRate = 1.5;
    };

   public:
    struct Config
    {
        friend Sieve;

       private:
        Config(size_t segment_size, size_t factor_size);

        void ComputeCloseness(float expansion);

        void Reserve();

        void ComputeTarget(const long_int & n);

        void ComputeSmallPrimes(const long_int & n);

        size_t segment_size_;
        size_t factor_size_;
        float target_;
        float closeness_;
        std::vector<size_t> primes_;
        std::vector<long_int> congruences_;
    };

   private:
    Sieve(const Config & config, const long_int & n);

   public:
    static Solution Solve(const long_int & n, const Config & config);
    static Solution SolveMultiThread(const long_int & n, const Config & config);
    static Config CreateConfig(const long_int & n);
    static Config CreateConfig(const long_int & n, size_t segment_size, size_t factor_size, float expansion_rate);

   private:
    void ComputeSieve(const Config & config, size_t left_border, size_t right_border);
    void ComputeSieveMultiThread(const Config & config, const long_int & n);
    Solution FindAllFactorizable(const Config & config, const long_int & n);
    void AddFactor(const std::optional<FactorSet> & factor, size_t i, Solution & solution) const;
    void AddPrimeInSieve(size_t i, size_t j, size_t p, size_t left_border, size_t right_border);
    static long_int ComputeTargetFunction(const long_int & r, const long_int & n, size_t i);
    size_t ComputeLeftBorder(size_t threads, size_t position) const;
    size_t ComputeRightBorder(size_t threads, size_t position) const;

   private:
    const long_int r_;
    std::vector<float> sieve_;
};

class QuadraticSieveFactorization : private FactorizationBase
{
   public:
    static FactorSet Factorize(const long_int & n, bool multi_thread = false);
    static FactorSet Factorize(const long_int & n, const Sieve::Config & config, bool multi_thread = false);
};

};  // namespace lpn
