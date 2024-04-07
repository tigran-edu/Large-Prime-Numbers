#pragma once

#include <optional>

#include "gaussian.hpp"

namespace lpn
{

class QuadraticCongruences
{
   public:
    static long_int SolvingQuadraticCongruences(const long_int & n, const long_int & p);

   private:
    static long_int Solve(const long_int & n, const long_int & h, const long_int & p);
    static long_int FindStartValue(const long_int & n, const long_int & p);
    static std::vector<bool> ToBinaryFormat(long_int val);
};

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

    struct Solution
    {
        long_int r;
        std::vector<size_t> primes;
        std::vector<size_t> positions;
        FactorSets factors;
    };

    static Solution Sieving(const long_int & n, const Config & config);
    static Config CreateConfig(const long_int & n);
    static Config CreateConfig(const long_int & n, size_t segment_size, size_t factor_size, float expansion_rate);

   private:
    static void AddLogAtSpecificPositions(size_t i, size_t j, size_t p, std::vector<float> & logs);
    static std::optional<FactorSet> TryToDecompose(const Config & config, size_t i, const long_int & n);
    static long_int ComputeTargetFunction(const long_int & r, const long_int & n, size_t i);
};

class QuadraticSieve
{
   public:
    static FactorSet Factorize(const long_int & n);
    static FactorSet Factorize(const long_int & n, const Sieve::Config & config);

   private:
    using Matrix = GaussianBasic::Matrix;
    using Line = GaussianBasic::Line;

   private:
    static bool IsPerfectSquare(const boost::dynamic_bitset<> & mask);
    static FactorSet FindFactor(const Sieve::Solution & sieve_solution, const Matrix & solutions, const long_int & n);
    static std::vector<size_t> GetParticipantsPositions(const Line & solution);
    static long_int ComputeX(const Sieve::Solution & solution, const std::vector<size_t> & positions,
                             const long_int & n);
    static long_int ComputeY(const Sieve::Solution & solution, const std::vector<size_t> & positions,
                             const long_int & n);
};

};  // namespace lpn
