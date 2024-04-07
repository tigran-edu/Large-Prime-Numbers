#pragma once

#include <optional>

#include "gaussian.hpp"

namespace lpn
{
struct SieveResult
{
    long_int r;
    std::vector<size_t> primes;
    std::vector<size_t> positions;
    FactorSets factors;
};

class QuadraticCongruences
{
   public:
    static long_int SolvingQuadraticCongruences(const long_int & n, const long_int & p);

   private:
    static long_int Solve(const long_int & n, const long_int & h, const long_int & p);
    static long_int FindStartValue(const long_int & n, const long_int & p);
    static std::vector<bool> Binary(long_int val);
};

class QuadraticSieve
{
   public:
    static FactorSet Factorize(const long_int & n);

   private:
    using Matrix = GaussianBasic::Matrix;
    using Line = GaussianBasic::Line;

   private:
    static long_int CheckResults(const SieveResult & result, const Line & solution, const long_int & n);
    static bool IsPerfectSquare(const boost::dynamic_bitset<> & mask);
    static FactorSet FindFactor(const SieveResult & result, const Matrix & solutions, const long_int & n);
    static std::vector<size_t> GetParticipantsPositions(const Line & solution);
    static long_int ComputeX(const SieveResult & result, const std::vector<size_t> & positions, const long_int & n);
    static long_int ComputeY(const SieveResult & result, const std::vector<size_t> & positions, const long_int & n);
};

struct Sieve
{
    static SieveResult Sieving(const long_int & n);

    struct Config
    {
        static Config CreateConfig(const long_int & n);

       private:
        void ComputeCloseness();

        void Reserve();

        void ComputeTarget(const long_int & n);

        void ComputeSmallPrimes(const long_int & n);

       public:
        float t{1.5};
        size_t m{50'000'000};
        size_t factor_size{2000};
        float target;
        float closeness;
        std::vector<size_t> primes;
        std::vector<long_int> congruences;
    };

   private:
    static void AddLogAtSpecificPositions(size_t i, size_t j, size_t p, std::vector<float> & logs);
    static std::optional<FactorSet> TryToDecompose(const Config & cf, size_t i, const long_int & n);
    static long_int ComputeTargetFunction(const long_int & r, const long_int & n, size_t i);
};

};  // namespace lpn
