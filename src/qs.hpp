#pragma once

#include <map>
#include <boost/multiprecision/cpp_bin_float.hpp>

#include "basic.hpp"

namespace lpn
{
using boost::multiprecision::cpp_bin_float_100;

class QuadraticSieve
{
   public:
    static Factor Factorize(long_int n);
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

struct Sieve
{
    struct Config
    {
        explicit Config(const long_int & n)
        {
            target = log10(cpp_bin_float_100(n)) / 2 + log10(cpp_bin_float_100(m));
            ComputeSmallPrimes(n);
        }

        void ComputeCloseness(const long_int & max_prime)
        {
            closeness = target - t * log10(cpp_bin_float_100(max_prime));
        }

        void ComputeSmallPrimes(const long_int & n)
        {
            size_t p = 3;
            primes.push_back(2);
            congruences.push_back(n % 2);
            while (primes.size() <= factor_size)
            {
                if (BasicIsPrime(p) && LegendreSymbol::Compute(n, p) == 1)
                {
                    primes.push_back(p);
                    congruences.push_back(QuadraticCongruences::SolvingQuadraticCongruences(n, p));
                }
                p++;
            }
            ComputeCloseness(primes.back());
        }

        static size_t Size(const long_int & n) { return n.str().size(); }

        cpp_bin_float_100 t{1.5};
        size_t m{5'000};
        size_t factor_size{28};
        cpp_bin_float_100 target;
        cpp_bin_float_100 closeness;
        std::vector<size_t> primes;
        std::vector<long_int> congruences;
    };

    static std::vector<long_int> Sieving(long_int n);
    static void Fill(size_t i, size_t p, std::vector<cpp_bin_float_100> & logs);
    static bool IsDecomposed(Config & cf, size_t i, const long_int & n);
    static long_int Eval(const long_int & r, const long_int & n, size_t i);
};

};  // namespace lpn
