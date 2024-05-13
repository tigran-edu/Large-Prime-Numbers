#include "qs.hpp"
#include "base.hpp"
#include "congruence.hpp"
#include "gaussian.hpp"
#include "thread_group.hpp"

namespace lpn
{
namespace math = boost::multiprecision;
using Config = Sieve::Config;

Config::Config(size_t segment_size, size_t factor_size, bool multi_thread)
    : segment_size_(segment_size), factor_size_(factor_size), multi_thread_(multi_thread)
{
}

Config Sieve::CreateConfig(const long_int & n, size_t segment_size, size_t factor_size, float expansion_rate,
                           bool multi_thread)
{
    Config config(segment_size, factor_size, multi_thread);
    config.Reserve();
    config.ComputePrimes(n);
    config.ComputeTarget(n);
    config.ComputeCloseness(expansion_rate);
    return config;
}

Config Sieve::CreateConfig(const long_int & n)
{
    return CreateConfig(n, BasicConfig::kDefSegmentSize, BasicConfig::kDefFactorSize, BasicConfig::kExpansionRate);
}

void Config::ComputeCloseness(float expansion_rate)
{
    closeness_ = target_ - expansion_rate * log10(float(primes_.back()));
}

void Config::Reserve()
{
    primes_.reserve(factor_size_);
    congruences_.reserve(factor_size_);
}

void Config::ComputeTarget(const long_int & n) { target_ = log10(double(n)) / 2 + log10(double(segment_size_)); }

void Config::ComputePrimes(const long_int & n)
{
    primes_.push_back(2);
    congruences_.push_back(n % 2);

    for (size_t p = 3; primes_.size() <= factor_size_; p += 2)
    {
        if (IsPrimeBasic(p) && ComputeLegendreSymbol(n, p) == 1)
        {
            primes_.push_back(p);
            congruences_.push_back(QuadraticCongruences::Solve(n, p));
        }
    }
}

FactorSet QuadraticSieveFactorization::Factorize(const long_int & n, const Sieve::Config & config)
{
    Solution solution = Sieve::Solve(n, config);
    GaussianBasic gs = GaussianBasic(solution.factors, solution.primes);
    Matrix matrix = gs.Solve();
    return FindFactor(solution, matrix, n);
}

FactorSet QuadraticSieveFactorization::Factorize(const long_int & n) { return Factorize(n, Sieve::CreateConfig(n)); }

Sieve::Sieve(const long_int & n, const Config & config)
    : r_(math::sqrt(n) - config.segment_size_), sieve_(2 * config.segment_size_, 0.0)
{
}

Solution Sieve::Solve(const long_int & n, const Config & config)
{
    Sieve sieve(n, config);
    if (config.multi_thread_)
    {
        sieve.ComputeSieveMultiThread(config, n);
    }
    else
    {
        sieve.ComputeSieve(config, 0, sieve.sieve_.size());
    }
    Solution solution = sieve.FindAllFactorizable(config, n);
    solution.primes = config.primes_;
    return solution;
}

void Sieve::ComputeSieve(const Config & config, size_t left_border, size_t right_border)
{
    for (size_t pos = 0; pos < config.primes_.size(); pos++)
    {
        size_t p = config.primes_[pos];
        // n = x^2 mod p; f(r) = r^2 - n; target f(r) == 0 mod p
        size_t i = (size_t)((config.congruences_[pos] + p - (r_ % p)) % p);
        size_t j = (size_t)((2 * p - config.congruences_[pos] - (r_ % p)) % p);
        AddPrimeInSieve(i, j, p, left_border, right_border);
    }
}

void Sieve::ComputeSieveMultiThread(const Config & config, const long_int & n)
{
    int threads = ThreadGroup::GetThreadAmount();
    auto group = ThreadGroup();
    for (size_t pos = 0; pos < threads; ++pos)
    {
        group.AddTask(
            [pos, config, this, threads]() mutable
            { this->ComputeSieve(config, ComputeLeftBorder(threads, pos), ComputeRightBorder(threads, pos)); });
    }
    group.ComputeAllTasks();
}

Solution Sieve::FindAllFactorizable(const Config & config, const long_int & n)
{
    Solution solution;
    for (size_t i = 0; i < sieve_.size(); ++i)
    {
        if (config.closeness_ <= sieve_[i] && solution.values.size() < 1.1 * config.factor_size_)
        {
            auto factor = TryToDecompose(config.primes_, ComputeTargetFunction(n, i));
            AddFactor(factor, i, solution);
        }
    }
    return solution;
}

void Sieve::AddFactor(const std::optional<FactorSet> & factor, size_t i, Solution & solution) const
{
    if (factor)
    {
        solution.values.push_back(i + r_);
        solution.factors.push_back(std::move(factor.value()));
    }
}

void Sieve::AddPrimeInSieve(size_t i, size_t j, size_t p, size_t left_border, size_t right_border)
{
    if (i < left_border)
    {
        i += ((left_border - i + p - 1) / p) * p;
    }
    if (j < left_border)
    {
        j += ((left_border - j + p - 1) / p) * p;
    }

    if (j < i)
    {
        std::swap(i, j);
    }

    float p_log = log10(float(p));

    if (i == j)
    {
        while (i < right_border)
        {
            sieve_[i] += p_log;
            i += p;
        }
        return;
    }

    while (j < right_border)
    {
        sieve_[i] += p_log;
        sieve_[j] += p_log;
        i += p;
        j += p;
    }

    if (i < right_border)
    {
        sieve_[i] += p_log;
    }
}

long_int Sieve::ComputeTargetFunction(const long_int & n, size_t i) const { return math::abs((r_ + i) * (r_ + i) - n); }

size_t Sieve::ComputeLeftBorder(size_t threads, size_t position) const
{
    return ((sieve_.size()) / threads) * position;
}

size_t Sieve::ComputeRightBorder(size_t threads, size_t position) const
{
    if (position + 1 == threads)
    {
        return sieve_.size();
    }
    return ((sieve_.size()) / threads) * (position + 1);
}

};  // namespace lpn
