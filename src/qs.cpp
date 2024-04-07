#include "qs.hpp"

#include <climits>
#include <cstdint>
#include <random>

namespace lpn
{
namespace math = boost::multiprecision;
using Config = Sieve::Config;

Config::Config(size_t segment_size, size_t factor_size) : segment_size_(segment_size), factor_size_(factor_size) {}

Config Sieve::CreateConfig(const long_int & n, size_t segment_size, size_t factor_size, float expansion_rate)
{
    Config config(segment_size, factor_size);
    config.Reserve();
    config.ComputeSmallPrimes(n);
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

void Config::ComputeSmallPrimes(const long_int & n)
{
    primes_.push_back(2);
    congruences_.push_back(n % 2);

    for (size_t p = 3; primes_.size() <= factor_size_; p += 2)
    {
        if (IsPrimeBasic(p) && ComputeLegendreSymbol(n, p) == 1)
        {
            primes_.push_back(p);
            congruences_.push_back(QuadraticCongruences::SolvingQuadraticCongruences(n, p));
        }
    }
}

FactorSet QuadraticSieve::Factorize(const long_int & n, const Sieve::Config & config)
{
    FactorSet factor;
    auto solution = Sieve::Sieving(n, config);
    GaussianBasic gs = GaussianBasic(solution.factors, solution.primes);
    auto solutions = gs.Solve();
    return FindFactor(solution, solutions, n);
}

FactorSet QuadraticSieve::Factorize(const long_int & n) { return Factorize(n, Sieve::CreateConfig(n)); }

bool QuadraticSieve::IsPerfectSquare(const boost::dynamic_bitset<> & mask) { return mask.none(); }

std::vector<size_t> QuadraticSieve::GetParticipantsPositions(const Line & solution)
{
    std::vector<size_t> positions;
    for (size_t i = 0; i < solution.participants.size(); ++i)
    {
        if (solution.participants[i])
        {
            positions.push_back(i);
        }
    }
    return positions;
}

long_int QuadraticSieve::ComputeX(const Sieve::Solution & solution, const std::vector<size_t> & positions,
                                  const long_int & n)
{
    long_int x = 1;
    for (auto pos : positions)
    {
        x = (x * (solution.r + solution.positions[pos])) % n;
    }
    return x;
}

long_int QuadraticSieve::ComputeY(const Sieve::Solution & solution, const std::vector<size_t> & positions,
                                  const long_int & n)
{
    long_int y = 1;
    FactorSet factor;
    for (auto pos : positions)
    {
        MergeFactorSets(factor, solution.factors[pos]);
    }
    for (auto [key, value] : factor)
    {
        assert(value % 2 == 0);
        y = (y * FastExponentiationWithMod(key, value / 2, n)) % n;
    }
    return y;
}

FactorSet QuadraticSieve::FindFactor(const Sieve::Solution & solution, const Matrix & solutions, const long_int & n)
{
    FactorSet factor;
    for (const auto & sol : solutions)
    {
        if (IsPerfectSquare(sol.mask))
        {
            auto positions = GetParticipantsPositions(sol);
            long_int x = ComputeX(solution, positions, n);
            long_int y = ComputeY(solution, positions, n);
            long_int gcd = lpn::gcd<long_int>(math::abs(x - y), n);
            if (gcd != 1 && gcd != n)
            {
                factor[gcd] = 1;
                factor[n / gcd] = 1;
                return factor;
            }
        }
    }
    factor[n] = 1;
    return factor;
}

long_int QuadraticCongruences::SolvingQuadraticCongruences(const long_int & n, const long_int & p)
{
    assert(p % 2 != 0);

    if (p % 4 == 3)
    {
        return FastExponentiationWithMod(n, (p + 1) / 4, p);
    }
    if (p % 8 == 5 && FastExponentiationWithMod(n, (p - 1) / 4, p) == 1)
    {
        return FastExponentiationWithMod(n, (p + 3) / 8, p);
    }
    if (p % 8 == 5 && FastExponentiationWithMod(n, (p - 1) / 4, p) == p - 1)
    {
        return (FastExponentiationWithMod(4 * n, (p + 3) / 8, p) * (p + 1) / 2) % p;
    }
    long_int h = FindStartValue(n, p);
    return (Solve(n, h, p) * (p + 1) / 2) % p;
}

long_int QuadraticCongruences::FindStartValue(const long_int & n, const long_int & p)
{
    static std::mt19937 gen{42};
    long_int h = 1;
    std::uniform_int_distribution<size_t> distrib(1, UINT64_MAX);
    while (ComputeLegendreSymbol(h * h - 4 * n, p) != -1)
    {
        h = distrib(gen);
    }
    return h;
}

long_int QuadraticCongruences::Solve(const long_int & n, const long_int & h, const long_int & p)
{
    long_int m = n;
    long_int v = h;
    long_int w = (h * h - 2 * n) % p;
    auto bitset = ToBinaryFormat((p + 1) / 2);

    for (size_t i = 1; i < bitset.size(); ++i)
    {
        long_int x = (v * w - h * m) % p;
        v = (v * v - 2 * m) % p;
        w = (w * w - 2 * n * m) % p;
        m = (m * m) % p;
        if (bitset[i])
        {
            w = x;
        }
        else
        {
            v = x;
            m = (n * m) % p;
        }
    }
    return (p + v) % p;
}

std::vector<bool> QuadraticCongruences::ToBinaryFormat(long_int val)
{
    std::vector<bool> bitset;
    while (val > 0)
    {
        bitset.push_back(val % 2 == 0);
        val >>= 1;
    }
    std::reverse(bitset.begin(), bitset.end());
    return bitset;
}

Sieve::Solution Sieve::Sieving(const long_int & n, const Config & config)
{
    Sieve::Solution solution;
    std::vector<float> logs(2 * config.segment_size_, 0.0);
    const long_int r = math::sqrt(n) - config.segment_size_;

    for (size_t pos = 0; pos < config.primes_.size(); ++pos)
    {
        size_t p = config.primes_[pos];
        // n = x^2 mod p; f(r) = r^2 - n; target f(r) == 0 mod p
        size_t i = (size_t)((config.congruences_[pos] + p - (r % p)) % p);
        size_t j = (size_t)((2 * p - config.congruences_[pos] - (r % p)) % p);
        AddLogAtSpecificPositions(i, j, p, logs);
    }
    for (size_t i = 0; i < logs.size(); ++i)
    {
        if (config.closeness_ <= logs[i] && logs[i] <= config.target_ &&
            solution.positions.size() < 1.15 * config.factor_size_)
        {
            auto factor = TryToDecompose(config, i, n);
            if (factor != std::nullopt)
            {
                solution.positions.push_back(i);
                solution.factors.push_back(std::move(factor.value()));
            }
        }
    }

    solution.primes = std::move(config.primes_);
    solution.r = r;
    return solution;
}

void Sieve::AddLogAtSpecificPositions(size_t i, size_t j, size_t p, std::vector<float> & logs)
{
    if (j < i)
    {
        std::swap(i, j);
    }

    float p_log = log10(float(p));

    if (i == j)
    {
        while (i < logs.size())
        {
            logs[i] += p_log;
            i += p;
        }
        return;
    }

    while (j < logs.size())
    {
        logs[i] += p_log;
        logs[j] += p_log;
        i += p;
        j += p;
    }

    if (i < logs.size())
    {
        logs[i] += p_log;
    }
}

std::optional<FactorSet> Sieve::TryToDecompose(const Config & config, size_t i, const long_int & n)
{
    long_int r = math::sqrt(n) - config.segment_size_;
    long_int value = ComputeTargetFunction(r, n, i);
    FactorSet factor;
    for (const auto p : config.primes_)
    {
        if (value % p == 0)
        {
            factor[p] += ExtractPowerFast(value, p);
        }
    }
    return value == 1 ? std::optional<FactorSet>(factor) : std::nullopt;
}

long_int Sieve::ComputeTargetFunction(const long_int & r, const long_int & n, size_t i)
{
    return math::abs((r + i) * (r + i) - n);
}

};  // namespace lpn
