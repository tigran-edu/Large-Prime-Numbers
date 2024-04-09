#include "qs.hpp"
#include "congruence.hpp"

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
            congruences_.push_back(QuadraticCongruences::Solve(n, p));
        }
    }
}

FactorSet QuadraticSieveFactorization::Factorize(const long_int & n, const Sieve::Config & config)
{
    FactorSet factor;
    auto solution = Sieve::Solve(n, config);
    GaussianBasic gs = GaussianBasic(solution.factors, solution.primes);
    auto matrix = gs.Solve();
    return FindFactor(solution, matrix, n);
}

FactorSet QuadraticSieveFactorization::Factorize(const long_int & n) { return Factorize(n, Sieve::CreateConfig(n)); }

bool QuadraticSieveFactorization::IsPerfectSquare(const boost::dynamic_bitset<> & mask) { return mask.none(); }

std::vector<size_t> QuadraticSieveFactorization::GetParticipantsPositions(const Line & line)
{
    std::vector<size_t> positions;
    for (size_t i = 0; i < line.participants.size(); ++i)
    {
        if (line.participants[i])
        {
            positions.push_back(i);
        }
    }
    return positions;
}

long_int QuadraticSieveFactorization::ComputeX(const Sieve::Solution & solution, const std::vector<size_t> & positions,
                                               const long_int & n)
{
    long_int x = 1;
    for (auto pos : positions)
    {
        x = (x * (solution.r + solution.positions[pos])) % n;
    }
    return x;
}

long_int QuadraticSieveFactorization::ComputeY(const Sieve::Solution & solution, const std::vector<size_t> & positions,
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

FactorSet QuadraticSieveFactorization::FindFactor(const Sieve::Solution & solution, const Matrix & matrix,
                                                  const long_int & n)
{
    FactorSet factor;
    for (const auto & line : matrix)
    {
        if (IsPerfectSquare(line.mask))
        {
            auto positions = GetParticipantsPositions(line);
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

Sieve::Solution Sieve::Solve(const long_int & n, const Config & config)
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
