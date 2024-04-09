#include "cfrac.hpp"
#include "aliases.hpp"
#include "basic.hpp"
#include "congruence.hpp"
#include "gaussian.hpp"

namespace lpn
{

namespace math = boost::multiprecision;

ContinuedFractions::State::State(const long_int & n)
    : n(n), sqrt_n(math::sqrt(n)), a(0), b0(0), b1(sqrt_n), c0(1), c1(n - sqrt_n * sqrt_n), p0(1), p1(sqrt_n)
{
}

void ContinuedFractions::State::Update()
{
    a = (sqrt_n + b1) / c1;
    b0 = std::exchange<long_int, long_int>(b1, a * c1 - b1);
    c0 = std::exchange<long_int, long_int>(c1, c0 + a * (b0 - b1));
    p0 = std::exchange<long_int, long_int>(p1, (p0 + a * p1) % n);
}

ContinuedFractions::Solution ContinuedFractions::Solve(const long_int & n, size_t factor_size)
{
    State state(n);
    Solution solution;
    solution.primes = FindQuadraticResiduePrimes(n, factor_size);
    while (solution.factors.size() < 1.1 * solution.primes.size())
    {
        state.Update();
        auto factor = TryToDecompose(solution.primes, state.c1);
        if (factor.has_value())
        {
            solution.factors.push_back(std::move(factor.value()));
            solution.xs.push_back(state.p1);
        }
    }
    return solution;
}

FactorSet ContinuedFractionsFactorization::Factorize(const long_int & n) { return Factorize(n, kBasicFactorSize); }

FactorSet ContinuedFractionsFactorization::Factorize(const long_int & n, size_t factor_size)
{
    auto solution = ContinuedFractions::Solve(n, factor_size);
    GaussianBasic gs = GaussianBasic(solution.factors, solution.primes);
    auto matrix = gs.Solve();
    return FindFactor(solution, matrix, n);
}

long_int ContinuedFractionsFactorization::ComputeX(const Solution & solution, const std::vector<size_t> & positions,
                                                   const long_int & n)
{
    long_int x = 1;
    for (auto pos : positions)
    {
        x = (x * solution.xs[pos]) % n;
    }
    return x;
}

long_int ContinuedFractionsFactorization::ComputeY(const Solution & solution, const std::vector<size_t> & positions,
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

bool ContinuedFractionsFactorization::IsPerfectSquare(const Line & line) { return line.IsMaskEmpty(); }

FactorSet ContinuedFractionsFactorization::FindFactor(const Solution & solution, const Matrix & matrix,
                                                      const long_int & n)
{
    FactorSet factor;
    for (const auto & line : matrix)
    {
        if (IsPerfectSquare(line))
        {
            auto positions = GaussianBasic::GetParticipantsPositions(line);
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

};  // namespace lpn
