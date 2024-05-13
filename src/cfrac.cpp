#include "cfrac.hpp"
#include "aliases.hpp"
#include "base.hpp"
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

Solution ContinuedFractions::Solve(const long_int & n, size_t factor_size)
{
    State state(n);
    Solution solution;
    solution.primes = FindQuadraticResiduePrimes(n, factor_size);
    while (solution.factors.size() < 1.1 * solution.primes.size())
    {
        state.Update();
        auto factor = TryToDecompose(solution.primes, state.c1);
        AddFactor(factor, state, solution);
    }
    return solution;
}

void ContinuedFractions::AddFactor(const std::optional<FactorSet> & factor, const State & state, Solution & solution)
{
    if (factor.has_value())
    {
        solution.factors.push_back(std::move(factor.value()));
        solution.values.push_back(state.p1);
    }
}

FactorSet ContinuedFractionsFactorization::Factorize(const long_int & n) { return Factorize(n, kBasicFactorSize); }

FactorSet ContinuedFractionsFactorization::Factorize(const long_int & n, size_t factor_size)
{
    auto solution = ContinuedFractions::Solve(n, factor_size);
    GaussianBasic gs = GaussianBasic(solution.factors, solution.primes);
    auto matrix = gs.Solve();
    return FindFactor(solution, matrix, n);
}

};  // namespace lpn
