#include "base.hpp"

namespace lpn
{

bool FactorizationBase::IsPerfectSquare(const Line & line) { return line.IsMaskEmpty(); }

std::vector<size_t> FactorizationBase::GetParticipantsPositions(const Line & line)
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

FactorSet FactorizationBase::FindFactor(const Solution & solution, const Matrix & matrix, const long_int & n)
{
    FactorSet factor;
    for (const auto & line : matrix)
    {
        if (IsPerfectSquare(line))
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

long_int FactorizationBase::ComputeY(const Solution & solution, const std::vector<size_t> & positions,
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

long_int FactorizationBase::ComputeX(const Solution & solution, const std::vector<size_t> & positions,
                                     const long_int & n)
{
    long_int x = 1;
    for (auto pos : positions)
    {
        x = (x * solution.values[pos]) % n;
    }
    return x;
}

};  // namespace lpn
