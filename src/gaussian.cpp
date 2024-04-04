#include "gaussian.hpp"

namespace lpn
{

GaussianBasic::GaussianBasic(const std::vector<FactorSet> & factors, const std::vector<size_t> & primes)
    : n_(factors.size()), m_(primes.size()), factors_(factors)

{
    for (size_t i = 0; i < n_; ++i)
    {
        matrix_.push_back(Bitset(n_, m_));
        for (size_t j = 0; j < m_; ++j)
        {
            auto iter = factors_[i].find(primes[j]);
            if (iter != factors_[i].end())
            {
                matrix_[i].mask[j] = bool(iter->second % 2);
            }
        }
        matrix_[i].participants[i] = true;
    }
}

std::vector<GaussianBasic::Bitset> GaussianBasic::Solve()
{
    for (size_t i = 0; i < n_; ++i)
    {
        size_t pos = 0;
        for (size_t j = 0; j < m_; ++j)
        {
            if (matrix_[i].mask[j])
            {
                pos = j;
                break;
            }
        }
        Add(pos, i);
    }
    return matrix_;
}

void GaussianBasic::Add(size_t pos, size_t line)
{
    for (size_t i = 0; i < n_; ++i)
    {
        if (line != i && matrix_[i].mask[pos])
        {
            matrix_[i] ^= matrix_[line];
        }
    }
}

void GaussianBasic::Print()
{
    for (size_t i = 0; i < n_; ++i)
    {
        std::cout << matrix_[i].mask << " " << matrix_[i].participants << '\n';
    }
    std::cout << std::endl;
}

};  // namespace lpn
