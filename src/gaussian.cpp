#include "gaussian.hpp"

namespace lpn
{
using Line = GaussianBasic::Line;
using Matrix = GaussianBasic::Matrix;

GaussianBasic::GaussianBasic(const FactorSets & factors, const std::vector<size_t> & primes)
    : n_(factors.size()), m_(primes.size()), matrix_(CreateMatrix(factors, primes)), factors_(factors)
{
}

Matrix GaussianBasic::CreateMatrix(const FactorSets & factors, const std::vector<size_t> & primes)
{
    size_t n = factors.size();
    size_t m = primes.size();

    Matrix matrix;
    for (size_t i = 0; i < n; ++i)
    {
        matrix.push_back(Line(n, m));
        for (size_t j = 0; j < m; ++j)
        {
            auto iter = factors[i].find(primes[j]);
            if (iter != factors[i].end())
            {
                matrix[i].mask[j] = bool(iter->second % 2);
            }
        }
        matrix[i].participants[i] = true;
    }
    return matrix;
}

size_t GaussianBasic::FindFirstNonZeroInLine(size_t line_pos)
{
    for (size_t col = 0; col < m_; ++col)
    {
        if (matrix_[line_pos].mask[col])
        {
            return col;
        }
    }
    return m_;
}

std::vector<Line> GaussianBasic::Solve()
{
    for (size_t i = 0; i < n_; ++i)
    {
        size_t col = FindFirstNonZeroInLine(i);
        if (col != m_)
        {
            Add(col, i);
        }
    }
    return matrix_;
}

void GaussianBasic::Add(size_t col, size_t line)
{
    for (size_t i = 0; i < n_; ++i)
    {
        if (line != i && matrix_[i].mask[col])
        {
            matrix_[i] ^= matrix_[line];
        }
    }
}

};  // namespace lpn
