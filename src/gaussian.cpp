#include "gaussian.hpp"
#include <algorithm>

namespace lpn
{
GaussianBasic::Line::Line(size_t n, size_t m) : participants(n), mask(m) {}

GaussianBasic::Line & GaussianBasic::Line::operator^=(Line & line)
{
    mask ^= line.mask;
    participants ^= line.participants;
    return *this;
}

bool GaussianBasic::Line::IsMaskEmpty() const { return mask.none(); }

GaussianBasic::GaussianBasic(const FactorSets & factors, const std::vector<size_t> & primes)
    : n_(factors.size()), m_(primes.size()), matrix_(CreateMatrix(factors, primes))
{
}

GaussianBasic::Matrix GaussianBasic::CreateMatrix(const FactorSets & factors, const std::vector<size_t> & primes)
{
    size_t n = factors.size();
    size_t m = primes.size();

    Matrix matrix = Matrix(n);
    for (size_t i = 0; i < n; ++i)
    {
        matrix[i] = Line(n, m);
        for (size_t j = 0; j < m; ++j)
        {
            auto iter = factors[i].find(primes[j]);
            if (iter != factors[i].end() && iter->second % 2 == 1)
            {
                matrix[i].mask[j] = true;
            }
        }
        matrix[i].participants[i] = true;
    }
    return matrix;
}

size_t GaussianBasic::FindFirstNonZeroInLine(size_t line_pos) const
{
    size_t col = matrix_[line_pos].mask.find_first();
    if (col == std::string::npos)
    {
        return m_;
    }
    return col;
}

GaussianBasic::Matrix GaussianBasic::Solve()
{
    for (size_t i = 0; i < n_; ++i)
    {
        size_t col = FindFirstNonZeroInLine(i);
        if (col != m_)
        {
            AddToAll(col, i);
        }
    }
    return matrix_;
}

void GaussianBasic::AddToAll(size_t col, size_t line)
{
    for (size_t i = 0; i < n_; ++i)
    {
        if (matrix_[i].mask[col] && i != line)
        {
            matrix_[i] ^= matrix_[line];
        }
    }
}

};  // namespace lpn
