#include <vector>
#include <boost/dynamic_bitset.hpp>
#include "basic.hpp"

namespace lpn
{
class GaussianBasic
{
   public:
    struct Bitset
    {
        Bitset(size_t n, size_t m) : participants(n), mask(m) {}

        void operator^=(Bitset & bitset)
        {
            mask ^= bitset.mask;
            participants ^= bitset.participants;
        }

        boost::dynamic_bitset<> participants;
        boost::dynamic_bitset<> mask;
    };

    GaussianBasic(const std::vector<Factor> & factors, const std::vector<size_t> & primes)
        : n_(factors.size()), m_(primes.size()), factors_(factors)

    {
        for (size_t i = 0; i < n_; ++i)
        {
            matrix_.push_back(Bitset(n_, m_));
            for (size_t j = 0; j < m_; ++j)
            {
                matrix_[i].mask[j] = (factors_[i][primes[j]]) % 2 == 1;
            }
            matrix_[i].participants[i] = true;
        }
    }

    std::vector<Bitset> Solve()
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
        // Print();
        return matrix_;
    }

    void Add(size_t pos, size_t line)
    {
        for (size_t i = 0; i < n_; ++i)
        {
            if (line != i && matrix_[i].mask[pos])
            {
                matrix_[i] ^= matrix_[line];
            }
        }
    }

    void Print()
    {
        for (size_t i = 0; i < n_; ++i)
        {
            std::cout << matrix_[i].mask << " " << matrix_[i].participants << '\n';
        }
        std::cout << std::endl;
    }

   private:
    size_t m_;
    size_t n_;
    std::vector<Factor> factors_;
    std::vector<Bitset> matrix_;
};
};  // namespace lpn
