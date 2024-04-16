#include <iostream>
#include "aliases.hpp"
#include "elliptic.hpp"
#include <chrono>
#include <boost/random.hpp>

// This file is used for testing purposes.
// Will be deleted in Final commit.

int main()
{
    using namespace lpn;  // NOLINT

    long_int n("231923777552225902498567");
    boost::mt19937 rng(43);
    FactorSet factor;
    while (factor.size() < 2)
    {
        auto start = std::chrono::high_resolution_clock::now();
        factor = EllipticCurveFactorization::Factorize(n, rng() % n, rng() % n, rng() % n, 1'000'000);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = duration_cast<std::chrono::microseconds>(end - start);
        std::cout << duration.count() / 1e6 << std::endl;
        for (auto & [key, value] : factor)
        {
            std::cout << value << " " << key << std::endl;
        }
    }
}
