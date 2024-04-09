#include <iostream>
#include "cfrac.hpp"
#include <chrono>

// This file is used for testing purposes.
// Will be deleted in Final commit.

int main()
{
    using namespace lpn;  // NOLINT

    long_int n("4482406424966880742829846540605971439398287609");
    auto start = std::chrono::high_resolution_clock::now();
    FactorSet factor = ContinuedFractionsFactorization::Factorize(n, 1000);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = duration_cast<std::chrono::microseconds>(end - start);
    std::cout << duration.count() / 1e6 << std::endl;
    for (auto & [key, value] : factor)
    {
        std::cout << value << " " << key << std::endl;
    }
}
