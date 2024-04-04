#include <iostream>
#include "qs.hpp"
#include <chrono>

// This file is used for testing purposes.
// Will be deleted in Final commit.

int main()
{
    using namespace lpn;  // NOLINT

    long_int n("59469489332848408438249254427481121839977");
    auto start = std::chrono::high_resolution_clock::now();
    FactorSet factor = QuadraticSieve::Factorize(n);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = duration_cast<std::chrono::microseconds>(end - start);
    std::cout << duration.count() / 1e6 << std::endl;
    for (auto & [key, value] : factor)
    {
        std::cout << value << " " << key << std::endl;
    }
}
