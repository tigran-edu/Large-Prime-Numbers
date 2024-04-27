#include <iostream>
#include "qs.hpp"
#include <thread>
#include "elliptic.hpp"
#include <chrono>
#include <boost/random.hpp>

#include <fstream>

// This file is used for testing purposes.
// Will be deleted in Final commit.

int main()
{
    using namespace lpn;  // NOLINT
    long_int n("804032999113025075736319759489");

    auto start = std::chrono::high_resolution_clock::now();
    auto config = Sieve::CreateConfig(n, 1600000, 700, 0.4);
    FactorSet factor = QuadraticSieveFactorization::Factorize(n, config);
    // auto end = std::chrono::high_resolution_clock::now();
    // auto duration = duration_cast<std::chrono::microseconds>(end - start);
    // std::cout << duration.count() / 1e6 << std::endl;
    // long_int n("231923777552225902498567");
    // boost::mt19937 rng(1233);
    // FactorSet factor;
    // auto start = std::chrono::high_resolution_clock::now();
    // while (factor.size() < 2)
    // {
    //     auto local_start = std::chrono::high_resolution_clock::now();
    //     factor = EllipticCurveFactorization::Factorize(n, rng() % n, rng() % n, rng() % n, 20'000);
    //     auto local_end = std::chrono::high_resolution_clock::now();
    //     auto duration = duration_cast<std::chrono::microseconds>(local_end - local_start);
    //     std::cout << duration.count() / 1e6 << std::endl;
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = duration_cast<std::chrono::microseconds>(end - start);
    std::cout << duration.count() / 1e6 << std::endl;
    for (auto & [key, value] : factor)
    {
        std::cout << value << " " << key << std::endl;
    }
    // }
    // auto end = std::chrono::high_resolution_clock::now();
    // auto duration = duration_cast<std::chrono::microseconds>(end - start);
    // std::cout << duration.count() / 1e6 << std::endl;
}
