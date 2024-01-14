#include <iostream>
#include "basic.h"
#include "factorization.h"

using boost::multiprecision::cpp_int;

int main()
{
    cpp_int n = 0;
    cpp_int c = 0;
    cpp_int max = 0;
    // cpp_int second = 0;
    std::cin >> n >> c >> max;
    // std::cin >> second;
    auto factor = lpn::RhoFactorization(n, c, max);
    for (auto & [dev, amount] : factor)
        std::cout << dev << " " << amount << '\n';
}
