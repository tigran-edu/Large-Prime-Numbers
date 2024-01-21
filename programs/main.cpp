#include <iostream>
#include "basic.h"
#include "factorization.h"

int main()
{
    using long_int = boost::multiprecision::cpp_int;

    long_int n = 0;
    long_int c = 0;
    long_int max = 0;
    // cpp_int second = 0;
    std::cin >> n >> c;
    // std::cin >> second;
    std::cout << lpn::gcd(n, c);
}
