#include <iostream>
#include "qs.hpp"

int main()
{
    using namespace lpn;  // NOLINT
    long_int n("59469489332848408438249254427481121839977");
    Factor factor = QuadraticSieve::Factorize(n);
    for (auto & [key, value] : factor)
    {
        std::cout << value << " " << key << std::endl;
    }
}
