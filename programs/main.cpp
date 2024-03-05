#include <iostream>
#include "basic.h"
#include "factorization.h"

int main()
{
    using long_int = boost::multiprecision::cpp_int;  // NOLINT

    long_int value("3456456345674323456765432345654323456543234565432345678976");
    long_int p( "8292399238022656234968869231588487594247212499421683182301788083552859865836334850585856382505692477");
    std::cout << -1 % 3;
}
