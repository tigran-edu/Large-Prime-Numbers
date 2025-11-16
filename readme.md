Large Prime Numbers
======
The course project was completed by student Andrian Tigran, HSE Faculty of Computer Science, Applied Mathematics and Information Science.
-----
## Build
The library is built using the [Cmake](https://cmake.org) version 3.22 or higher. To use the LPN, the pre-installed library [Boost](https://www.boost.org/) version 1.84.0 or higher is also needed.
## Basic Concepts
All the main aliases are located in the file [aliases.hpp](https://github.com/tigran-edu/Large-Prime-Numbers/blob/main/src/aliases.hpp).
* **FactorSet** - set of divisors of a number with powers
* **FactorSets** - set consisting of FactorSet
* **long_int** - type with an infinite number of digits

# Factorization algorithms

## Trial Division
The basic algo located in the [basic.hpp](https://github.com/tigran-edu/Large-Prime-Numbers/blob/main/src/basic.hpp#L86C11-L86C25). The function FactorizeBasic from `basic.hpp` provides a basic factorization algorithm.
```c++
FactorSet FactorizeBasic(long_int a);
```

## Pollard Rho 
The method is located in the file [rho.npp](https://github.com/tigran-edu/Large-Prime-Numbers/blob/main/src/rho.hpp#L10C1-L11C1). The class `RhoFactorization` provides 2 static public methods:
```c++
static FactorSet Factorize(const long_int & n, size_t starting_point);
static FactorSet Factorize(const long_int & n, const long_int & c, size_t max_iter, size_t frequency,
                             size_t max_attemp, size_t starting_point);
```
The first one uses basic configuration:
```c++
struct BasicConfig
{
    kDefCValue = 3;
    kPrimes = {2, 3, 5, 7, 11, 13};
    kDefMaxIterValue = 1000000;
    kDefFreqValue = 10;
    kDefMaxAttempValue = 100;
};
```

If the user requires more precise tuning, there is the second public method available that allow the algorithm to be fully customized for a specific input number n.


## Quadratic Sieve
The method is located in the file [qs.hpp](https://github.com/tigran-edu/Large-Prime-Numbers/blob/main/src/qs.hpp#L70C7-L70C34). The class `QuadraticSieveFactorization` provides 2 static public methods:
```c++
static FactorSet Factorize(const long_int & n);
static FactorSet Factorize(const long_int & n, const Sieve::Config & config);
```
The first one uses basic configuration of the Sieve:
```c++
struct BasicConfig
{
    kDefSegmentSize = 50'000'000;
    kDefFactorSize = 2000;
    kExpansionRate = 1.5;
};
```
If the user requires more precise tuning, there is the second public method available that allow the algorithm to be fully customized by a config.
### Config
* `segment_size` - the size of sieve segment, segment - half of the sieve size
* `factor_size` -  the size of the set of prime numbers used in the algorithm
* `multi_thread` - enable parallel
* `expansion_rate` - the coefficient that sets the minimum limit of the value in the sieve
  
The config is created through a static method in the class `Sieve`:
```c++
  static Config CreateConfig(const long_int & n, size_t segment_size, size_t factor_size, float expansion_rate,
                             bool multi_thread = false);
```

## Continued Fraction
The method is located in the file [cfrac.hpp](https://github.com/tigran-edu/Large-Prime-Numbers/blob/main/src/cfrac.hpp#L38). The class `ContinuedFractionsFactorization` provides 2 static public methods:
```c++
static FactorSet Factorize(const long_int & n);
static FactorSet Factorize(const long_int & n, size_t factor_size);
```
The first one uses default factor_size = 4000

## Elliptic Curve
The method is located in the file [elliptic.hpp](https://github.com/tigran-edu/Large-Prime-Numbers/blob/main/src/elliptic.hpp#L8). The class `EllipticCurveFactorization` provides static public method:
```c++
static FactorSet Factorize(const long_int & n, const long_int & a, const long_int & x, const long_int & y,
                           size_t max_iter);
```
### Parameters
* `a` - elliptic curve parameter
*  `x` и `y` - coordinate on the curve
*  `max_iter` - maximum number of iterations

# Checking numbers for primality
The library LPN provides 2 methods for the checking numbers for primality:
- PseudoPrimeTest
- StrongPseudoPrimeTest
  
Methods are located in the file [basic.hpp](https://github.com/tigran-edu/Large-Prime-Numbers/blob/main/src/basic.hpp).
```c++
template <typename Container>
bool IsPseudoPrime(const long_int & p, const Container & primes);

template <typename Container>
bool IsStrongPseudoPrime(const long_int & p, const Container & primes);
```
