Курсовой проект Large Prime Numbers
======
Курсовой проект выполнен студентом Андряном Тиграном, ВШЭ ФКН ПМИ БПМИ215.
-----
Данный раздел является инструкцией по использованию библиотеки LPN.
## Сборка
Сборка библиотеки происходит при помощи утилиты [Cmake](https://cmake.org) версии 3.22 и выше. Для использования библиотеки также понадобится предустановленная библиотека [Boost](https://www.boost.org/) версии 1.84.0.
## Основные понятия
Все основные aliases находятся в файле [aliases.hpp](https://github.com/tigran-edu/Large-Prime-Numbers/blob/main/src/aliases.hpp).
* **FactorSet** - набор делителей числа со степенями
* **FactorSets** - набор из FactorSet
* **long_int** - тип целых чисел с бесконечным количеством знаков

# Алгоритмы факторизации

## Trial Division
Базовые алгоритмы находятся в файле [basic.hpp](https://github.com/tigran-edu/Large-Prime-Numbers/blob/main/src/basic.hpp#L86C11-L86C25). Функция FactorizeBasic из `basic.hpp` предоставляет базовый алгоритм факторизации.
```c++
FactorSet FactorizeBasic(long_int a);
```

## Pollard Rho 
Данный алгоритм располагается в файле [rho.npp](https://github.com/tigran-edu/Large-Prime-Numbers/blob/main/src/rho.hpp#L10C1-L11C1). Класс `RhoFactorization` предоставляет 2 публичных статических метода:
```c++
static FactorSet Factorize(const long_int & n, size_t starting_point);
static FactorSet Factorize(const long_int & n, const long_int & c, size_t max_iter, size_t frequency,
                             size_t max_attemp, size_t starting_point);
```
Первая функция подразумевает использование базовой конфигурации алгоритма:
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
Если пользователь желает более точной настройки, то для этого предоставлен 2 публичный метод, который позволяет полностью настроить алгоритм под определенное вводимое число n.


## Quadratic Sieve
Алгоритм представлен в файле [qs.hpp](https://github.com/tigran-edu/Large-Prime-Numbers/blob/main/src/qs.hpp#L70C7-L70C34). Класс `QuadraticSieveFactorization` предоставляет 2 публичных статических метода:
```c++
static FactorSet Factorize(const long_int & n);
static FactorSet Factorize(const long_int & n, const Sieve::Config & config);
```
Первый метод подразумевает использование базового конфига решета:
```c++
struct BasicConfig
{
    kDefSegmentSize = 50'000'000;
    kDefFactorSize = 2000;
    kExpansionRate = 1.5;
};
```
Если же пользователь хочет более тонкой настройки, то для этого предоставлен 2 публичный метод, который принимает конфиг.
### Config
* `segment_size` - размер сегмента решета, сегмент - половина решета
* `factor_size` - размер множетсва простых чисел, используемых в алгоритме
* `multi_thread` - булевая переменная, при значении True использует многопоточное решение решета
* `expansion_rate` - коэффициент, устанавливающий минимальную границу значения в решете
  
Конфиг создается через статических метод в классе `Sieve`:
```c++
  static Config CreateConfig(const long_int & n, size_t segment_size, size_t factor_size, float expansion_rate,
                             bool multi_thread = false);
```

## Continued Fraction
Алгоритм расположен в файле [cfrac.hpp](https://github.com/tigran-edu/Large-Prime-Numbers/blob/main/src/cfrac.hpp#L38). Класс `ContinuedFractionsFactorization` предоставляет 2 публичных статических метода:
```c++
static FactorSet Factorize(const long_int & n);
static FactorSet Factorize(const long_int & n, size_t factor_size);
```
Первый метод использует базовый размер множества простых чисел - 4000.

## Elliptic Curve
Алгоритм расположен в файле [elliptic.hpp](https://github.com/tigran-edu/Large-Prime-Numbers/blob/main/src/elliptic.hpp#L8). Класс `EllipticCurveFactorization` предоставляет публичный статический метод:
```c++
static FactorSet Factorize(const long_int & n, const long_int & a, const long_int & x, const long_int & y,
                           size_t max_iter);
```
### Параметры
* `a` - параметр эллиптической кривой
*  `x` и `y` - координата на кривой
*  `max_iter` - максимальное количество итераций алгоритма

# Проверка чисел на простоту
Библиотека LPN предоставляет 2 метода проверки чисел на простоту:
- PseudoPrimeTest
- StrongPseudoPrimeTest
  
Оба алгоритма находятся в файле [basic.hpp](https://github.com/tigran-edu/Large-Prime-Numbers/blob/main/src/basic.hpp).
```c++
template <typename Container>
bool IsPseudoPrime(const long_int & p, const Container & primes);

template <typename Container>
bool IsStrongPseudoPrime(const long_int & p, const Container & primes);
```
