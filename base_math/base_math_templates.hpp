#ifndef __BASE_MATH_TEMPLATES_HPP__
#define __BASE_MATH_TEMPLATES_HPP__

#include "base_math_macros.hpp"

#include <array>

namespace Base {
namespace Math {

template <std::size_t... Values> struct ValueArgumentSizeList {};

template <int... Values> struct ValueArgumentIntList {};

template <std::size_t N, typename List> struct AppendSize;

template <std::size_t N, std::size_t... Values>
struct AppendSize<N, ValueArgumentSizeList<Values...>> {
  using type = ValueArgumentSizeList<Values..., N>;
};

template <int N, typename List> struct AppendInt;

template <int N, int... Values>
struct AppendInt<N, ValueArgumentIntList<Values...>> {
  using type = ValueArgumentIntList<Values..., N>;
};

template <std::size_t N> struct Factorial {
  static constexpr std::size_t value = N * Factorial<N - 1>::value;
};

template <> struct Factorial<0> {
  static constexpr std::size_t value = static_cast<std::size_t>(1);
};

template <std::size_t N> struct Factorial_2 {
  static constexpr std::size_t value =
      static_cast<std::size_t>(2) * Factorial_2<N - 1>::value;
};

template <> struct Factorial_2<0> {
  static constexpr std::size_t value = static_cast<std::size_t>(1);
};

/* exp Maclaurin factor */
namespace ExpMaclaurinFactor {

template <std::size_t N> struct MakeList {
  using type = typename AppendSize<Factorial<(N - 1)>::value,
                                   typename MakeList<N - 1>::type>::type;
};

template <> struct MakeList<1> {
  using type = ValueArgumentSizeList<Factorial<1>::value>;
};

template <std::size_t... Values>
constexpr std::array<double, sizeof...(Values)>
to_array(ValueArgumentSizeList<Values...>) {
  return {static_cast<double>(static_cast<double>(1) /
                              static_cast<double>(Values))...};
}

} // namespace ExpMaclaurinFactor

template <int N> struct EvenOddSign {
  static constexpr int value = (N % static_cast<int>(2) == static_cast<int>(0))
                                   ? static_cast<int>(1)
                                   : static_cast<int>(-1);
};

/* cos Maclaurin factor */
namespace CosMaclaurinFactor {

template <int N> struct MakeList {
  using type = typename AppendInt<(static_cast<int>(Factorial<(2 * N)>::value) *
                                   EvenOddSign<(N + 1)>::value),
                                  typename MakeList<N - 1>::type>::type;
};

template <> struct MakeList<1> {
  using type =
      ValueArgumentIntList<(static_cast<int>(2) * Factorial<1>::value)>;
};

template <int... Values>
constexpr std::array<double, sizeof...(Values)>
to_array(ValueArgumentIntList<Values...>) {
  return {static_cast<double>(static_cast<double>(2) /
                              static_cast<double>(Values))...};
}

} // namespace CosMaclaurinFactor

/* sin Maclaurin factor */
namespace SinMaclaurinFactor {

template <int N> struct MakeList {
  using type =
      typename AppendInt<(static_cast<int>(Factorial<(2 * N - 1)>::value) *
                          EvenOddSign<(N + 1)>::value),
                         typename MakeList<N - 1>::type>::type;
};

template <> struct MakeList<1> {
  using type = ValueArgumentIntList<(Factorial<1>::value)>;
};

template <int... Values>
constexpr std::array<double, sizeof...(Values)>
to_array(ValueArgumentIntList<Values...>) {
  return {static_cast<double>(static_cast<double>(1) /
                              static_cast<double>(Values))...};
}

} // namespace SinMaclaurinFactor

} // namespace Math
} // namespace Base

#endif // __BASE_MATH_TEMPLATES_HPP__
