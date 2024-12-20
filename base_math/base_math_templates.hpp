#ifndef BASE_MATH_TEMPLATES_HPP
#define BASE_MATH_TEMPLATES_HPP

#include <array>

namespace Base {
namespace Math {

template <std::size_t... Values> struct DoubleValueList {};

template <std::size_t N, typename List> struct Append;

template <std::size_t N, std::size_t... Values>
struct Append<N, DoubleValueList<Values...>> {
  using type = DoubleValueList<Values..., N>;
};

template <std::size_t N> struct Factorial {
  static constexpr std::size_t value = N * Factorial<N - 1>::value;
};

template <> struct Factorial<0> {
  static constexpr std::size_t value = static_cast<std::size_t>(1);
};

template <std::size_t N> struct MakeExpMcloughlinFactorList {
  using type =
      typename Append<Factorial<(N - 1)>::value,
                      typename MakeExpMcloughlinFactorList<N - 1>::type>::type;
};

template <> struct MakeExpMcloughlinFactorList<1> {
  using type = DoubleValueList<Factorial<1>::value>;
};

template <std::size_t... Values>
constexpr std::array<double, sizeof...(Values)>
toArray(DoubleValueList<Values...>) {
  return {static_cast<double>(Values)...};
}

} // namespace Math
} // namespace Base

#endif // BASE_MATH_TEMPLATES_HPP
