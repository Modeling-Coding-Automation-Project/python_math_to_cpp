#ifndef BASE_MATH_TEMPLATES_HPP
#define BASE_MATH_TEMPLATES_HPP

#include <array>

namespace Base {
namespace Math {

template <std::size_t... Values> struct ValueArgumentList {};

template <std::size_t N, typename List> struct Append;

template <std::size_t N, std::size_t... Values>
struct Append<N, ValueArgumentList<Values...>> {
  using type = ValueArgumentList<Values..., N>;
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
  using type = ValueArgumentList<Factorial<1>::value>;
};

template <std::size_t... Values>
constexpr std::array<double, sizeof...(Values)>
to_exp_mcloughlin_factor_array(ValueArgumentList<Values...>) {
  return {static_cast<double>(static_cast<double>(1) /
                              static_cast<double>(Values))...};
}

} // namespace Math
} // namespace Base

#endif // BASE_MATH_TEMPLATES_HPP
