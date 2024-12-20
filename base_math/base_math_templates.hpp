#ifndef BASE_MATH_TEMPLATES_HPP
#define BASE_MATH_TEMPLATES_HPP

#include <array>

namespace Base {
namespace Math {

template <int... Values> struct IntList {};

template <int N, typename List> struct Append;

template <int N, int... Values> struct Append<N, IntList<Values...>> {
  using type = IntList<Values..., N>;
};

template <int N> struct MakeIntList {
  using type = typename Append<N, typename MakeIntList<N - 1>::type>::type;
};

template <> struct MakeIntList<0> {
  using type = IntList<0>;
};

template <int... Values> void printList(IntList<Values...>) {
  int arr[] = {(std::cout << Values << " ", 0)...};
  (void)arr;
  std::cout << std::endl;
}

template <int... Values>
constexpr std::array<double, sizeof...(Values)> toArray(IntList<Values...>) {
  return {static_cast<double>(Values)...};
}

} // namespace Math
} // namespace Base

#endif // BASE_MATH_TEMPLATES_HPP
