/********************************************************************************
@file base_math_templates.hpp
@brief
Provides C++ template metaprogramming utilities for compile-time mathematical
computations and type list manipulations. Includes templates for generating
and manipulating lists of compile-time values, computing factorials, generating
Maclaurin series coefficients for exponential, sine, and cosine functions,
and handling alternating sign logic. These utilities enable efficient and
type-safe mathematical code generation at compile time.
********************************************************************************/
#ifndef __BASE_MATH_TEMPLATES_HPP__
#define __BASE_MATH_TEMPLATES_HPP__

#include "base_math_macros.hpp"

#include <array>

namespace Base {
namespace Math {

template <std::size_t... Values> struct ValueArgumentSizeList {};

template <int... Values> struct ValueArgumentIntList {};

template <std::size_t N, typename List> struct AppendSize;

/**
 * @brief Appends a new size value to an existing ValueArgumentSizeList.
 *
 * This template specialization defines a member type `type` which represents
 * a new ValueArgumentSizeList with the value `N` appended to the end of the
 * existing list of values (`Values...`).
 *
 * @tparam N The size value to append.
 * @tparam Values The existing list of size values.
 */
template <std::size_t N, std::size_t... Values>
struct AppendSize<N, ValueArgumentSizeList<Values...>> {
  using type = ValueArgumentSizeList<Values..., N>;
};

template <int N, typename List> struct AppendInt;

/**
 * @brief Appends an integer value N to a ValueArgumentIntList.
 *
 * This template specialization defines a type alias 'type' that represents
 * a new ValueArgumentIntList with the integer N appended to the end of the
 * existing list of integer values (Values...).
 *
 * @tparam N The integer value to append.
 * @tparam Values The existing list of integer values.
 */
template <int N, int... Values>
struct AppendInt<N, ValueArgumentIntList<Values...>> {
  using type = ValueArgumentIntList<Values..., N>;
};

/**
 * @brief Compile-time computation of factorial using template metaprogramming.
 *
 * This template struct recursively computes the factorial of a given
 * non-negative integer N at compile time. The computed value is accessible via
 * the static constexpr member `value`.
 *
 * @tparam N The non-negative integer whose factorial is to be computed.
 *
 * @note A template specialization for N = 0 should be provided to terminate the
 * recursion.
 *
 * @example
 * constexpr std::size_t fact5 = Factorial<5>::value; // fact5 == 120
 */
template <std::size_t N> struct Factorial {
  static constexpr std::size_t value = N * Factorial<N - 1>::value;
};

/**
 * @brief Template specialization for computing the factorial of 0.
 *
 * This specialization defines the base case for the Factorial template,
 * setting the value to 1, as 0! is mathematically defined to be 1.
 */
template <> struct Factorial<0> {
  static constexpr std::size_t value = static_cast<std::size_t>(1);
};

/**
 * @brief Compile-time computation of the factorial of N, multiplied by 2 at
 * each step.
 *
 * This template struct recursively computes the value of 2 * (2 * ... * (2 *
 * Factorial_2<0>::value)), effectively calculating 2^N * Factorial_2<0>::value
 * at compile time.
 *
 * @tparam N The non-negative integer for which the value is computed.
 * @note A specialization for the base case (e.g., N = 0) should be provided
 * elsewhere.
 */
template <std::size_t N> struct Factorial_2 {
  static constexpr std::size_t value =
      static_cast<std::size_t>(2) * Factorial_2<N - 1>::value;
};

/**
 * @brief Template specialization for computing the factorial of 0, multiplied
 * by 2.
 *
 * This specialization defines the base case for the Factorial_2 template,
 * setting the value to 1, as 2^0 * Factorial_2<0>::value is mathematically
 * defined to be 1.
 */
template <> struct Factorial_2<0> {
  static constexpr std::size_t value = static_cast<std::size_t>(1);
};

/* exp Maclaurin factor */

namespace ExpMaclaurinFactor {

/**
 * @brief Template metafunction to recursively construct a type list.
 *
 * This struct template, `MakeList`, generates a type list of size `N` at
 * compile time. It uses recursion to append the factorial value of `(N - 1)` to
 * the list generated for `N - 1`.
 *
 * @tparam N The size of the list to generate.
 *
 * @note Requires the definitions of `AppendSize`, `Factorial`, and the base
 * case specialization for `MakeList<0>`.
 */
template <std::size_t N> struct MakeList {
  using type = typename AppendSize<Factorial<(N - 1)>::value,
                                   typename MakeList<N - 1>::type>::type;
};

/**
 * @brief Template specialization of MakeList for N = 1.
 *
 * This specialization defines the type alias 'type' as ValueArgumentSizeList
 * instantiated with the value of Factorial<1>::value.
 *
 * Typically used in template metaprogramming to generate a type list
 * corresponding to the factorial of 1 (which is 1).
 */
template <> struct MakeList<1> {
  using type = ValueArgumentSizeList<Factorial<1>::value>;
};

/**
 * @brief Converts a ValueArgumentSizeList of compile-time integer values to a
 * std::array of doubles, where each element is the reciprocal (1.0 / value) of
 * the corresponding input value.
 *
 * @tparam Values Variadic template parameter pack representing integer values.
 * @param  Unnamed parameter of type ValueArgumentSizeList<Values...> used for
 * template deduction.
 * @return constexpr std::array<double, sizeof...(Values)> containing the
 * reciprocals of the input values.
 *
 * @note Each value in Values must be non-zero to avoid division by zero at
 * compile time.
 */
template <std::size_t... Values>
constexpr std::array<double, sizeof...(Values)>
to_array(ValueArgumentSizeList<Values...>) {
  return {static_cast<double>(static_cast<double>(1) /
                              static_cast<double>(Values))...};
}

} // namespace ExpMaclaurinFactor

/**
 * @brief Compile-time computation of alternating sign based on integer parity.
 *
 * This template struct computes a value of 1 if the template parameter N is
 * even, and -1 if N is odd. The result is available as a static constexpr
 * integer member `value`.
 *
 * @tparam N The integer whose parity determines the sign.
 *
 * @note Useful for template metaprogramming where alternating signs are needed,
 *       such as in series expansions or mathematical computations at compile
 * time.
 *
 * @example
 *   static_assert(EvenOddSign<2>::value == 1, "2 is even, so value should be
 * 1"); static_assert(EvenOddSign<3>::value == -1, "3 is odd, so value should be
 * -1");
 */
template <int N> struct EvenOddSign {
  static constexpr int value = (N % static_cast<int>(2) == static_cast<int>(0))
                                   ? static_cast<int>(1)
                                   : static_cast<int>(-1);
};

/* cos Maclaurin factor */

namespace CosMaclaurinFactor {

/**
 * @brief Generates a compile-time list of integers using template
 * metaprogramming.
 *
 * This template recursively constructs a type list by appending an integer
 * value at each step. The value appended at each recursion is calculated as:
 *   (Factorial<(2 * N)>::value) * EvenOddSign<(N + 1)>::value
 * where:
 *   - Factorial is assumed to be a template that computes the factorial of its
 * argument at compile time.
 *   - EvenOddSign is assumed to be a template that provides a sign (+1 or -1)
 * based on whether its argument is even or odd.
 *   - AppendInt is assumed to be a template that appends an integer to a type
 * list.
 *
 * @tparam N The current integer in the recursion, determines the length and
 * values of the list.
 * @note Requires the definitions of Factorial, EvenOddSign, and AppendInt
 * templates.
 */
template <int N> struct MakeList {
  using type = typename AppendInt<(static_cast<int>(Factorial<(2 * N)>::value) *
                                   EvenOddSign<(N + 1)>::value),
                                  typename MakeList<N - 1>::type>::type;
};

/**
 * @brief Template specialization for the base case of MakeList.
 *
 * This specialization defines the type alias 'type' as ValueArgumentIntList
 * instantiated with the value of (2 * Factorial<1>::value).
 *
 * Typically used in template metaprogramming to generate a type list
 * corresponding to the factorial of 1, multiplied by 2.
 */
template <> struct MakeList<1> {
  using type =
      ValueArgumentIntList<(static_cast<int>(2) * Factorial<1>::value)>;
};

/**
 * @brief Converts a ValueArgumentIntList of integer template parameters to a
 * std::array of doubles, where each element is computed as 2 divided by the
 * corresponding integer value.
 *
 * @tparam Values Variadic list of integer template parameters.
 * @param ValueArgumentIntList<Values...> A type representing a list of integer
 * values.
 * @return constexpr std::array<double, sizeof...(Values)> An array containing
 * the result of 2.0 / Value for each input value.
 */
template <int... Values>
constexpr std::array<double, sizeof...(Values)>
to_array(ValueArgumentIntList<Values...>) {
  return {static_cast<double>(static_cast<double>(2) /
                              static_cast<double>(Values))...};
}

} // namespace CosMaclaurinFactor

/* sin Maclaurin factor */

namespace SinMaclaurinFactor {

/**
 * @brief Recursively constructs a type list of integers using template
 * metaprogramming.
 *
 * For a given integer N, this template generates a type list where each element
 * is computed as: (Factorial<(2 * N - 1)>::value) * EvenOddSign<(N + 1)>::value
 * The list is built by prepending the computed value for N to the list
 * generated for N-1.
 *
 * @tparam N The number of elements to generate in the list.
 * @note Requires the definitions of Factorial, EvenOddSign, and AppendInt
 * templates.
 */
template <int N> struct MakeList {
  using type =
      typename AppendInt<(static_cast<int>(Factorial<(2 * N - 1)>::value) *
                          EvenOddSign<(N + 1)>::value),
                         typename MakeList<N - 1>::type>::type;
};

/**
 * @brief Template specialization of MakeList for N=1.
 *
 * This specialization defines the type alias 'type' as ValueArgumentIntList
 * instantiated with the value of Factorial<1>::value (which is 1).
 *
 * Typically used in template metaprogramming to generate a list containing
 * the factorial of 1 as a compile-time constant.
 */
template <> struct MakeList<1> {
  using type = ValueArgumentIntList<(Factorial<1>::value)>;
};

/**
 * @brief Converts a ValueArgumentIntList of integer template parameters to a
 * std::array of doubles, where each element is the reciprocal (1.0 / Value) of
 * the corresponding integer value.
 *
 * @tparam Values Variadic list of integer template parameters.
 * @param ValueArgumentIntList<Values...> A type representing a list of integer
 * values.
 * @return constexpr std::array<double, sizeof...(Values)> An array containing
 * the reciprocals of the input values as doubles.
 */
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
