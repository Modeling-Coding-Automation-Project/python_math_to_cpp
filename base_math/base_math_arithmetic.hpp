/**
 * @file base_math_arithmetic.hpp
 * @brief Provides basic arithmetic utility functions for mathematical
 * operations in the Base::Math namespace.
 *
 * This header defines a set of constexpr constants and template functions to
 * perform basic arithmetic operations, such as absolute value and modulus, with
 * optional use of standard math library functions. The implementation is
 * designed to be generic and work with various numeric types, supporting both
 * custom and standard math operations based on compile-time flags.
 *
 * @namespace Base
 *     Root namespace for the base library components.
 *
 * @namespace Base::Math
 *     Contains mathematical utility functions and constants for arithmetic
 * operations.
 */
#ifndef __BASE_MATH_ARITHMETIC_HPP__
#define __BASE_MATH_ARITHMETIC_HPP__

#include "base_math_macros.hpp"

#include "base_math_mathematical_constants.hpp"

#include <cstddef>

#ifdef __BASE_MATH_USE_STD_MATH__
#include <cmath>
#else  // __BASE_MATH_USE_STD_MATH__
#endif // __BASE_MATH_USE_STD_MATH__

namespace Base {
namespace Math {

constexpr double ARITHMETIC_DIVISION_MIN = 1.0e-10;

constexpr std::size_t MOD_REPEAT_MAX_NUMBER = 100;

/* abs */
/**
 * @brief Computes the absolute value of a given number.
 *
 * This function template returns the absolute value of the input parameter `x`.
 * If `x` is less than zero, it returns `-x`; otherwise, it returns `x`.
 *
 * @tparam T Numeric type supporting comparison and negation.
 * @param x The value whose absolute value is to be computed.
 * @return The absolute value of `x`.
 */
template <typename T> inline T abs_base_math(const T &x) {
  return x < static_cast<T>(0) ? -x : x;
}

/**
 * @brief Computes the absolute value of the given input.
 *
 * This function returns the absolute value of the input parameter `x`.
 * If the macro `__BASE_MATH_USE_STD_MATH__` is defined, it uses `std::abs`.
 * Otherwise, it uses the custom implementation `Base::Math::abs_base_math`.
 *
 * @tparam T The type of the input value.
 * @param x The value whose absolute value is to be computed.
 * @return The absolute value of `x`.
 */
template <typename T> inline T abs(const T &x) {

#ifdef __BASE_MATH_USE_STD_MATH__
  return std::abs(x);
#else  // __BASE_MATH_USE_STD_MATH__
  return Base::Math::abs_base_math(x);
#endif // __BASE_MATH_USE_STD_MATH__
}

/* mod */

/**
 * @brief Computes the modulus of two values with integer casting and custom
 * division minimum.
 *
 * This function calculates the modulus (remainder) of x divided by y, with
 * additional logic:
 * - If y is less than division_min, the function returns 0.
 * - The absolute value of y is used for the division.
 * - The quotient is computed by dividing x by abs(y), casting the result to
 * int, and then back to T.
 * - The result is calculated as x minus the product of the quotient and abs(y).
 *
 * @tparam T Numeric type supporting arithmetic operations and casting.
 * @param x The dividend.
 * @param y The divisor.
 * @param division_min The minimum value for the divisor; if y < division_min,
 * returns 0.
 * @return The modulus of x and y, or 0 if y < division_min.
 */
template <typename T>
inline T mod_with_casting_integer(const T &x, const T &y,
                                  const T &division_min) {
  if (y < division_min) {
    return static_cast<T>(0);
  }

  T abs_y = Base::Math::abs_base_math(y);
  T quotient = static_cast<T>(static_cast<int>(x / abs_y));
  T result = x - quotient * abs_y;

  return result;
}

/**
 * @brief Computes the modulus (remainder) of two values.
 *
 * This function returns the result of the modulus operation between `x` and
 * `y`. If the macro `__BASE_MATH_USE_STD_MATH__` is defined, it uses
 * `std::fmod` for floating-point types. Otherwise, it uses a custom
 * implementation via `Base::Math::mod_with_casting_integer`.
 *
 * @tparam T The type of the operands (should support arithmetic operations).
 * @param x The dividend.
 * @param y The divisor.
 * @return The remainder after dividing `x` by `y`.
 */
template <typename T> inline T mod(const T &x, const T &y) {

#ifdef __BASE_MATH_USE_STD_MATH__
  return std::fmod(x, y);
#else  // __BASE_MATH_USE_STD_MATH__
  return Base::Math::mod_with_casting_integer(
      x, y, static_cast<T>(Base::Math::ARITHMETIC_DIVISION_MIN));
#endif // __BASE_MATH_USE_STD_MATH__
}

} // namespace Math
} // namespace Base

#endif // __BASE_MATH_ARITHMETIC_HPP__
