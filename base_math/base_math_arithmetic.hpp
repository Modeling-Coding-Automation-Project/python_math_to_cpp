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
#include <cstdint>
#include <cstring>
#include <type_traits>

#ifdef __BASE_MATH_USE_STD_MATH__
#include <cmath>
#else // __BASE_MATH_USE_STD_MATH__

#ifdef __NEVER_USE_CMATH_BUT_REQUIRES_IEEE_754_STANDARD__
#else // __NEVER_USE_CMATH_BUT_REQUIRES_IEEE_754_STANDARD__
#include <cmath>
#endif // __NEVER_USE_CMATH_BUT_REQUIRES_IEEE_754_STANDARD__

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

/**
 * @brief Bitwise cast between float and uint32_t.
 *
 * These functions perform a bitwise cast between a `float` and a `uint32_t`
 * without changing the bit representation. They use `std::memcpy` to safely
 * copy the bits from one type to another.
 *
 * @param x The float value to be cast to uint32_t.
 * @return The uint32_t representation of the float's bit pattern.
 */
inline uint32_t bitcast_u32(const float &x) {
  uint32_t u;
  std::memcpy(&u, &x, sizeof(u));
  return u;
}

/**
 * @brief Bitwise cast between uint32_t and float.
 *
 * These functions perform a bitwise cast between a `uint32_t` and a `float`
 * without changing the bit representation. They use `std::memcpy` to safely
 * copy the bits from one type to another.
 *
 * @param u The uint32_t value to be cast to float.
 * @return The float representation of the uint32_t's bit pattern.
 */
inline float bitcast_f32(const uint32_t &u) {
  float x;
  std::memcpy(&x, &u, sizeof(x));
  return x;
}

/**
 * @brief Count leading zeros in a 32-bit unsigned integer.
 *
 * This function counts the number of leading zeros in the binary
 * representation of a 32-bit unsigned integer `x`. The input `x` must be
 * non-zero, as the behavior for `clz(0)` is undefined.
 *
 * @param x The 32-bit unsigned integer to count leading zeros for (must be
 * non-zero).
 * @return The number of leading zeros in `x`.
 */
inline int clz32(uint32_t x) {
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_clz(x);
#else
  // fallback (slow but correct)
  int n = 0;
  if ((x & static_cast<uint32_t>(0xFFFF0000)) == 0) {
    n += 16;
    x <<= 16;
  }
  if ((x & static_cast<uint32_t>(0xFF000000)) == 0) {
    n += 8;
    x <<= 8;
  }
  if ((x & static_cast<uint32_t>(0xF0000000)) == 0) {
    n += 4;
    x <<= 4;
  }
  if ((x & static_cast<uint32_t>(0xC0000000)) == 0) {
    n += 2;
    x <<= 2;
  }
  if ((x & static_cast<uint32_t>(0x80000000)) == 0) {
    n += 1;
  }
  return n;
#endif
}

/**
 * @brief Fast implementation of frexpf for float.
 *
 * This function decomposes a floating-point number `x` into its normalized
 * fraction and exponent. The fraction is returned, and the exponent is stored
 * in the location pointed to by `out_exp`. The function handles special cases
 * such as NaN, infinity, zero, and subnormal numbers.
 * @param x The floating-point number to decompose.
 * @param out_exp Pointer to an integer where the exponent will be stored.
 * @return The normalized fraction of `x`.
 */
inline float fast_frexpf(float x, int *out_exp) {
  uint32_t ux = bitcast_u32(x);

  const uint32_t sign = ux & static_cast<uint32_t>(0x80000000);
  const uint32_t exp = (ux >> 23) & static_cast<uint32_t>(0xFF);
  const uint32_t mant = ux & static_cast<uint32_t>(0x7FFFFF);

  // NaN/Inf
  if (exp == static_cast<uint32_t>(0xFF)) {
    if (out_exp)
      *out_exp = 0;
    return x;
  }

  // Zero or Subnormal
  if (exp == static_cast<uint32_t>(0)) {
    if (mant == static_cast<uint32_t>(0)) { // ±0
      if (out_exp)
        *out_exp = 0;
      return x;
    }

    // Subnormal: mant != 0
    // mant is 23 bits. Find the shift to align the leading 1 to bit 22.
    // English: p (0..22) is the position of the highest bit of mant, shift = 22
    // - p
    // p = 31 - clz(mant)
    const int p = 31 - clz32(mant); // mant != 0 なのでOK
    const int shift = 22 - p;       // 0..22
    const uint32_t norm_mant =
        (mant << shift) & static_cast<uint32_t>(0x7FFFFF);

    // Derivation of exp_out:
    // Subnormal numbers effectively have exp = -126 (when exponent bits are 0)
    // When normalized by shifting left 'shift' times, the value is multiplied
    // by 2^shift, so the exponent needs to be decreased by 'shift'.
    // Additionally, to make frac in [0.5,1), add 1 (same as normalized
    // numbers).
    // => out_exp = -125 - shift
    if (out_exp)
      *out_exp = -125 - shift;

    // Normalized numbers have exponent bits set to 126 (biased)
    const uint32_t frac_bits =
        sign | (static_cast<uint32_t>(126) << 23) | norm_mant;
    return bitcast_f32(frac_bits);
  }

  // Normal: 1..254
  if (out_exp)
    *out_exp = (int)exp - 126; // e - 127 + 1
  const uint32_t frac_bits = sign | (static_cast<uint32_t>(126) << 23) | mant;
  return bitcast_f32(frac_bits);
}

inline float frexpf(float x, int *exp) {
#ifdef __BASE_MATH_USE_STD_MATH__
  return std::frexp(x, exp);
#else  // __BASE_MATH_USE_STD_MATH__
  return fast_frexpf(x, exp);
#endif // __BASE_MATH_USE_STD_MATH__
}

inline float fast_ldexp_float(float x, int exp) {
  uint32_t ux = bitcast_u32(x);

  const uint32_t sign = ux & static_cast<uint32_t>(0x80000000);
  uint32_t e = (ux >> 23) & static_cast<uint32_t>(0xFF);
  uint32_t mant = ux & static_cast<uint32_t>(0x7FFFFF);

  // NaN / Inf
  if (e == static_cast<uint32_t>(0xFF)) {
    return x;
  }

  // Zero
  if (e == static_cast<uint32_t>(0) && mant == static_cast<uint32_t>(0)) {
    return x;
  }

  // ---- Normalized number ----
  if (e != static_cast<uint32_t>(0)) {
    int new_e = (int)e + exp;

    if ((unsigned)new_e >= 1u && new_e <= 254) {
      return bitcast_f32(sign | ((uint32_t)new_e << 23) | mant);
    }

    if (new_e >= 255) {
      return bitcast_f32(sign | static_cast<uint32_t>(0x7F800000)); // Inf
    }

    // Underflow -> Subnormal or 0
    // Add hidden 1 to mantissa
    mant |= static_cast<uint32_t>(0x800000);
    int shift = 1 - new_e; // right shift count

    if (shift >= 24) {
      return bitcast_f32(sign);
    }

    mant >>= shift;
    return bitcast_f32(sign | mant);
  }

  // ---- Subnormal ----
  // Normalize and then process recursively
  // mant != 0
  int p = 31 - clz32(mant);
  int shift = 22 - p;
  mant <<= shift;

  // Normalized exponent after normalization is -126 - shift
  int new_exp = exp - shift - 126;

  // Reconstruct
  int final_e = new_exp + 127;
  mant &= static_cast<uint32_t>(0x7FFFFF);

  if (final_e >= 255) {
    return bitcast_f32(sign | static_cast<uint32_t>(0x7F800000));
  }
  if (final_e <= 0) {
    if (final_e <= -23) {
      return bitcast_f32(sign);
    }
    mant |= static_cast<uint32_t>(0x800000);
    mant >>= (1 - final_e);
    return bitcast_f32(sign | mant);
  }

  return bitcast_f32(sign | (static_cast<uint32_t>(final_e) << 23) | mant);
}

/**
 * @brief Computes the ldexp (load exponent) for double precision floating-point
 * numbers.
 *
 * This function multiplies a double-precision floating-point number `x` by
 * 2 raised to the power of `exp`. It handles special cases such as NaN,
 * infinity, zero, normalized, and subnormal numbers.
 *
 * @param x The double-precision floating-point number.
 * @param exp The exponent to which 2 is raised.
 * @return The result of `x * (2 ^ exp)`.
 */
inline uint64_t bitcast_u64(const double &x) {
  uint64_t u;
  std::memcpy(&u, &x, sizeof(u));
  return u;
}

/**
 * @brief Bitwise cast between uint64_t and double.
 *
 * These functions perform a bitwise cast between a `uint64_t` and a `double`
 * without changing the bit representation. They use `std::memcpy` to safely
 * copy the bits from one type to another.
 *
 * @param u The uint64_t value to be cast to double.
 * @return The double representation of the uint64_t's bit pattern.
 */
inline double bitcast_f64(const uint64_t &u) {
  double x;
  std::memcpy(&x, &u, sizeof(x));
  return x;
}

/**
 * @brief Count leading zeros in a 64-bit unsigned integer.
 *
 * This function counts the number of leading zeros in the binary
 * representation of a 64-bit unsigned integer `x`. The input `x` must be
 * non-zero, as the behavior for `clz(0)` is undefined.
 *
 * @param x The 64-bit unsigned integer to count leading zeros for (must be
 * non-zero).
 * @return The number of leading zeros in `x`.
 */
inline int clz64(uint64_t x) {
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_clzll(x);
#else
  int n = 0;
  if ((x & static_cast<uint64_t>(0xFFFFFFFF00000000)) == 0) {
    n += 32;
    x <<= 32;
  }
  if ((x & static_cast<uint64_t>(0xFFFF000000000000)) == 0) {
    n += 16;
    x <<= 16;
  }
  if ((x & static_cast<uint64_t>(0xFF00000000000000)) == 0) {
    n += 8;
    x <<= 8;
  }
  if ((x & static_cast<uint64_t>(0xF000000000000000)) == 0) {
    n += 4;
    x <<= 4;
  }
  if ((x & static_cast<uint64_t>(0xC000000000000000)) == 0) {
    n += 2;
    x <<= 2;
  }
  if ((x & static_cast<uint64_t>(0x8000000000000000)) == 0) {
    n += 1;
  }
  return n;
#endif
}

inline double fast_ldexp_double(double x, int exp) {
  uint64_t ux = bitcast_u64(x);

  const uint64_t sign = ux & static_cast<uint64_t>(0x8000000000000000);
  uint64_t e = (ux >> 52) & static_cast<uint64_t>(0x7FF);
  uint64_t mant = ux & static_cast<uint64_t>(0xFFFFFFFFFFFFF);

  // NaN / Inf
  if (e == static_cast<uint64_t>(0x7FF)) {
    return x;
  }

  // Zero
  if (e == static_cast<uint64_t>(0) && mant == static_cast<uint64_t>(0)) {
    return x;
  }

  // ---- Normalized number --
  if (e != static_cast<uint64_t>(0)) {
    int new_e = (int)e + exp;

    // Normal range
    if ((unsigned)new_e >= 1u && new_e <= 2046) {
      return bitcast_f64(sign | ((uint64_t)new_e << 52) | mant);
    }

    // Overflow
    if (new_e >= 2047) {
      return bitcast_f64(sign | static_cast<uint64_t>(0x7FF0000000000000));
    }

    // Underflow -> subnormal / zero
    mant |= (static_cast<uint64_t>(1) << 52); // hidden bit
    int shift = 1 - new_e;                    // right shift count

    if (shift >= 53) {
      return bitcast_f64(sign); // ±0
    }

    mant >>= shift;
    return bitcast_f64(sign | mant);
  }

  // ---- Subnormal ----
  // mant != 0
  int p = 63 - clz64(mant); // highest set bit
  int shift = 52 - p;

  mant <<= shift;

  // Normalized exponent after normalization is -1022 - shift
  int new_exp = exp - shift - 1022;
  int final_e = new_exp + 1023;

  mant &= static_cast<uint64_t>(0xFFFFFFFFFFFFF);

  if (final_e >= 2047) {
    return bitcast_f64(sign | static_cast<uint64_t>(0x7FF0000000000000));
  }

  if (final_e <= 0) {
    if (final_e <= -52) {
      return bitcast_f64(sign);
    }
    mant |= (static_cast<uint64_t>(1) << 52);
    mant >>= (1 - final_e);
    return bitcast_f64(sign | mant);
  }

  return bitcast_f64(sign | ((uint64_t)final_e << 52) | mant);
}

/**
 * @brief Computes the ldexp (load exponent) for floating-point numbers.
 *
 * This function multiplies a floating-point number `x` by 2 raised to the
 * power of `exp`. It handles special cases such as NaN, infinity, zero,
 * normalized, and subnormal numbers. The implementation varies based on the
 * type of `x` (float or double) and compile-time flags.
 *
 * @tparam T The type of the floating-point number (float or double).
 * @param x The floating-point number.
 * @param exp The exponent to which 2 is raised.
 * @return The result of `x * (2 ^ exp)`.
 */
template <typename T>
typename std::enable_if<std::is_same<T, double>::value, T>::type
ldexp(T x, int exp) {
#ifdef __BASE_MATH_USE_STD_MATH__
  return std::ldexp(x, exp);
#else  // __BASE_MATH_USE_STD_MATH__
  return fast_ldexp_double(x, exp);
#endif // __BASE_MATH_USE_STD_MATH__
}

/**
 * @brief Computes the ldexp (load exponent) for single-precision floating-point
 * numbers.
 *
 * This function multiplies a single-precision floating-point number `x` by
 * 2 raised to the power of `exp`. It handles special cases such as NaN,
 * infinity, zero, normalized, and subnormal numbers.
 *
 * @param x The single-precision floating-point number.
 * @param exp The exponent to which 2 is raised.
 * @return The result of `x * (2 ^ exp)`.
 */
template <typename T>
typename std::enable_if<std::is_same<T, float>::value, T>::type ldexp(T x,
                                                                      int exp) {
#ifdef __BASE_MATH_USE_STD_MATH__
  return std::ldexp(x, exp);
#else  // __BASE_MATH_USE_STD_MATH__
  return fast_ldexp_float(x, exp);
#endif // __BASE_MATH_USE_STD_MATH__
}

} // namespace Math
} // namespace Base

#endif // __BASE_MATH_ARITHMETIC_HPP__
