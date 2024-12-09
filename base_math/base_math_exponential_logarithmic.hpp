#ifndef BASE_MATH_EXPONENTIAL_LOGARITHMIC_HPP
#define BASE_MATH_EXPONENTIAL_LOGARITHMIC_HPP

#include "base_math_arithmetic.hpp"
#include "base_math_macros.hpp"
#include "base_math_mathematical_constants.hpp"
#include "base_math_utility.hpp"

#include <cstddef>
#include <cstring>

#ifdef BASE_MATH_USE_STD_MATH
#include <cmath>
#else  // BASE_MATH_USE_STD_MATH
#endif // BASE_MATH_USE_STD_MATH

namespace Base {
namespace Math {

constexpr double EXPONENTIAL_LOGARITHMIC_DIVISION_MIN = 1.0e-10;

#ifdef BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS
constexpr std::size_t SQRT_REPEAT_NUMBER = 4;
constexpr std::size_t EXP_REPEAT_NUMBER = 4;
constexpr std::size_t EXP2_REPEAT_NUMBER = 4;
constexpr std::size_t LOG_REPEAT_NUMBER = 5;
#else  // BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS
constexpr std::size_t SQRT_REPEAT_NUMBER = 5;
constexpr std::size_t EXP_REPEAT_NUMBER = 7;
constexpr std::size_t EXP2_REPEAT_NUMBER = 8;
constexpr std::size_t LOG_REPEAT_NUMBER = 7;
#endif // BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS

constexpr double EXP_INPUT_MAX = 87.0;
constexpr double EXP_INPUT_MIN = -87.0;
constexpr double EXP2_INPUT_MAX = 126.0;
constexpr double EXP2_INPUT_MIN = -126.0;

constexpr double LOG_OUTPUT_MIN = -1.0e20;
constexpr double LOG_SCALE_FACTOR = 5.0;
constexpr double LOG_OF_LOG_SCALE_FACTOR =
    1.609437912434100; // ln(LOG_SCALE_FACTOR)
constexpr std::size_t LOG_SCALING_LOOP_MAX =
    54; // ln(1e38) / ln(LOG_SCALE_FACTOR)

/* sqrt */
template <typename T> inline T sqrt_base_math(const T &x) {
  T result = static_cast<T>(0);

  if (x <= static_cast<T>(0)) {
    /* Do Nothing. */
  } else {
    T guess = x / static_cast<T>(2);

    for (std::size_t i = 0; i < Base::Math::SQRT_REPEAT_NUMBER; ++i) {
      guess = (guess + x / guess) / static_cast<T>(2);
    }

    result = guess;
  }

  return result;
}

/*************************************************
*  Detail: Calculates 1 / sqrt(x) within an error margin of 0.2%.
*          Explanation on Wikipedia:
*          https://en.wikipedia.org/wiki/Fast_inverse_square_root
*          This function made significant achievements in 3D rendering in the
early game industry.
*          The function here is based on the wiki function with variable names
changed.
*          Also, *(float *)& has been changed to memcpy.
*          The calculation speed is about 4 times faster than math.h's sqrt
(depending on the environment).
*          Inputting 0.0F outputs 1.98e+19, but avoid inputting negative values
as much as possible.
*          If the absolute value is small, 0 will be output, but inputting
values less than -1.0F may result in inf.
**************************************************/
template <typename T> inline T fast_inverse_square_root(T input) {

  float y = static_cast<float>(0);

  if (input <= static_cast<T>(0)) {
    /* Do Nothing. */
  } else {
    long i = static_cast<long>(0);
    float x = static_cast<float>(0);

    x = static_cast<float>(input) * static_cast<float>(0.5);
    y = static_cast<float>(input);
    std::memcpy(&i, &y, 4);
    i = static_cast<long>(0x5f3759df) - (i >> 1);
    std::memcpy(&y, &i, 4);
    y = y * (static_cast<float>(1.5) - (x * y * y));
  }

  return static_cast<T>(y);
}

template <typename T> inline T fast_square_root(T input) {

  if (input <= EXPONENTIAL_LOGARITHMIC_DIVISION_MIN) {
    return static_cast<T>(0);
  } else {
    return static_cast<T>(1) / Base::Math::fast_inverse_square_root(input);
  }
}

template <typename T> inline T sqrt(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::sqrt(x);
#else
#ifdef BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS

#ifdef BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD
  return Base::Math::fast_square_root(x);
#else  // BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD
  return Base::Math::sqrt_base_math(x);
#endif // BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD

#else // BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS
  return Base::Math::sqrt_base_math(x);

#endif // BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS
#endif // BASE_MATH_USE_STD_MATH
}

/* exp */
template <typename T> inline T exp_base_math(const T &x) {

  T x_wrapped = x;

  if (x > static_cast<T>(Base::Math::EXP_INPUT_MAX)) {
    x_wrapped = static_cast<T>(Base::Math::EXP_INPUT_MAX);
  } else if (x < static_cast<T>(Base::Math::EXP_INPUT_MIN)) {
    x_wrapped = static_cast<T>(Base::Math::EXP_INPUT_MIN);
  }

  int n = static_cast<int>(x_wrapped / static_cast<T>(Base::Math::LN_2));
  T remainder = x_wrapped - n * static_cast<T>(Base::Math::LN_2);

  T result = static_cast<T>(1);
  T term = static_cast<T>(1);

  for (std::size_t i = 1; i < Base::Math::EXP_REPEAT_NUMBER; ++i) {
    term *= remainder / static_cast<T>(i);
    result += term;
  }

  if (n > 0) {
    for (int i = 0; i < n; ++i) {
      result *= static_cast<T>(2);
    }
  } else {
    for (int i = 0; i < -n; ++i) {
      result *= static_cast<T>(0.5);
    }
  }

  return result;
}

template <typename T> inline T exp(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::exp(x);
#else
  return Base::Math::exp_base_math(x);
#endif
}

/* exp2 */
template <typename T> inline T exp2_base_math(const T &x) {

  T x_wrapped = x;

  if (x > static_cast<T>(Base::Math::EXP2_INPUT_MAX)) {
    x_wrapped = static_cast<T>(Base::Math::EXP2_INPUT_MAX);
  } else if (x < static_cast<T>(Base::Math::EXP2_INPUT_MIN)) {
    x_wrapped = static_cast<T>(Base::Math::EXP2_INPUT_MIN);
  }

  int n = static_cast<int>(x_wrapped);
  T remainder = x_wrapped - static_cast<T>(n);

  T result = static_cast<T>(1);
  T term = static_cast<T>(1);

  for (std::size_t i = 1; i < Base::Math::EXP2_REPEAT_NUMBER; ++i) {
    term *= static_cast<T>(Base::Math::LN_2) * remainder / static_cast<T>(i);
    result += term;
  }

  if (n > 0) {
    for (int i = 0; i < n; ++i) {
      result *= static_cast<T>(2);
    }
  } else {
    for (int i = 0; i < -n; ++i) {
      result *= static_cast<T>(0.5);
    }
  }

  return result;
}

template <typename T> inline T exp2(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::exp2(x);
#else
  return Base::Math::exp2_base_math(x);
#endif
}

/* log */
template <typename T> inline T log_base_math(const T &x) {
  T result = static_cast<T>(0);

  if (x <= static_cast<T>(0)) {
    result = static_cast<T>(Base::Math::LOG_OUTPUT_MIN);
  } else {

    T scaled_x = x;
    int scale = 0;
    std::size_t loop_count = 0;

    for (; scaled_x > Base::Math::LOG_SCALE_FACTOR &&
           loop_count < Base::Math::LOG_SCALING_LOOP_MAX;
         ++scale, ++loop_count) {
      scaled_x /= Base::Math::LOG_SCALE_FACTOR;
    }

    loop_count = 0;
    for (; scaled_x < static_cast<T>(1) / Base::Math::LOG_SCALE_FACTOR &&
           loop_count < Base::Math::LOG_SCALING_LOOP_MAX;
         --scale, ++loop_count) {
      scaled_x *= Base::Math::LOG_SCALE_FACTOR;
    }

    T guess = scaled_x - static_cast<T>(1);
    T exp_guess = static_cast<T>(0);
    for (std::size_t i = 0; i < Base::Math::LOG_REPEAT_NUMBER; ++i) {
      exp_guess = std::exp(guess);
      guess = guess - (exp_guess - scaled_x) / exp_guess;
    }

    result = guess + static_cast<T>(scale) *
                         static_cast<T>(Base::Math::LOG_OF_LOG_SCALE_FACTOR);
  }

  return result;
}

template <typename T> inline T log(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::log(x);
#else
  return Base::Math::log_base_math(x);
#endif
}

/* log2 */
template <typename T> inline T log2_base_math(const T &x) {

  return Base::Math::log_base_math(x) / static_cast<T>(Base::Math::LN_2);
}

template <typename T> inline T log2(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::log2(x);
#else
  return Base::Math::log2_base_math(x);
#endif
}

/* log10 */
template <typename T> inline T log10_base_math(const T &x) {

  return Base::Math::log_base_math(x) / static_cast<T>(Base::Math::LN_10);
}

template <typename T> inline T log10(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::log10(x);
#else
  return Base::Math::log10_base_math(x);
#endif
}

/* pow */
template <typename T> inline T pow_base_math(const T &x, const T &y) {
  return Base::Math::exp_base_math(y * Base::Math::log_base_math(x));
}

template <typename T> inline T pow(const T &x, const T &y) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::pow(x, y);
#else
  return Base::Math::pow_base_math(x, y);
#endif
}

} // namespace Math
} // namespace Base

#endif // BASE_MATH_EXPONENTIAL_LOGARITHMIC_HPP
