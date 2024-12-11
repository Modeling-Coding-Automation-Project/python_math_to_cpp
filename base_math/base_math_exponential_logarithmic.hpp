#ifndef BASE_MATH_EXPONENTIAL_LOGARITHMIC_HPP
#define BASE_MATH_EXPONENTIAL_LOGARITHMIC_HPP

#include "base_math_arithmetic.hpp"
#include "base_math_macros.hpp"
#include "base_math_mathematical_constants.hpp"
#include "base_math_utility.hpp"

#include <cmath>
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
constexpr std::size_t SQRT_REPEAT_NUMBER = 0;
constexpr std::size_t EXP_REPEAT_NUMBER = 4;
constexpr std::size_t EXP2_REPEAT_NUMBER = 4;
constexpr std::size_t LOG_REPEAT_NUMBER = 5;
#else  // BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS
constexpr std::size_t SQRT_REPEAT_NUMBER = 1;
constexpr std::size_t EXP_REPEAT_NUMBER = 7;
constexpr std::size_t EXP2_REPEAT_NUMBER = 8;
constexpr std::size_t LOG_REPEAT_NUMBER = 7;
#endif // BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS

constexpr std::size_t SQRT_REPEAT_NUMBER_MOSTLY_ACCURATE = 2;

constexpr float ONE_AND_SQRT2_VEC[2] = {1.0F, 1.4142136F};

constexpr double EXP_INPUT_MAX = 87.0;
constexpr double EXP_OUTPUT_MAX = 6.076030225056872e+37;
constexpr double EXP_INPUT_MIN = -87.0;
constexpr double EXP_OUTPUT_MIN = 1.6458114310822737e-38;

constexpr double EXP2_INPUT_MAX = 126.0;
constexpr double EXP2_OUTPUT_MAX = 8.507059173023462e+37;
constexpr double EXP2_INPUT_MIN = -126.0;
constexpr double EXP2_OUTPUT_MIN = 1.1754943508222875e-38;

constexpr double LOG_OUTPUT_MIN = -1.0e20;
constexpr double LOG_SCALE_FACTOR = 5.0;
constexpr double LOG_OF_LOG_SCALE_FACTOR =
    1.609437912434100; // ln(LOG_SCALE_FACTOR)
constexpr std::size_t LOG_SCALING_LOOP_MAX =
    54; // ln(1e38) / ln(LOG_SCALE_FACTOR)

/* inverse sqrt */
/*************************************************
*  Function: fast_inverse_square_root
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

template <typename T, int N> struct RsqrtBaseMathLoop {
  static void compute(T x_T, T &result) {
    result *= static_cast<T>(1.5) - x_T * result * result;
    RsqrtBaseMathLoop<T, N - 1>::compute(x_T, result);
  }
};

template <typename T> struct RsqrtBaseMathLoop<T, 1> {
  static void compute(T x_T, T &result) {
    result *= static_cast<T>(1.5) - x_T * result * result;
  }
};

template <typename T> struct RsqrtBaseMathLoop<T, 0> {
  static void compute(T x_T, T &result) {
    /* Do Nothing. */
    static_cast<void>(x_T);
    static_cast<void>(result);
  }
};

template <typename T, std::size_t LOOP_NUMBER>
inline T rsqrt_base_math(const T &x, const T &division_min) {
  T result = static_cast<T>(0);

  T x_wrapped = x;
  if (x < division_min) {
    x_wrapped = division_min;
  }

  float x_float = static_cast<float>(x_wrapped);

  int e = 0;
  T h = static_cast<T>(0);
  float r = 1.8284271F - 0.82842712F * std::frexpf(x_float, &e);

  r = Base::Math::ldexp(
      r * Base::Math::ONE_AND_SQRT2_VEC[e & static_cast<int>(0x00000001)],
      -e >> 1);

  result = static_cast<T>(r);

  h = x * result * result;
  result *= static_cast<T>(1.875) -
            h * (static_cast<T>(1.25) - h * static_cast<T>(0.375));

  RsqrtBaseMathLoop<T, LOOP_NUMBER>::compute(x_wrapped * static_cast<T>(0.5),
                                             result);

  return result;
}

template <typename T> inline T rsqrt(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return static_cast<T>(1) / std::sqrt(x);
#else
#ifdef BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS

#ifdef BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD
  return Base::Math::fast_inverse_square_root(x);
#else  // BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD
  return Base::Math::rsqrt_base_math<T, Base::Math::SQRT_REPEAT_NUMBER>(
      x, static_cast<T>(Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN));
#endif // BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD

#else // BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS
  return Base::Math::rsqrt_base_math<T, Base::Math::SQRT_REPEAT_NUMBER>(
      x, static_cast<T>(Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN));

#endif // BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS
#endif // BASE_MATH_USE_STD_MATH
}

/* sqrt */
template <typename T, std::size_t LOOP_NUMBER>
inline T sqrt_base_math(const T &x) {

  T x_wrapped = x;

  if (x < static_cast<T>(Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN)) {
    x_wrapped =
        static_cast<T>(Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN);
  }

  return x_wrapped *
         Base::Math::rsqrt_base_math<T, LOOP_NUMBER>(
             x_wrapped,
             static_cast<T>(Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN));
}

template <typename T, std::size_t LOOP_NUMBER>
inline T sqrt_base_math(const T &x, const T &division_min) {

  T x_wrapped = x;

  if (x < division_min) {
    x_wrapped = division_min;
  }

  return x_wrapped *
         Base::Math::rsqrt_base_math<T, LOOP_NUMBER>(x_wrapped, division_min);
}

template <typename T> inline T fast_square_root(T input) {

  if (input <= EXPONENTIAL_LOGARITHMIC_DIVISION_MIN) {
    return static_cast<T>(0);
  } else {
    return input * Base::Math::fast_inverse_square_root(input);
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
  return Base::Math::sqrt_base_math<T, Base::Math::SQRT_REPEAT_NUMBER>(x);
#endif // BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD

#else // BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS
  return Base::Math::sqrt_base_math<T, Base::Math::SQRT_REPEAT_NUMBER>(x);

#endif // BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS
#endif // BASE_MATH_USE_STD_MATH
}

/* exp */
template <typename T, std::size_t LOOP_MAX, std::size_t N>
struct ExpIterationLoop {
  static void compute(T &result, T &term, const T &remainder) {
    term *= remainder / static_cast<T>(LOOP_MAX - N);
    result += term;

    ExpIterationLoop<T, LOOP_MAX, N - 1>::compute(result, term, remainder);
  }
};

template <typename T, std::size_t LOOP_MAX>
struct ExpIterationLoop<T, LOOP_MAX, 0> {
  static void compute(T &result, T &term, const T &remainder) {
    /* Do Nothing. */
    static_cast<void>(result);
    static_cast<void>(term);
    static_cast<void>(remainder);
  }
};

template <typename T, std::size_t LOOP_NUMBER>
inline T exp_base_math(const T &x) {

  T result = static_cast<T>(1);

  if (x > static_cast<T>(Base::Math::EXP_INPUT_MAX)) {
    result = static_cast<T>(Base::Math::EXP_OUTPUT_MAX);
  } else if (x < static_cast<T>(Base::Math::EXP_INPUT_MIN)) {
    result = static_cast<T>(Base::Math::EXP_OUTPUT_MIN);
  } else {

    int n = static_cast<int>(x / static_cast<T>(Base::Math::LN_2));
    T remainder = x - n * static_cast<T>(Base::Math::LN_2);
    T term = static_cast<T>(1);

    ExpIterationLoop<T, LOOP_NUMBER, LOOP_NUMBER - 1>::compute(result, term,
                                                               remainder);

    result = Base::Math::ldexp(result, n);
  }

  return result;
}

template <typename T> inline T exp(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::exp(x);
#else
  return Base::Math::exp_base_math<T, Base::Math::EXP_REPEAT_NUMBER>(x);
#endif
}

/* exp2 */
template <typename T, std::size_t LOOP_MAX, std::size_t N>
struct Exp2IterationLoop {
  static void compute(T &result, T &term, const T &remainder) {
    term *= static_cast<T>(Base::Math::LN_2) * remainder /
            static_cast<T>(LOOP_MAX - N);
    result += term;

    Exp2IterationLoop<T, LOOP_MAX, N - 1>::compute(result, term, remainder);
  }
};

template <typename T, std::size_t LOOP_MAX>
struct Exp2IterationLoop<T, LOOP_MAX, 0> {
  static void compute(T &result, T &term, const T &remainder) {
    /* Do Nothing. */
    static_cast<void>(result);
    static_cast<void>(term);
    static_cast<void>(remainder);
  }
};

template <typename T, std::size_t LOOP_NUMBER>
inline T exp2_base_math(const T &x) {

  T result = static_cast<T>(1);

  if (x > static_cast<T>(Base::Math::EXP2_INPUT_MAX)) {
    result = static_cast<T>(Base::Math::EXP2_OUTPUT_MAX);
  } else if (x < static_cast<T>(Base::Math::EXP2_INPUT_MIN)) {
    result = static_cast<T>(Base::Math::EXP2_OUTPUT_MIN);
  } else {

    int n = static_cast<int>(x);
    T remainder = x - static_cast<T>(n);

    T term = static_cast<T>(1);

    Exp2IterationLoop<T, LOOP_NUMBER, LOOP_NUMBER - 1>::compute(result, term,
                                                                remainder);

    result = Base::Math::ldexp(result, n);
  }

  return result;
}

template <typename T> inline T exp2(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::exp2(x);
#else
  return Base::Math::exp2_base_math<T, Base::Math::EXP2_REPEAT_NUMBER>(x);
#endif
}

/* log */
template <typename T, std::size_t N> struct LogIterationLoop {
  static void compute(T &exp_guess, T &guess, const T &scaled_x) {
    exp_guess = std::exp(guess);
    guess = guess - (exp_guess - scaled_x) / exp_guess;

    LogIterationLoop<T, N - 1>::compute(exp_guess, guess, scaled_x);
  }
};

template <typename T> struct LogIterationLoop<T, 0> {
  static void compute(T &exp_guess, T &guess, const T &scaled_x) {
    /* Do Nothing. */
    static_cast<void>(exp_guess);
    static_cast<void>(guess);
    static_cast<void>(scaled_x);
  }
};

template <typename T, std::size_t LOOP_NUMBER>
inline T log_base_math(const T &x) {
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

    LogIterationLoop<T, LOOP_NUMBER - 1>::compute(exp_guess, guess, scaled_x);

    result = guess + static_cast<T>(scale) *
                         static_cast<T>(Base::Math::LOG_OF_LOG_SCALE_FACTOR);
  }

  return result;
}

template <typename T> inline T log(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::log(x);
#else
  return Base::Math::log_base_math<T, Base::Math::LOG_REPEAT_NUMBER>(x);
#endif
}

/* log2 */
template <typename T, std::size_t LOOP_NUMBER>
inline T log2_base_math(const T &x) {

  return Base::Math::log_base_math<T, LOOP_NUMBER>(x) /
         static_cast<T>(Base::Math::LN_2);
}

template <typename T> inline T log2(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::log2(x);
#else
  return Base::Math::log2_base_math<T, Base::Math::LOG_REPEAT_NUMBER>(x);
#endif
}

/* log10 */
template <typename T, std::size_t LOOP_NUMBER>
inline T log10_base_math(const T &x) {

  return Base::Math::log_base_math<T, LOOP_NUMBER>(x) /
         static_cast<T>(Base::Math::LN_10);
}

template <typename T> inline T log10(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::log10(x);
#else
  return Base::Math::log10_base_math<T, Base::Math::LOG_REPEAT_NUMBER>(x);
#endif
}

/* pow */
template <typename T, std::size_t EXP_LOOP_NUMBER, std::size_t LOG_LOOP_NUMBER>
inline T pow_base_math(const T &x, const T &y) {
  return Base::Math::exp_base_math<T, EXP_LOOP_NUMBER>(
      y * Base::Math::log_base_math<T, LOG_LOOP_NUMBER>(x));
}

template <typename T> inline T pow(const T &x, const T &y) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::pow(x, y);
#else
  return Base::Math::pow_base_math<T, Base::Math::EXP_REPEAT_NUMBER,
                                   Base::Math::LOG_REPEAT_NUMBER>(x, y);
#endif
}

} // namespace Math
} // namespace Base

#endif // BASE_MATH_EXPONENTIAL_LOGARITHMIC_HPP
