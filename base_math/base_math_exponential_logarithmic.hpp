#ifndef BASE_MATH_EXPONENTIAL_LOGARITHMIC_HPP
#define BASE_MATH_EXPONENTIAL_LOGARITHMIC_HPP

#include "base_math_arithmetic.hpp"
#include "base_math_macros.hpp"
#include "base_math_mathematical_constants.hpp"
#include "base_math_templates.hpp"

#include <cmath>
#include <cstddef>
#include <cstring>

namespace Base {
namespace Math {

constexpr double EXPONENTIAL_LOGARITHMIC_DIVISION_MIN = 1.0e-10;

#ifdef BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS
constexpr std::size_t SQRT_REPEAT_NUMBER = 0;
constexpr std::size_t EXP_REPEAT_NUMBER = 5;
constexpr std::size_t EXP2_REPEAT_NUMBER = 4;
constexpr std::size_t LOG_REPEAT_NUMBER = 6;

constexpr int SQRT_EXTRACTION_DOUBLE_REPEAT_NUMBER = -21;
constexpr int SQRT_EXTRACTION_FLOAT_REPEAT_NUMBER = -6;
#else  // BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS
constexpr std::size_t SQRT_REPEAT_NUMBER = 1;
constexpr std::size_t EXP_REPEAT_NUMBER = 9;
constexpr std::size_t EXP2_REPEAT_NUMBER = 8;
constexpr std::size_t LOG_REPEAT_NUMBER = 7;

constexpr int SQRT_EXTRACTION_DOUBLE_REPEAT_NUMBER = -7;
constexpr int SQRT_EXTRACTION_FLOAT_REPEAT_NUMBER = 8;
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

constexpr unsigned long long TABLE_FOR_EXP[16] = {
    0x059B0D3158574ull, // 2^( 1 /32)-1
    0x11301D0125B51ull, // 2^( 3 /32)-1
    0x1D4873168B9AAull, // 2^( 5 /32)-1
    0x29E9DF51FDEE1ull, // 2^( 7 /32)-1
    0x371A7373AA9CBull, // 2^( 9 /32)-1
    0x44E086061892Dull, // 2^( 11 /32)-1
    0x5342B569D4F82ull, // 2^( 13 /32)-1
    0x6247EB03A5585ull, // 2^( 15 /32)-1
    0x71F75E8EC5F74ull, // 2^( 17 /32)-1
    0x82589994CCE13ull, // 2^( 19 /32)-1
    0x93737B0CDC5E5ull, // 2^( 21 /32)-1
    0xA5503B23E255Dull, // 2^( 23 /32)-1
    0xB7F76F2FB5E47ull, // 2^( 25 /32)-1
    0xCB720DCEF9069ull, // 2^( 27 /32)-1
    0xDFC97337B9B5Full, // 2^( 29 /32)-1
    0xF50765B6E4540ull, // 2^( 31 /32)-1
};

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
template <typename T>
inline T fast_inverse_square_root(const T &input, const T &division_min) {

  float y = static_cast<float>(0);

  T input_wrapped = input;
  if (input < division_min) {
    input_wrapped = division_min;
  }

  if (input <= static_cast<T>(0)) {
    /* Do Nothing. */
  } else {
    long i = static_cast<long>(0);
    float x = static_cast<float>(0);

    x = static_cast<float>(input_wrapped) * static_cast<float>(0.5);
    y = static_cast<float>(input_wrapped);
    std::memcpy(&i, &y, 4);
    i = static_cast<long>(0x5f3759df) - (i >> 1);
    std::memcpy(&y, &i, 4);
    y = y * (static_cast<float>(1.5) - (x * y * y));
  }

  return static_cast<T>(y);
}

/* sqrt extraction double */
template <int Loop_Limit, int I, int I_Limit>
struct SqrtExtractionDoubleFirstLoop {
  static void execute(unsigned long long &c, unsigned long long &m,
                      unsigned long long &y, unsigned long long &a) {
    c = static_cast<unsigned long>(
        (y << static_cast<int>(1) | static_cast<unsigned long>(1)) <= (m >> I));
    a = (a << static_cast<int>(1)) | c;
    y = (y << static_cast<int>(1)) | c;
    m -= (c * y) << I;
    y += c;

    SqrtExtractionDoubleFirstLoop<Loop_Limit, (I - 2),
                                  (I - 2 - Loop_Limit)>::execute(c, m, y, a);
  }
};

template <int Loop_Limit, int I_Limit>
struct SqrtExtractionDoubleFirstLoop<Loop_Limit, -2, I_Limit> {
  static void execute(unsigned long long &c, unsigned long long &m,
                      unsigned long long &y, unsigned long long &a) {
    static_cast<void>(c);
    static_cast<void>(m);
    static_cast<void>(y);
    static_cast<void>(a);
    // Base case: do nothing
  }
};

template <int Loop_Limit, int I>
struct SqrtExtractionDoubleFirstLoop<Loop_Limit, I, -2> {
  static void execute(unsigned long long &c, unsigned long long &m,
                      unsigned long long &y, unsigned long long &a) {
    static_cast<void>(c);
    static_cast<void>(m);
    static_cast<void>(y);
    static_cast<void>(a);
    // Base case: do nothing
  }
};

template <int Loop_Limit>
struct SqrtExtractionDoubleFirstLoop<Loop_Limit, -2, -2> {
  static void execute(unsigned long long &c, unsigned long long &m,
                      unsigned long long &y, unsigned long long &a) {
    static_cast<void>(c);
    static_cast<void>(m);
    static_cast<void>(y);
    static_cast<void>(a);
    // Base case: do nothing
  }
};

template <int I> struct SqrtExtractionDoubleSecondLoop {
  static void execute(unsigned long long &c, unsigned long long &m,
                      unsigned long long &y, unsigned long long &a) {
    m <<= static_cast<int>(2);
    c = static_cast<unsigned long long>(
        (y << static_cast<int>(1) | static_cast<unsigned long long>(1)) <= m);
    a = (a << static_cast<int>(1)) | c;
    y = (y << static_cast<int>(1)) | c;
    m -= (c * y);
    y += c;
    SqrtExtractionDoubleSecondLoop<I - 1>::execute(c, m, y, a);
  }
};

template <> struct SqrtExtractionDoubleSecondLoop<0> {
  static void execute(unsigned long long &c, unsigned long long &m,
                      unsigned long long &y, unsigned long long &a) {
    m <<= static_cast<int>(2);
    c = static_cast<unsigned long long>(
        (y << static_cast<int>(1) | static_cast<unsigned long long>(1)) <= m);
    a = (a << static_cast<int>(1)) | c;
    y = (y << static_cast<int>(1)) | c;
    m -= (c * y);
    y += c;
  }
};

template <> struct SqrtExtractionDoubleSecondLoop<-1> {
  static void execute(unsigned long long &c, unsigned long long &m,
                      unsigned long long &y, unsigned long long &a) {
    static_cast<void>(c);
    static_cast<void>(m);
    static_cast<void>(y);
    static_cast<void>(a);
    // Base case: do nothing
  }
};

/*************************************************
 *  Function: sqrt_extraction_double
 *  Detail: calculates sqrt with extraction of mantissa.
 *  Input value depends on the IEEE 754 standard.
 *  When In_r is specified as the maximum 26,
 *  it completely matches math.h's sqrt().
 *  When In_r is 0, the error is within 8.0e^-7 [%].
 *  When In_r is -10, the error is within 8.0e^-3 [%].
 *  When In_r is -20, the error is within 8 [%].
 *  At -26, since almost no calculation is performed, the error is very large.
 **************************************************/
template <int EXTRACTION_DOUBLE_REPEAT_NUMBER>
inline double sqrt_extraction_double(const double &value) {

  double result = static_cast<double>(0);
  unsigned short n =
      static_cast<unsigned short>(0); // n is the exponent of value
  unsigned long long m =
      static_cast<unsigned long long>(0); // m is the mantissa of value
  unsigned long long Temp;

  unsigned long long a = static_cast<unsigned long long>(0);
  unsigned long long c = static_cast<unsigned long long>(0);
  unsigned long long y = static_cast<unsigned long long>(0);

  std::memcpy(&Temp, &value, static_cast<std::size_t>(8));

  // To calculate the square root, the square root of 2 generated by the power
  // of 2 is assigned to the mantissa, and the value is adjusted.
  // The exponent is 2 to the power of (n - 1024), and when the square root is
  // taken, it becomes 2 to the power of (n/2 - 512).
  // In the following, the case is divided into even and odd numbers so that n
  // does not become a fraction when divided by 2.
  // When n is odd
  // The mantissa is (m + 2^52), and the exponent is (n - 1) / 2 + 511
  // When n is even
  // The mantissa is (m + 2^52) * 2, and the exponent is n / 2 + 511
  // However, since 1 is added to the exponent due to the rounding up of the
  // mantissa later, 510 is added here
  n = static_cast<unsigned short>(
      ((static_cast<unsigned long long>(0x7FFFFFFFFFFFFFFF) & Temp) >>
       static_cast<int>(52)));

  m = (static_cast<unsigned long long>(0x10000000000000) +
       (static_cast<unsigned long long>(0xFFFFFFFFFFFFF) & Temp))
      << (static_cast<unsigned short>(1) & (~n));

  // the exponent part to Temp
  Temp = static_cast<unsigned long long>(
             (((n + static_cast<unsigned short>(
                        (static_cast<unsigned short>(1) & n))) >>
               static_cast<int>(1)) +
              static_cast<unsigned short>(510)))
         << static_cast<int>(52);

  SqrtExtractionDoubleFirstLoop<
      ((static_cast<int>(-2) * EXTRACTION_DOUBLE_REPEAT_NUMBER -
        static_cast<int>(2)) *
       (EXTRACTION_DOUBLE_REPEAT_NUMBER < static_cast<int>(-1))),
      54,
      (54 - ((static_cast<int>(-2) * EXTRACTION_DOUBLE_REPEAT_NUMBER -
              static_cast<int>(2)) *
             (EXTRACTION_DOUBLE_REPEAT_NUMBER <
              static_cast<int>(-1))))>::execute(c, m, y, a);

  SqrtExtractionDoubleSecondLoop<(
      (EXTRACTION_DOUBLE_REPEAT_NUMBER *
       static_cast<int>(
           (EXTRACTION_DOUBLE_REPEAT_NUMBER >= static_cast<int>(0)))) +
      static_cast<int>(-1) *
          (static_cast<int>((EXTRACTION_DOUBLE_REPEAT_NUMBER <
                             static_cast<int>(0)))))>::execute(c, m, y, a);

  // Round the least significant digit (1 carry up).
  // Therefore, the above square root calculation is performed one more time.
  // The smaller the number of iterations, the greater the left shift amount.
  // Add the mantissa to Temp.
  Temp += (((a + (static_cast<unsigned long long>(1) & a)) >> 1)
           << (static_cast<int>(26) - EXTRACTION_DOUBLE_REPEAT_NUMBER));

  std::memcpy(&result, &Temp, static_cast<std::size_t>(8));

  return result;
}

/* sqrt extraction float */
template <int Loop_Limit, int I, int I_Limit>
struct SqrtExtractionFloatFirstLoop {
  static void execute(unsigned long &c, unsigned long &m, unsigned long &y,
                      unsigned long &a) {
    c = static_cast<unsigned long>(
        (y << static_cast<int>(1) | static_cast<unsigned long>(1)) <= (m >> I));
    a = (a << static_cast<int>(1)) | c;
    y = (y << static_cast<int>(1)) | c;
    m -= (c * y) << I;
    y += c;

    SqrtExtractionFloatFirstLoop<Loop_Limit, (I - 2),
                                 (I - 2 - Loop_Limit)>::execute(c, m, y, a);
  }
};

template <int Loop_Limit, int I_Limit>
struct SqrtExtractionFloatFirstLoop<Loop_Limit, -1, I_Limit> {
  static void execute(unsigned long &c, unsigned long &m, unsigned long &y,
                      unsigned long &a) {
    static_cast<void>(c);
    static_cast<void>(m);
    static_cast<void>(y);
    static_cast<void>(a);
    // Base case: do nothing
  }
};

template <int Loop_Limit, int I>
struct SqrtExtractionFloatFirstLoop<Loop_Limit, I, -1> {
  static void execute(unsigned long &c, unsigned long &m, unsigned long &y,
                      unsigned long &a) {
    static_cast<void>(c);
    static_cast<void>(m);
    static_cast<void>(y);
    static_cast<void>(a);
    // Base case: do nothing
  }
};

template <int Loop_Limit>
struct SqrtExtractionFloatFirstLoop<Loop_Limit, -1, -1> {
  static void execute(unsigned long &c, unsigned long &m, unsigned long &y,
                      unsigned long &a) {
    static_cast<void>(c);
    static_cast<void>(m);
    static_cast<void>(y);
    static_cast<void>(a);
    // Base case: do nothing
  }
};

template <int I> struct SqrtExtractionFloatSecondLoop {
  static void execute(unsigned long &c, unsigned long &m, unsigned long &y,
                      unsigned long &a) {
    m <<= static_cast<int>(2);
    c = static_cast<unsigned long>(
        (y << static_cast<int>(1) | static_cast<unsigned long>(1)) <= m);
    a = (a << static_cast<int>(1)) | c;
    y = (y << static_cast<int>(1)) | c;
    m -= (c * y);
    y += c;
    SqrtExtractionFloatSecondLoop<I - 1>::execute(c, m, y, a);
  }
};

template <> struct SqrtExtractionFloatSecondLoop<0> {
  static void execute(unsigned long &c, unsigned long &m, unsigned long &y,
                      unsigned long &a) {
    m <<= static_cast<int>(2);
    c = static_cast<unsigned long>(
        (y << static_cast<int>(1) | static_cast<unsigned long>(1)) <= m);
    a = (a << static_cast<int>(1)) | c;
    y = (y << static_cast<int>(1)) | c;
    m -= (c * y);
    y += c;
  }
};

template <> struct SqrtExtractionFloatSecondLoop<-1> {
  static void execute(unsigned long &c, unsigned long &m, unsigned long &y,
                      unsigned long &a) {
    static_cast<void>(c);
    static_cast<void>(m);
    static_cast<void>(y);
    static_cast<void>(a);
    // Base case: do nothing
  }
};

/*************************************************
 *  Function: sqrt_extraction_float
 *  Detail: calculates sqrt with extraction of mantissa.
 *  check the comment of double version.
 **************************************************/
template <int EXTRACTION_FLOAT_REPEAT_NUMBER>
inline float sqrt_extraction_float(const float &value) {

  float result = static_cast<float>(0);
  unsigned short n = static_cast<short>(0); // n is the exponent of value
  unsigned long m = static_cast<long>(0);   // m is the mantissa of value
  unsigned long Temp = static_cast<long>(0);

  unsigned long a = static_cast<unsigned long>(0);
  unsigned long c = static_cast<unsigned long>(0);
  unsigned long y = static_cast<unsigned long>(0);

  std::memcpy(&Temp, &value, static_cast<std::size_t>(4));

  // The exponent is 2 to the power of  (n - 128), and when the square root is
  // taken, it becomes 2 to the power of (n/2 - 64).
  // When n is odd
  // The mantissa is (m + 2^23), and the exponent is (n - 1) / 2 + 63
  // When n is even
  // The mantissa is (m + 2^23) * 2, and the exponent is n / 2 + 63
  // However, since 1 is added to the exponent due to the rounding up of the
  // mantissa later, 62 is added here
  n = static_cast<unsigned short>((
      (static_cast<unsigned long>(0x7FFFFFFF) & Temp) >> static_cast<int>(23)));

  m = (static_cast<unsigned long long>(0x800000) +
       (static_cast<unsigned long long>(0x7FFFFF) & Temp))
      << (static_cast<unsigned short>(1) & (~n));

  // the exponent part to Temp
  Temp = static_cast<unsigned long>(
             (((n + static_cast<unsigned short>(
                        (static_cast<unsigned short>(1) & n))) >>
               static_cast<int>(1)) +
              static_cast<unsigned short>(62)))
         << static_cast<int>(23);

  // Use the square root calculation to calculate the square root of m, which
  // has been updated above. The resulting square root has the number of digits
  // for the number of iterations.
  SqrtExtractionFloatFirstLoop<
      ((static_cast<int>(-2) * EXTRACTION_FLOAT_REPEAT_NUMBER -
        static_cast<int>(2)) *
       (EXTRACTION_FLOAT_REPEAT_NUMBER < static_cast<int>(-1))),
      25,
      (25 - ((static_cast<int>(-2) * EXTRACTION_FLOAT_REPEAT_NUMBER -
              static_cast<int>(2)) *
             (EXTRACTION_FLOAT_REPEAT_NUMBER <
              static_cast<int>(-1))))>::execute(c, m, y, a);

  y <<= static_cast<int>(1);

  SqrtExtractionFloatSecondLoop<(
      (EXTRACTION_FLOAT_REPEAT_NUMBER *
       static_cast<int>(
           (EXTRACTION_FLOAT_REPEAT_NUMBER >= static_cast<int>(0)))) +
      static_cast<int>(-1) *
          (static_cast<int>((EXTRACTION_FLOAT_REPEAT_NUMBER <
                             static_cast<int>(0)))))>::execute(c, m, y, a);

  // Round the least significant digit (1 carry up).
  // Therefore, the above square root calculation is performed one more time.
  // The smaller the number of iterations, the greater the left shift amount.
  // Add the mantissa to Temp.
  Temp += (((a + (static_cast<unsigned long>(1) & a)) >> 1)
           << (static_cast<int>(12) - EXTRACTION_FLOAT_REPEAT_NUMBER));

  std::memcpy(&result, &Temp, static_cast<std::size_t>(4));

  return result;
}

/* rsqrt loop */
template <typename T, int N> struct RsqrtNewtonLoop {
  static void compute(T x_T, T &result) {
    result *= static_cast<T>(1.5) - x_T * result * result;
    RsqrtNewtonLoop<T, N - 1>::compute(x_T, result);
  }
};

template <typename T> struct RsqrtNewtonLoop<T, 1> {
  static void compute(T x_T, T &result) {
    result *= static_cast<T>(1.5) - x_T * result * result;
  }
};

template <typename T> struct RsqrtNewtonLoop<T, 0> {
  static void compute(T x_T, T &result) {
    /* Do Nothing. */
    static_cast<void>(x_T);
    static_cast<void>(result);
  }
};

template <typename T, std::size_t LOOP_NUMBER>
inline T rsqrt_newton_method(const T &x, const T &division_min) {
  T result = static_cast<T>(0);

  T x_wrapped = x;
  if (x < division_min) {
    x_wrapped = division_min;
  }

  float x_float = static_cast<float>(x_wrapped);

  int e = 0;
  T h = static_cast<T>(0);
  float r = 1.8284271F - 0.82842712F * frexpf(x_float, &e);

  r = ldexp(r * Base::Math::ONE_AND_SQRT2_VEC[e & static_cast<int>(0x00000001)],
            -e >> 1);

  result = static_cast<T>(r);

  h = x * result * result;
  result *= static_cast<T>(1.875) -
            h * (static_cast<T>(1.25) - h * static_cast<T>(0.375));

  RsqrtNewtonLoop<T, LOOP_NUMBER>::compute(x_wrapped * static_cast<T>(0.5),
                                           result);

  return result;
}

template <typename T> inline T rsqrt(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return static_cast<T>(1) / std::sqrt(x);
#else
#ifdef BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS

#ifdef BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD
  return Base::Math::fast_inverse_square_root(
      x, static_cast<T>(Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN));
#else  // BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD
  return Base::Math::rsqrt_newton_method<T, Base::Math::SQRT_REPEAT_NUMBER>(
      x, static_cast<T>(Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN));
#endif // BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD

#else // BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS
  return Base::Math::rsqrt_newton_method<T, Base::Math::SQRT_REPEAT_NUMBER>(
      x, static_cast<T>(Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN));

#endif // BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS
#endif // BASE_MATH_USE_STD_MATH
}

/* sqrt */
template <typename T, std::size_t LOOP_NUMBER>
inline T sqrt_newton_method(const T &x) {

  T x_wrapped = x;

  if (x < static_cast<T>(Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN)) {
    x_wrapped =
        static_cast<T>(Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN);
  }

  return x_wrapped *
         Base::Math::rsqrt_newton_method<T, LOOP_NUMBER>(
             x_wrapped,
             static_cast<T>(Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN));
}

template <typename T, std::size_t LOOP_NUMBER>
inline T sqrt_newton_method(const T &x, const T &division_min) {

  T x_wrapped = x;

  if (x < division_min) {
    x_wrapped = division_min;
  }

  return x_wrapped * Base::Math::rsqrt_newton_method<T, LOOP_NUMBER>(
                         x_wrapped, division_min);
}

template <typename T> inline T fast_square_root(const T &input) {

  T x_wrapped = input;

  if (input <
      static_cast<T>(Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN)) {
    x_wrapped =
        static_cast<T>(Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN);
  }

  return x_wrapped *
         Base::Math::fast_inverse_square_root(
             x_wrapped,
             static_cast<T>(Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN));
}

template <typename T> inline T sqrt(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::sqrt(x);
#else
#ifdef BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS

#ifdef BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD
  return Base::Math::fast_square_root(x);
#else  // BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD
  return Base::Math::sqrt_newton_method<T, Base::Math::SQRT_REPEAT_NUMBER>(x);
#endif // BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD

#else // BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS
  return Base::Math::sqrt_newton_method<T, Base::Math::SQRT_REPEAT_NUMBER>(x);

#endif // BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS
#endif // BASE_MATH_USE_STD_MATH
}

/* exp Mcloughlin Expansion */
template <typename T, std::size_t LOOP_MAX, std::size_t N>
struct ExpMcloughlinIterationLoop {
  static void compute(T &result, T &term, const T &remainder) {
    term *= remainder / static_cast<T>(LOOP_MAX - N);
    result += term;

    ExpMcloughlinIterationLoop<T, LOOP_MAX, N - 1>::compute(result, term,
                                                            remainder);
  }
};

template <typename T, std::size_t LOOP_MAX>
struct ExpMcloughlinIterationLoop<T, LOOP_MAX, 0> {
  static void compute(T &result, T &term, const T &remainder) {
    /* Do Nothing. */
    static_cast<void>(result);
    static_cast<void>(term);
    static_cast<void>(remainder);
  }
};

template <typename T, std::size_t LOOP_NUMBER>
inline T exp_mcloughlin_expansion(const T &x) {

  T result = static_cast<T>(1);

  if (x > static_cast<T>(Base::Math::EXP_INPUT_MAX)) {
    result = static_cast<T>(Base::Math::EXP_OUTPUT_MAX);
  } else if (x < static_cast<T>(Base::Math::EXP_INPUT_MIN)) {
    result = static_cast<T>(Base::Math::EXP_OUTPUT_MIN);
  } else {

    int n = static_cast<int>(x / static_cast<T>(Base::Math::LN_2));
    T remainder = x - n * static_cast<T>(Base::Math::LN_2);
    T term = static_cast<T>(1);

    ExpMcloughlinIterationLoop<T, LOOP_NUMBER, LOOP_NUMBER - 1>::compute(
        result, term, remainder);

    result = ldexp(result, n);
  }

  return result;
}

template <typename T> inline T exp(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::exp(x);
#else
  return Base::Math::exp_mcloughlin_expansion<T, Base::Math::EXP_REPEAT_NUMBER>(
      x);
#endif
}

/* exp Mcloughlin Expansion with table */
template <std::size_t MCLOUGHLIN_EXPANSION_REPEAT_NUMBER>
double exp_mcloughlin_expansion_with_table(double x) {
  double result = static_cast<double>(1);

  using EXP_MCLOUGHLIN_FACTOR_List = typename MakeExpMcloughlinFactorList<
      MCLOUGHLIN_EXPANSION_REPEAT_NUMBER>::type;

  constexpr auto EXP_MCLOUGHLIN_FACTOR =
      Base::Math::to_exp_mcloughlin_factor_array(EXP_MCLOUGHLIN_FACTOR_List{});

  if (x > static_cast<double>(Base::Math::EXP_INPUT_MAX)) {
    result = static_cast<double>(Base::Math::EXP_OUTPUT_MAX);
  } else if (x < static_cast<double>(Base::Math::EXP_INPUT_MIN)) {
    result = static_cast<double>(Base::Math::EXP_OUTPUT_MIN);
  } else {

    double y = EXP_MCLOUGHLIN_FACTOR[MCLOUGHLIN_EXPANSION_REPEAT_NUMBER - 1];

    double z = static_cast<double>(0);
    double r = static_cast<double>(0);
    int q = static_cast<int>(0);
    unsigned long long w = static_cast<unsigned long long>(0);

    z = x * (static_cast<double>(16) / static_cast<double>(Base::Math::LN_2));
    q = static_cast<int>(z) - (x < static_cast<double>(0));
    r = x -
        ((q << static_cast<int>(1)) + static_cast<int>(1)) *
            (static_cast<double>(Base::Math::LN_2) / static_cast<double>(32));
    w = static_cast<unsigned long long>(static_cast<int>(1023) +
                                        static_cast<int>(q >> 4))
            << static_cast<int>(52) ^
        static_cast<unsigned long long>(
            Base::Math::TABLE_FOR_EXP[q & static_cast<int>(0xF)]);

    std::memcpy(&z, &w, static_cast<int>(8));

    // for (int i = 0; i < static_cast<int>(7); ++i) {
    for (int i = static_cast<int>(0);
         i < static_cast<int>(MCLOUGHLIN_EXPANSION_REPEAT_NUMBER); ++i) {
      y = y * r +
          EXP_MCLOUGHLIN_FACTOR[MCLOUGHLIN_EXPANSION_REPEAT_NUMBER - i - 1];
    }

    result = y * z;
  }

  return result;
}

/* exp2 */
template <typename T, std::size_t LOOP_MAX, std::size_t N>
struct Exp2NewtonIterationLoop {
  static void compute(T &result, T &term, const T &remainder) {
    term *= static_cast<T>(Base::Math::LN_2) * remainder /
            static_cast<T>(LOOP_MAX - N);
    result += term;

    Exp2NewtonIterationLoop<T, LOOP_MAX, N - 1>::compute(result, term,
                                                         remainder);
  }
};

template <typename T, std::size_t LOOP_MAX>
struct Exp2NewtonIterationLoop<T, LOOP_MAX, 0> {
  static void compute(T &result, T &term, const T &remainder) {
    /* Do Nothing. */
    static_cast<void>(result);
    static_cast<void>(term);
    static_cast<void>(remainder);
  }
};

template <typename T, std::size_t LOOP_NUMBER>
inline T exp2_mcloughlin_expansion(const T &x) {

  T result = static_cast<T>(1);

  if (x > static_cast<T>(Base::Math::EXP2_INPUT_MAX)) {
    result = static_cast<T>(Base::Math::EXP2_OUTPUT_MAX);
  } else if (x < static_cast<T>(Base::Math::EXP2_INPUT_MIN)) {
    result = static_cast<T>(Base::Math::EXP2_OUTPUT_MIN);
  } else {

    int n = static_cast<int>(x);
    T remainder = x - static_cast<T>(n);

    T term = static_cast<T>(1);

    Exp2NewtonIterationLoop<T, LOOP_NUMBER, LOOP_NUMBER - 1>::compute(
        result, term, remainder);

    result = ldexp(result, n);
  }

  return result;
}

template <typename T> inline T exp2(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::exp2(x);
#else
  return Base::Math::exp2_mcloughlin_expansion<T,
                                               Base::Math::EXP2_REPEAT_NUMBER>(
      x);
#endif
}

/* log */
template <typename T, std::size_t N> struct LogNewtonIterationLoop {
  static void compute(T &exp_guess, T &guess, const T &scaled_x) {
    exp_guess = Base::Math::exp(guess);
    guess = guess - (exp_guess - scaled_x) / exp_guess;

    LogNewtonIterationLoop<T, N - 1>::compute(exp_guess, guess, scaled_x);
  }
};

template <typename T> struct LogNewtonIterationLoop<T, 0> {
  static void compute(T &exp_guess, T &guess, const T &scaled_x) {
    /* Do Nothing. */
    static_cast<void>(exp_guess);
    static_cast<void>(guess);
    static_cast<void>(scaled_x);
  }
};

template <typename T, std::size_t LOOP_NUMBER>
inline T log_newton_method(const T &x) {
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

    LogNewtonIterationLoop<T, LOOP_NUMBER - 1>::compute(exp_guess, guess,
                                                        scaled_x);

    result = guess + static_cast<T>(scale) *
                         static_cast<T>(Base::Math::LOG_OF_LOG_SCALE_FACTOR);
  }

  return result;
}

template <typename T> inline T log(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::log(x);
#else
  return Base::Math::log_newton_method<T, Base::Math::LOG_REPEAT_NUMBER>(x);
#endif
}

/* log2 */
template <typename T, std::size_t LOOP_NUMBER>
inline T log2_newton_method(const T &x) {

  return Base::Math::log_newton_method<T, LOOP_NUMBER>(x) /
         static_cast<T>(Base::Math::LN_2);
}

template <typename T> inline T log2(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::log2(x);
#else
  return Base::Math::log2_newton_method<T, Base::Math::LOG_REPEAT_NUMBER>(x);
#endif
}

/* log10 */
template <typename T, std::size_t LOOP_NUMBER>
inline T log10_newton_method(const T &x) {

  return Base::Math::log_newton_method<T, LOOP_NUMBER>(x) /
         static_cast<T>(Base::Math::LN_10);
}

template <typename T> inline T log10(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::log10(x);
#else
  return Base::Math::log10_newton_method<T, Base::Math::LOG_REPEAT_NUMBER>(x);
#endif
}

/* pow */
template <typename T, std::size_t EXP_LOOP_NUMBER, std::size_t LOG_LOOP_NUMBER>
inline T pow_base_math(const T &x, const T &y) {
  return Base::Math::exp_mcloughlin_expansion<T, EXP_LOOP_NUMBER>(
      y * Base::Math::log_newton_method<T, LOG_LOOP_NUMBER>(x));
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
