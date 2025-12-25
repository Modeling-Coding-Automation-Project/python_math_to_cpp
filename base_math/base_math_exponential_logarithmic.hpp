/**
 * @file base_math_exponential_logarithmic.hpp
 * @brief High-performance mathematical functions for exponentials, logarithms,
 * and square roots.
 *
 * This header provides a collection of mathematical functions and algorithms
 * for computing exponentials, logarithms, square roots, and related operations.
 * The implementations include fast approximations, Maclaurin expansions,
 * Newton's method, and IEEE 754-dependent algorithms, offering a balance
 * between speed and accuracy. The code is designed to be configurable via
 * macros to select between standard library functions and custom
 * high-performance routines.
 *
 * Main features:
 * - Fast inverse square root (Quake III algorithm)
 * - Square root extraction for float and double using bit manipulation
 * - Exponential and logarithmic functions using Maclaurin expansion and lookup
 * tables
 * - Newton's method for exponentials, logarithms, and square roots
 * - Support for both float and double types
 * - Configurable accuracy and performance via preprocessor macros
 *
 * Notable classes and namespaces:
 * - namespace Base::Math: Contains all mathematical functions and constants.
 * - namespace SqrtExtractionDouble: Helper templates for double-precision
 * square root extraction.
 * - namespace SqrtExtractionFloat: Helper templates for single-precision square
 * root extraction.
 * - namespace RsqrtNewton: Templates for Newton's method in reciprocal square
 * root calculation.
 * - namespace ExpMaclaurinIteration: Templates for Maclaurin expansion of
 * exponential functions.
 * - namespace ExpMcLoughlinExpansion: Templates for Maclaurin expansion with
 * lookup tables.
 * - namespace Exp2NewtonIteration: Templates for Maclaurin expansion of exp2
 * functions.
 * - struct LogNewtonIterationLoop: Template for Newton's method in logarithm
 * calculation.
 */
#ifndef __BASE_MATH_EXPONENTIAL_LOGARITHMIC_HPP__
#define __BASE_MATH_EXPONENTIAL_LOGARITHMIC_HPP__

#include "base_math_macros.hpp"

#include "base_math_arithmetic.hpp"
#include "base_math_mathematical_constants.hpp"
#include "base_math_templates.hpp"

#include <cmath>
#include <cstddef>
#include <cstring>

namespace Base {
namespace Math {

constexpr double EXPONENTIAL_LOGARITHMIC_DIVISION_MIN = 1.0e-10;

#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
constexpr std::size_t SQRT_REPEAT_NUMBER = 0;
constexpr std::size_t EXP_REPEAT_NUMBER = 5;
constexpr std::size_t EXP2_REPEAT_NUMBER = 4;
constexpr std::size_t LOG_REPEAT_NUMBER = 6;

constexpr int SQRT_EXTRACTION_REPEAT_NUMBER = 6;
constexpr std::size_t EXP_MACLAURIN_WITH_TABLE_REPEAT_NUMBER = 1;

#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
constexpr std::size_t SQRT_REPEAT_NUMBER = 1;
constexpr std::size_t EXP_REPEAT_NUMBER = 9;
constexpr std::size_t EXP2_REPEAT_NUMBER = 8;
constexpr std::size_t LOG_REPEAT_NUMBER = 7;

constexpr int SQRT_EXTRACTION_REPEAT_NUMBER = 19;

constexpr std::size_t EXP_MACLAURIN_WITH_TABLE_REPEAT_NUMBER = 4;

#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__

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

constexpr unsigned long long TABLE_FOR_EXP_DOUBLE[16] = {
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

constexpr unsigned long TABLE_FOR_EXP_FLOAT[16] = {
    0x02CD86ul, // 2^( 1 /32)-1
    0x08980Eul, // 2^( 3 /32)-1
    0x0EA439ul, // 2^( 5 /32)-1
    0x14F4EFul, // 2^( 7 /32)-1
    0x1B8D39ul, // 2^( 9 /32)-1
    0x227043ul, // 2^( 11 /32)-1
    0x29A15Aul, // 2^( 13 /32)-1
    0x3123F5ul, // 2^( 15 /32)-1
    0x38FBAFul, // 2^( 17 /32)-1
    0x412C4Cul, // 2^( 19 /32)-1
    0x49B9BDul, // 2^( 21 /32)-1
    0x52A81Dul, // 2^( 23 /32)-1
    0x5BFBB7ul, // 2^( 25 /32)-1
    0x65B906ul, // 2^( 27 /32)-1
    0x6FE4B9ul, // 2^( 29 /32)-1
    0x7A83B2ul, // 2^( 31 /32)-1
};

constexpr std::size_t EXP_MACLAURIN_FACTOR_MAX_SIZE = 7;

using EXP_MACLAURIN_FACTOR_LIST =
    typename ExpMaclaurinFactor::MakeList<EXP_MACLAURIN_FACTOR_MAX_SIZE>::type;

constexpr auto EXP_MACLAURIN_FACTOR =
    ExpMaclaurinFactor::to_array(EXP_MACLAURIN_FACTOR_LIST{});

constexpr double TABLE_FOR_LOG_DOUBLE[17] = {
    0.0,                      // log( 16 /16)
    0.0606246218164348425806, // log( 17 /16)
    0.1177830356563834545387, // log( 18 /16)
    0.17185025692665922234,   // log( 19 /16)
    0.2231435513142097557662, // log( 20 /16)
    0.2719337154836417588316, // log( 21 /16)
    0.3184537311185346158102, // log( 22 /16)
    0.3629054936893684531378, // log( 23 /16)
    0.405465108108164381978,  // log( 24 /16)
    0.4462871026284195115325, // log( 25 /16)
    0.4855078157817008078017, // log( 26 /16)
    0.5232481437645478365168, // log( 27 /16)
    0.5596157879354226862708, // log( 28 /16)
    0.5947071077466927895143, // log( 29 /16)
    0.6286086594223741377443, // log( 30 /16)
    0.6613984822453650082602, // log( 31 /16)
    0.6931471805599453094172, // log( 32 /16)
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
    std::memcpy(&i, &y, static_cast<std::size_t>(4));
    i = static_cast<long>(0x5f3759df) - (i >> 1);
    std::memcpy(&y, &i, static_cast<std::size_t>(4));
    y = y * (static_cast<float>(1.5) - (x * y * y));
  }

  return static_cast<T>(y);
}

/* sqrt extraction double */
namespace SqrtExtractionDouble {

template <int Loop_Limit, int I, int I_Limit> struct FirstLoop {
  /**
   * @brief Executes a single iteration of a specialized mathematical algorithm
   * involving bitwise operations.
   *
   * This static method performs a sequence of bitwise manipulations and
   * arithmetic operations on the provided unsigned long long references. It
   * updates the values of `c`, `m`, `y`, and `a` based on the current state,
   * and then recursively calls the `FirstLoop` template with updated template
   * parameters.
   *
   * @param c Reference to an unsigned long long variable, used as a control or
   * carry flag in the computation.
   * @param m Reference to an unsigned long long variable, typically
   * representing the modulus or a working value.
   * @param y Reference to an unsigned long long variable, used as an
   * accumulator or intermediate result.
   * @param a Reference to an unsigned long long variable, used as an
   * accumulator or result.
   *
   * @tparam Loop_Limit The loop limit for the recursive template call.
   * @tparam I The current bit index or iteration count.
   */
  static void execute(unsigned long long &c, unsigned long long &m,
                      unsigned long long &y, unsigned long long &a) {
    c = static_cast<unsigned long>(
        (y << static_cast<int>(1) | static_cast<unsigned long>(1)) <= (m >> I));
    a = (a << static_cast<int>(1)) | c;
    y = (y << static_cast<int>(1)) | c;
    m -= (c * y) << I;
    y += c;

    FirstLoop<Loop_Limit, (I - 2), (I - 2 - Loop_Limit)>::execute(c, m, y, a);
  }
};

template <int Loop_Limit, int I_Limit>
struct FirstLoop<Loop_Limit, -2, I_Limit> {
  /**
   * @brief Executes the base case operation for the exponential/logarithmic
   * function.
   *
   * This static method serves as a base case and does not perform any
   * operations on the provided parameters. All parameters are marked as unused
   * to avoid compiler warnings.
   *
   * @param c Reference to an unsigned long long value (unused).
   * @param m Reference to an unsigned long long value (unused).
   * @param y Reference to an unsigned long long value (unused).
   * @param a Reference to an unsigned long long value (unused).
   */
  static void execute(unsigned long long &c, unsigned long long &m,
                      unsigned long long &y, unsigned long long &a) {
    static_cast<void>(c);
    static_cast<void>(m);
    static_cast<void>(y);
    static_cast<void>(a);
    // Base case: do nothing
  }
};

/**
 * @brief Specialization of the FirstLoop struct for the case when the third
 * template parameter is -2.
 *
 * This specialization provides a static execute function that takes four
 * unsigned long long references as parameters. In this base case, the function
 * does nothing except explicitly casting each parameter to void to suppress
 * unused variable warnings.
 *
 * @tparam Loop_Limit The loop limit parameter (not used in this
 * specialization).
 * @tparam I The current loop index (not used in this specialization).
 *
 * @param c Reference to an unsigned long long variable (not used).
 * @param m Reference to an unsigned long long variable (not used).
 * @param y Reference to an unsigned long long variable (not used).
 * @param a Reference to an unsigned long long variable (not used).
 */
template <int Loop_Limit, int I> struct FirstLoop<Loop_Limit, I, -2> {
  static void execute(unsigned long long &c, unsigned long long &m,
                      unsigned long long &y, unsigned long long &a) {
    static_cast<void>(c);
    static_cast<void>(m);
    static_cast<void>(y);
    static_cast<void>(a);
    // Base case: do nothing
  }
};

/**
 * @brief Specialization of the FirstLoop template for the base case.
 *
 * This specialization is triggered when the template parameters are
 * <Loop_Limit, -2, -2>. It represents the base case for the recursive
 * FirstLoop structure, where no further action is required.
 *
 * The execute function takes four unsigned long long references as parameters,
 * but does nothing with them, effectively terminating the recursion.
 *
 * @tparam Loop_Limit The loop limit parameter (unused in this specialization).
 * @param c Reference to an unsigned long long variable (unused).
 * @param m Reference to an unsigned long long variable (unused).
 * @param y Reference to an unsigned long long variable (unused).
 * @param a Reference to an unsigned long long variable (unused).
 */
template <int Loop_Limit> struct FirstLoop<Loop_Limit, -2, -2> {
  static void execute(unsigned long long &c, unsigned long long &m,
                      unsigned long long &y, unsigned long long &a) {
    static_cast<void>(c);
    static_cast<void>(m);
    static_cast<void>(y);
    static_cast<void>(a);
    // Base case: do nothing
  }
};

template <int I> struct SecondLoop {
  /**
   * @brief Executes a single iteration of a specialized bitwise algorithm,
   * updating the provided variables.
   *
   * This static method performs a sequence of bitwise and arithmetic operations
   * on the input references. It is typically used as part of a recursive or
   * loop-based algorithm (such as in multi-precision arithmetic, logarithmic or
   * exponential calculations). The function updates the values of `c`, `m`,
   * `y`, and `a` in place, and then recursively calls the next iteration via
   * `SecondLoop<I - 1>::execute`.
   *
   * @param[out] c Reference to an unsigned long long variable, used as a carry
   * or flag in the computation.
   * @param[in,out] m Reference to an unsigned long long variable, typically
   * representing a modulus or accumulator.
   * @param[in,out] y Reference to an unsigned long long variable, often
   * representing a partial result or operand.
   * @param[in,out] a Reference to an unsigned long long variable, often
   * representing an accumulating result.
   */
  static void execute(unsigned long long &c, unsigned long long &m,
                      unsigned long long &y, unsigned long long &a) {
    m <<= static_cast<int>(2);
    c = static_cast<unsigned long long>(
        (y << static_cast<int>(1) | static_cast<unsigned long long>(1)) <= m);
    a = (a << static_cast<int>(1)) | c;
    y = (y << static_cast<int>(1)) | c;
    m -= (c * y);
    y += c;
    SecondLoop<I - 1>::execute(c, m, y, a);
  }
};

template <> struct SecondLoop<0> {
  /**
   * @brief Performs a single iteration of a bitwise algorithm, updating the
   * provided variables.
   *
   * This static function manipulates four unsigned long long references (`c`,
   * `m`, `y`, `a`) using bitwise operations. It is typically used in algorithms
   * involving exponentiation, logarithms, or root extraction in integer
   * arithmetic.
   *
   * @param c Reference to an unsigned long long used as a temporary or carry
   * variable.
   * @param m Reference to an unsigned long long representing the main operand,
   * shifted and reduced.
   * @param y Reference to an unsigned long long, updated with the result of the
   * current iteration.
   * @param a Reference to an unsigned long long, accumulates the result through
   * bitwise operations.
   */
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

template <> struct SecondLoop<-1> {
  /**
   * @brief Executes the base case operation for the exponential/logarithmic
   * function.
   *
   * This static method serves as a base case in a recursive or templated
   * mathematical operation. It takes four unsigned long long references as
   * parameters and performs no operation on them, effectively acting as a
   * no-op. All parameters are explicitly marked as unused to avoid compiler
   * warnings.
   *
   * @param c Reference to an unsigned long long parameter (unused).
   * @param m Reference to an unsigned long long parameter (unused).
   * @param y Reference to an unsigned long long parameter (unused).
   * @param a Reference to an unsigned long long parameter (unused).
   */
  static void execute(unsigned long long &c, unsigned long long &m,
                      unsigned long long &y, unsigned long long &a) {
    static_cast<void>(c);
    static_cast<void>(m);
    static_cast<void>(y);
    static_cast<void>(a);
    // Base case: do nothing
  }
};

} // namespace SqrtExtractionDouble

/**
 * @brief Computes the square root of a double-precision floating-point value
 * using bit manipulation and iterative extraction.
 *
 * This function extracts the exponent and mantissa from the IEEE 754
 * representation of the input value, then reconstructs the square root by
 * adjusting the exponent and performing iterative refinement on the mantissa.
 * The number of extraction/repeat iterations is controlled by the template
 * parameter EXTRACTION_DOUBLE_REPEAT_NUMBER.
 *
 *  When REPEAT_NUMBER is specified as the maximum 26,
 *  it completely matches math.h's sqrt().
 *  When REPEAT_NUMBER is 0, the error is within 8.0e^-7 [%].
 *  When REPEAT_NUMBER is -10, the error is within 8.0e^-3 [%].
 *  When REPEAT_NUMBER is -20, the error is within 8 [%].
 *  At -26, since almost no calculation is performed, the error is very large.
 *
 * @tparam EXTRACTION_DOUBLE_REPEAT_NUMBER Number of extraction iterations for
 * refining the square root calculation.
 * @param value The input double value for which the square root is to be
 * computed.
 * @return The computed square root of the input value.
 *
 * @note This implementation uses low-level bit manipulation and custom
 * iterative loops for mantissa extraction. It is intended for specialized use
 * cases where standard library sqrt may not be suitable. The function assumes
 * IEEE 754 double-precision representation.
 */
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

  SqrtExtractionDouble::FirstLoop<
      ((static_cast<int>(-2) * EXTRACTION_DOUBLE_REPEAT_NUMBER -
        static_cast<int>(2)) *
       (EXTRACTION_DOUBLE_REPEAT_NUMBER < static_cast<int>(-1))),
      54,
      (54 - ((static_cast<int>(-2) * EXTRACTION_DOUBLE_REPEAT_NUMBER -
              static_cast<int>(2)) *
             (EXTRACTION_DOUBLE_REPEAT_NUMBER <
              static_cast<int>(-1))))>::execute(c, m, y, a);

  SqrtExtractionDouble::SecondLoop<(
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
namespace SqrtExtractionFloat {

template <int Loop_Limit, int I, int I_Limit> struct FirstLoop {
  /**
   * @brief Executes a single iteration of a specialized bitwise algorithm.
   *
   * This static method performs a sequence of bitwise and arithmetic operations
   * on the provided unsigned long references. It manipulates the values of `c`,
   * `m`, `y`, and `a` according to a specific algorithm, likely related to
   * exponential or logarithmic computations. The method also recursively
   * invokes the `FirstLoop` template with updated template parameters.
   *
   * @param c Reference to an unsigned long variable, used as a bitwise flag and
   * updated in-place.
   * @param m Reference to an unsigned long variable, representing a mutable
   * state in the algorithm.
   * @param y Reference to an unsigned long variable, updated through bitwise
   * operations.
   * @param a Reference to an unsigned long variable, accumulating results
   * through bitwise operations.
   */
  static void execute(unsigned long &c, unsigned long &m, unsigned long &y,
                      unsigned long &a) {
    c = static_cast<unsigned long>(
        (y << static_cast<int>(1) | static_cast<unsigned long>(1)) <= (m >> I));
    a = (a << static_cast<int>(1)) | c;
    y = (y << static_cast<int>(1)) | c;
    m -= (c * y) << I;
    y += c;

    FirstLoop<Loop_Limit, (I - 2), (I - 2 - Loop_Limit)>::execute(c, m, y, a);
  }
};

template <int Loop_Limit, int I_Limit>
struct FirstLoop<Loop_Limit, -1, I_Limit> {

  /**
   * @brief Executes an operation on the provided unsigned long references.
   *
   * This is a base case implementation that performs no operation on the input
   * parameters. All parameters are explicitly marked as unused to avoid
   * compiler warnings.
   *
   * @param c Reference to an unsigned long value (unused).
   * @param m Reference to an unsigned long value (unused).
   * @param y Reference to an unsigned long value (unused).
   * @param a Reference to an unsigned long value (unused).
   */
  static void execute(unsigned long &c, unsigned long &m, unsigned long &y,
                      unsigned long &a) {
    static_cast<void>(c);
    static_cast<void>(m);
    static_cast<void>(y);
    static_cast<void>(a);
    // Base case: do nothing
  }
};

template <int Loop_Limit, int I> struct FirstLoop<Loop_Limit, I, -1> {
  /**
   * @brief Executes an operation on the provided unsigned long references.
   *
   * This is a base case implementation that performs no operation on the input
   * parameters. All parameters are explicitly marked as unused to avoid
   * compiler warnings.
   *
   * @param c Reference to an unsigned long value (unused).
   * @param m Reference to an unsigned long value (unused).
   * @param y Reference to an unsigned long value (unused).
   * @param a Reference to an unsigned long value (unused).
   */
  static void execute(unsigned long &c, unsigned long &m, unsigned long &y,
                      unsigned long &a) {
    static_cast<void>(c);
    static_cast<void>(m);
    static_cast<void>(y);
    static_cast<void>(a);
    // Base case: do nothing
  }
};

template <int Loop_Limit> struct FirstLoop<Loop_Limit, -1, -1> {
  /**
   * @brief Executes an operation on the provided unsigned long references.
   *
   * This is a base case implementation that performs no operation on the input
   * parameters. All parameters are explicitly marked as unused to avoid
   * compiler warnings.
   *
   * @param c Reference to an unsigned long value (unused).
   * @param m Reference to an unsigned long value (unused).
   * @param y Reference to an unsigned long value (unused).
   * @param a Reference to an unsigned long value (unused).
   */
  static void execute(unsigned long &c, unsigned long &m, unsigned long &y,
                      unsigned long &a) {
    static_cast<void>(c);
    static_cast<void>(m);
    static_cast<void>(y);
    static_cast<void>(a);
    // Base case: do nothing
  }
};

template <int I> struct SecondLoop {
  /**
   * @brief Executes a single iteration of the SecondLoop algorithm.
   *
   * This static method performs bitwise and arithmetic operations on the
   * provided unsigned long references, updating their values according to the
   * algorithm's logic. It shifts and manipulates the variables `m`, `c`, `y`,
   * and `a`, and recursively calls the next iteration of SecondLoop.
   *
   * @param c Reference to an unsigned long used as a control or carry variable.
   * @param m Reference to an unsigned long, typically representing a modulus or
   * intermediate value.
   * @param y Reference to an unsigned long, used as an accumulator or
   * intermediate result.
   * @param a Reference to an unsigned long, used as an accumulator or result
   * variable.
   */
  static void execute(unsigned long &c, unsigned long &m, unsigned long &y,
                      unsigned long &a) {
    m <<= static_cast<int>(2);
    c = static_cast<unsigned long>(
        (y << static_cast<int>(1) | static_cast<unsigned long>(1)) <= m);
    a = (a << static_cast<int>(1)) | c;
    y = (y << static_cast<int>(1)) | c;
    m -= (c * y);
    y += c;
    SecondLoop<I - 1>::execute(c, m, y, a);
  }
};

template <> struct SecondLoop<0> {
  /**
   * @brief Performs a single step of a bitwise algorithm, updating the provided
   * variables.
   *
   * This static function manipulates the input references using bitwise
   * operations and arithmetic. It is likely used in algorithms involving
   * exponentiation, logarithms, or similar mathematical computations where
   * bitwise manipulation is required.
   *
   * @param c Reference to an unsigned long variable, used as a flag or result
   * of a comparison.
   * @param m Reference to an unsigned long variable, shifted and decremented
   * based on computation.
   * @param y Reference to an unsigned long variable, updated with bitwise
   * operations and incremented.
   * @param a Reference to an unsigned long variable, updated with bitwise
   * operations.
   */
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

template <> struct SecondLoop<-1> {
  /**
   * @brief Executes the base case operation for the exponential/logarithmic
   * function.
   *
   * This static method is intended as a base case and performs no operation on
   * the provided unsigned long references. All parameters are explicitly marked
   * as unused to avoid compiler warnings.
   *
   * @param c Reference to an unsigned long value (unused).
   * @param m Reference to an unsigned long value (unused).
   * @param y Reference to an unsigned long value (unused).
   * @param a Reference to an unsigned long value (unused).
   */
  static void execute(unsigned long &c, unsigned long &m, unsigned long &y,
                      unsigned long &a) {
    static_cast<void>(c);
    static_cast<void>(m);
    static_cast<void>(y);
    static_cast<void>(a);
    // Base case: do nothing
  }
};

} // namespace SqrtExtractionFloat

/**
 * @brief Computes the square root of a floating-point value using bit
 * manipulation and iterative extraction.
 *
 * This function implements a custom square root calculation for 32-bit
 * floating-point numbers (float) by decomposing the input value into its
 * exponent and mantissa, manipulating the bits directly, and performing
 * iterative refinement using template metaprogramming for the number of
 * extraction steps.
 *
 * @tparam EXTRACTION_FLOAT_REPEAT_NUMBER Number of iterations for the
 * extraction loop, controlling the precision.
 * @param value The input floating-point value to compute the square root of.
 * @return The computed square root of the input value as a float.
 *
 * @note This function is intended for specialized use cases where standard
 * library sqrt may not be suitable, such as environments requiring bit-level
 * control or custom precision/iteration management.
 * @warning The function assumes IEEE 754 single-precision float representation
 * and may not handle special cases (e.g., negative values, NaN, infinity) as
 * the standard sqrt function does.
 */
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
  SqrtExtractionFloat::FirstLoop<
      ((static_cast<int>(-2) * EXTRACTION_FLOAT_REPEAT_NUMBER -
        static_cast<int>(2)) *
       (EXTRACTION_FLOAT_REPEAT_NUMBER < static_cast<int>(-1))),
      25,
      (25 - ((static_cast<int>(-2) * EXTRACTION_FLOAT_REPEAT_NUMBER -
              static_cast<int>(2)) *
             (EXTRACTION_FLOAT_REPEAT_NUMBER <
              static_cast<int>(-1))))>::execute(c, m, y, a);

  y <<= static_cast<int>(1);

  SqrtExtractionFloat::SecondLoop<(
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

/* sqrt extraction */
/**
 * @brief Computes the square root extraction for floating point values of type
 * double, using a specified number of extraction repetitions.
 *
 * This function template is enabled only for the type double. It delegates the
 * computation to Base::Math::sqrt_extraction_double, passing the number of
 * repetitions as (EXTRACTION_FloatingPoint_REPEAT_NUMBER - 26).
 *
 * @tparam T The floating point type (must be double for this specialization).
 * @tparam EXTRACTION_FloatingPoint_REPEAT_NUMBER The number of extraction
 * repetitions to use.
 * @param x The input value for which the square root extraction is computed.
 * @return The result of the square root extraction as a value of type T.
 */
template <typename T, int EXTRACTION_FloatingPoint_REPEAT_NUMBER>
typename std::enable_if<std::is_same<T, double>::value, T>::type
sqrt_extraction_FloatingPoint(const T &x) {

  return Base::Math::sqrt_extraction_double<
      static_cast<int>(EXTRACTION_FloatingPoint_REPEAT_NUMBER) -
      static_cast<int>(26)>(x);
}

/**
 * @brief Computes the square root extraction for floating point types (float
 * specialization).
 *
 * This template function is enabled only when the type T is float. It calls the
 * Base::Math::sqrt_extraction_float function with a template parameter
 * calculated as (EXTRACTION_FloatingPoint_REPEAT_NUMBER - 12), passing the
 * input value x.
 *
 * @tparam T The floating point type (must be float for this specialization).
 * @tparam EXTRACTION_FloatingPoint_REPEAT_NUMBER The repeat number used for
 * extraction, must be >= 12.
 * @param x The input value for which the square root extraction is performed.
 * @return The result of the square root extraction as a value of type T.
 */
template <typename T, int EXTRACTION_FloatingPoint_REPEAT_NUMBER>
typename std::enable_if<std::is_same<T, float>::value, T>::type
sqrt_extraction_FloatingPoint(const T &x) {

  return Base::Math::sqrt_extraction_float<
      static_cast<int>(EXTRACTION_FloatingPoint_REPEAT_NUMBER) -
      static_cast<int>(12)>(x);
}

/**
 * @brief Computes the square root of a value using a custom extraction method.
 *
 * This function template delegates the computation to
 * Base::Math::sqrt_extraction_FloatingPoint, allowing for a configurable
 * number of extraction iterations via the
 * EXTRACTION_FloatingPoint_REPEAT_NUMBER template parameter.
 *
 * @tparam T The numeric type of the input value (e.g., float, double).
 * @tparam EXTRACTION_FloatingPoint_REPEAT_NUMBER The number of extraction
 * iterations to perform.
 * @param x The value whose square root is to be computed.
 * @return The computed square root of x.
 */
template <typename T, int EXTRACTION_FloatingPoint_REPEAT_NUMBER>
inline T sqrt_extraction(const T &x) {

  return Base::Math::sqrt_extraction_FloatingPoint<
      T, EXTRACTION_FloatingPoint_REPEAT_NUMBER>(x);
}

/* rsqrt newton method loop */
namespace RsqrtNewton {

template <typename T, int N> struct Loop {
  /**
   * @brief Performs a single iteration of the Newton-Raphson method to refine
   * the approximation of the reciprocal square root of x_T.
   *
   * This static function updates the value of 'result' in place using the
   * Newton-Raphson formula: result = result * (1.5 - x_T * result * result) It
   * then recursively calls the next iteration via Loop<T, N - 1>::compute.
   *
   * @tparam T The numeric type (e.g., float, double).
   * @tparam N The number of iterations (handled by the Loop template).
   * @param x_T The input value for which the reciprocal square root is being
   * approximated.
   * @param result Reference to the current approximation, which will be
   * updated.
   */
  static void compute(T x_T, T &result) {
    result *= static_cast<T>(1.5) - x_T * result * result;
    Loop<T, N - 1>::compute(x_T, result);
  }
};

template <typename T> struct Loop<T, 1> {
  /**
   * @brief Performs a single iteration of the Newton-Raphson method to
   * approximate the reciprocal square root of x_T.
   *
   * This function updates the value of 'result' using the Newton-Raphson
   * formula: result = result * (1.5 - x_T * result * result) which improves the
   * approximation of 1 / sqrt(x_T).
   *
   * @tparam T Numeric type (e.g., float, double).
   * @param x_T The input value for which the reciprocal square root is being
   * approximated.
   * @param result Reference to the current approximation of the reciprocal
   * square root; will be updated in place.
   */
  static void compute(T x_T, T &result) {
    result *= static_cast<T>(1.5) - x_T * result * result;
  }
};

template <typename T> struct Loop<T, 0> {
  /**
   * @brief Computes a value based on the input x_T and stores the result in
   * result.
   *
   * This is a placeholder implementation that currently does nothing.
   *
   * @tparam T The data type of the input and output.
   * @param x_T The input value.
   * @param result Reference to the variable where the computed result should be
   * stored.
   */
  static void compute(T x_T, T &result) {
    /* Do Nothing. */
    static_cast<void>(x_T);
    static_cast<void>(result);
  }
};

} // namespace RsqrtNewton

/* rsqrt newton method */
/**
 * @brief Computes the reciprocal square root of a value using the
 * Newton-Raphson method.
 *
 * This function estimates 1/sqrt(x) for the input value `x` using an initial
 * approximation and a specified number of Newton-Raphson refinement iterations.
 * The function ensures numerical stability by clamping `x` to a minimum value
 * (`division_min`) if necessary.
 *
 * @tparam T           The floating-point type of the input and output (e.g.,
 * float, double).
 * @tparam LOOP_NUMBER The number of Newton-Raphson refinement iterations to
 * perform.
 * @param x            The input value for which to compute the reciprocal
 * square root.
 * @param division_min The minimum value to use for division to avoid
 * instability.
 * @return T           The estimated reciprocal square root of `x`.
 */
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
  float r = 1.8284271F - 0.82842712F * fast_frexpf(x_float, &e);

  r = fast_ldexp(
      r * Base::Math::ONE_AND_SQRT2_VEC[e & static_cast<int>(0x00000001)],
      -e >> 1);

  result = static_cast<T>(r);

  h = x_wrapped * result * result;
  result *= static_cast<T>(1.875) -
            h * (static_cast<T>(1.25) - h * static_cast<T>(0.375));

  RsqrtNewton::Loop<T, LOOP_NUMBER>::compute(x_wrapped * static_cast<T>(0.5),
                                             result);

  return result;
}

/* rsqrt */
/**
 * @brief Computes the reciprocal square root of a value with division-by-zero
 * protection.
 *
 * This function calculates 1 / sqrt(x) for the given input value `x`, using
 * different algorithms depending on compile-time configuration macros:
 * - If `__BASE_MATH_USE_STD_MATH__` is defined, it uses the standard library's
 * `std::sqrt`.
 * - If `__BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD__` is defined,
 * it uses a custom square root extraction algorithm.
 * - Otherwise, it uses the Newton-Raphson method for reciprocal square root
 * calculation.
 *
 * If `x` is less than `division_min`, the function uses `division_min` instead
 * to avoid division by zero or taking the square root of a non-positive number.
 *
 * @tparam T Numeric type of the input and output values.
 * @param x The input value for which to compute the reciprocal square root.
 * @param division_min The minimum value to use for `x` to prevent division by
 * zero.
 * @return The reciprocal square root of `x` (or `division_min` if `x <
 * division_min`).
 */
template <typename T> inline T rsqrt(const T &x, const T &division_min) {

#ifdef __BASE_MATH_USE_STD_MATH__

  T x_wrapped = x;

  if (x < division_min) {
    x_wrapped = division_min;
  }

  return static_cast<T>(1) / std::sqrt(x_wrapped);

#else // __BASE_MATH_USE_STD_MATH__

#ifdef __BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD__

  T x_wrapped = x;

  if (x < division_min) {
    x_wrapped = division_min;
  }

  return static_cast<T>(1) /
         Base::Math::sqrt_extraction<T,
                                     Base::Math::SQRT_EXTRACTION_REPEAT_NUMBER>(
             x_wrapped);

#else // __BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD__

  return Base::Math::rsqrt_newton_method<T, Base::Math::SQRT_REPEAT_NUMBER>(
      x, division_min);

#endif // __BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD__

#endif // __BASE_MATH_USE_STD_MATH__
}

/**
 * @brief Computes the reciprocal of the square root of the input value.
 *
 * This function returns 1 / sqrt(x) for the given input x, using the
 * Base::Math::rsqrt implementation. It also applies a minimum division
 * threshold defined by Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN
 * to ensure numerical stability.
 *
 * @tparam T Numeric type of the input value.
 * @param x The value for which to compute the reciprocal square root.
 * @return The reciprocal of the square root of x.
 */
template <typename T> inline T rsqrt(const T &x) {

  return Base::Math::rsqrt<T>(
      x, static_cast<T>(Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN));
}

/* sqrt newton method */
/**
 * @brief Computes the square root of a value using the Newton-Raphson method.
 *
 * This function wraps the input value to ensure it is not less than a
 * predefined minimum (Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN). It
 * then computes the square root by multiplying the wrapped value by the result
 * of the reciprocal square root calculation using Newton's method.
 *
 * @tparam T The numeric type of the input value.
 * @tparam LOOP_NUMBER The number of iterations for the Newton-Raphson method.
 * @param x The input value for which the square root is to be computed.
 * @return The computed square root of the input value.
 */
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

/**
 *  When REPEAT_NUMBER is specified as the maximum 26,
 *  it completely matches math.h's sqrt().
 *  When REPEAT_NUMBER is 0, the error is within 8.0e^-7 [%].
 *  When REPEAT_NUMBER is -10, the error is within 8.0e^-3 [%].
 *  When REPEAT_NUMBER is -20, the error is within 8 [%].
 *  At -26, since almost no calculation is performed, the error is very large.
 */
template <typename T, std::size_t LOOP_NUMBER>
inline T sqrt_newton_method(const T &x, const T &division_min) {

  T x_wrapped = x;

  if (x < division_min) {
    x_wrapped = division_min;
  }

  return x_wrapped * Base::Math::rsqrt_newton_method<T, LOOP_NUMBER>(
                         x_wrapped, division_min);
}

/**
 * @brief Computes an approximate square root of the input value using a fast
 * inverse square root algorithm.
 *
 * This function wraps the input value to ensure it is not less than a
 * predefined minimum threshold
 * (`Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN`). If the input is below
 * this threshold, it is clamped to the minimum value. The function then
 * calculates the square root by multiplying the (possibly clamped) input by the
 * result of `Base::Math::fast_inverse_square_root`, which is expected to
 * provide a fast approximation of the inverse square root.
 *
 * @tparam T Numeric type of the input value.
 * @param input The value for which to compute the square root.
 * @return Approximate square root of the input value.
 */
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

/* sqrt */
/**
 * @brief Computes the square root of a given value with a minimum division
 * threshold.
 *
 * This function calculates the square root of the input value `x` using
 * different algorithms depending on compile-time macros:
 * - If `__BASE_MATH_USE_STD_MATH__` is defined, it uses `std::sqrt` from the
 * standard library.
 * - If `__BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD__` is defined,
 * it uses a custom IEEE 754-dependent algorithm
 * (`Base::Math::sqrt_extraction`).
 * - Otherwise, it uses the Newton-Raphson method
 * (`Base::Math::sqrt_newton_method`).
 *
 * If the input value `x` is less than `division_min`, the function uses
 * `division_min` instead to avoid division by very small numbers or zero.
 *
 * @tparam T Numeric type of the input and output values.
 * @param x The value whose square root is to be computed.
 * @param division_min The minimum threshold for the input value to prevent
 * division by small numbers.
 * @return The computed square root of `x`, or of `division_min` if `x <
 * division_min`.
 */
template <typename T> inline T sqrt(const T &x, const T &division_min) {

#ifdef __BASE_MATH_USE_STD_MATH__

  T x_wrapped = x;

  if (x < division_min) {
    x_wrapped = division_min;
  }

  return std::sqrt(x_wrapped);

#else // __BASE_MATH_USE_STD_MATH__

#ifdef __BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD__

  T x_wrapped = x;

  if (x < division_min) {
    x_wrapped = division_min;
  }

  return Base::Math::sqrt_extraction<T,
                                     Base::Math::SQRT_EXTRACTION_REPEAT_NUMBER>(
      x_wrapped);

#else // __BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD__

  return Base::Math::sqrt_newton_method<T, Base::Math::SQRT_REPEAT_NUMBER>(
      x, division_min);

#endif // __BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD__

#endif // __BASE_MATH_USE_STD_MATH__
}

/**
 * @brief Computes the square root of the given value.
 *
 * This function is a wrapper around Base::Math::sqrt<T>, providing
 * a convenient interface for calculating the square root of a value of type T.
 * It uses a predefined minimum division constant for exponential and
 * logarithmic operations.
 *
 * @tparam T The numeric type of the input value.
 * @param x The value whose square root is to be computed.
 * @return The square root of the input value x.
 */
template <typename T> inline T sqrt(const T &x) {

  return Base::Math::sqrt<T>(
      x, static_cast<T>(Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN));
}

/* exp Maclaurin Expansion */
namespace ExpMaclaurinIteration {

/**
 * @brief Recursive template struct to perform a compile-time loop for
 * mathematical computations.
 *
 * This struct defines a static compute function that recursively updates the
 * result and term using the provided remainder. The recursion is controlled by
 * the template parameter N, decrementing it at each step until the base case is
 * reached (typically specialized for N = 0).
 *
 * @tparam T        The numeric type used for computation (e.g., float, double).
 * @tparam LOOP_MAX The maximum number of loop iterations (used for
 * normalization).
 * @tparam N        The current iteration index (recursively decremented).
 *
 * @param result    Reference to the accumulated result of the computation.
 * @param term      Reference to the current term in the computation, updated
 * each iteration.
 * @param remainder The value used to update the term at each step.
 */
template <typename T, std::size_t LOOP_MAX, std::size_t N> struct Loop {
  static void compute(T &result, T &term, const T &remainder) {
    term *= remainder / static_cast<T>(LOOP_MAX - N);
    result += term;

    Loop<T, LOOP_MAX, N - 1>::compute(result, term, remainder);
  }
};

/**
 * @brief Specialization of the Loop struct for the base case where the loop
 * counter is 0.
 *
 * This specialization provides a static compute function that effectively does
 * nothing, serving as the termination point for recursive template
 * instantiations of Loop. All parameters are marked as unused to avoid compiler
 * warnings.
 *
 * @tparam T The type of the values being processed.
 * @tparam LOOP_MAX The maximum number of loop iterations (not used in this
 * specialization).
 */
template <typename T, std::size_t LOOP_MAX> struct Loop<T, LOOP_MAX, 0> {
  static void compute(T &result, T &term, const T &remainder) {
    /* Do Nothing. */
    static_cast<void>(result);
    static_cast<void>(term);
    static_cast<void>(remainder);
  }
};

} // namespace ExpMaclaurinIteration

/**
 * @brief Computes the exponential function e^x using the Maclaurin series
 * expansion.
 *
 * This function approximates the exponential of a given value `x` using a
 * Maclaurin series expansion up to `LOOP_NUMBER` terms. It handles large and
 * small values of `x` by clamping the result to predefined maximum and minimum
 * output values.
 *
 * The computation splits `x` into an integer multiple of ln(2) and a remainder,
 * computes the Maclaurin series for the remainder, and then scales the result
 * by 2^n using ldexp for improved numerical stability.
 *
 * @tparam T The floating-point type to use for computation (e.g., float,
 * double).
 * @tparam LOOP_NUMBER The number of terms to use in the Maclaurin series
 * expansion.
 * @param x The exponent value for which to compute e^x.
 * @return The approximate value of e^x.
 */
template <typename T, std::size_t LOOP_NUMBER>
inline T exp_maclaurin_expansion(const T &x) {

  T result = static_cast<T>(1);

  if (x > static_cast<T>(Base::Math::EXP_INPUT_MAX)) {
    result = static_cast<T>(Base::Math::EXP_OUTPUT_MAX);
  } else if (x < static_cast<T>(Base::Math::EXP_INPUT_MIN)) {
    result = static_cast<T>(Base::Math::EXP_OUTPUT_MIN);
  } else {

    int n = static_cast<int>(x / static_cast<T>(Base::Math::LN_2));
    T remainder = x - n * static_cast<T>(Base::Math::LN_2);
    T term = static_cast<T>(1);

    ExpMaclaurinIteration::Loop<T, LOOP_NUMBER, LOOP_NUMBER - 1>::compute(
        result, term, remainder);

    result = fast_ldexp(result, n);
  }

  return result;
}

/* exp Maclaurin Expansion with table */
namespace ExpMcLoughlinExpansion {

/**
 * @brief Recursive template structure to perform a looped computation for
 * Maclaurin series expansion.
 *
 * This struct recursively computes a value using the Maclaurin series
 * coefficients for the exponential function. At each recursion step, it updates
 * the value of `y` by multiplying it with `r` and adding the corresponding
 * coefficient from `Base::Math::EXP_MACLAURIN_FACTOR`. The recursion continues
 * until the base case is reached.
 *
 * @tparam T The numeric type used for computation (e.g., float, double).
 * @tparam N The current step in the recursion, corresponding to the Maclaurin
 * series term.
 *
 * @param y Reference to the value being computed and updated at each step.
 * @param r The value used in the multiplication at each recursion step.
 */
template <typename T, std::size_t N> struct Loop {
  static void compute(T &y, T &r) {
    y = y * r + static_cast<T>(Base::Math::EXP_MACLAURIN_FACTOR[N - 1]);
    Loop<T, N - 1>::compute(y, r);
  }
};

template <typename T> struct Loop<T, 0> {
  /**
   * @brief Computes a result based on the input parameters y and r.
   *
   * This static function is intended to perform a computation using the
   * references y and r. In the current implementation, the function does
   * nothing and simply suppresses unused variable warnings.
   *
   * @tparam T The type of the input parameters.
   * @param y Reference to the first input parameter.
   * @param r Reference to the second input parameter.
   */
  static void compute(T &y, T &r) {
    static_cast<void>(y);
    static_cast<void>(r);
    // Do nothing
  }
};

} // namespace ExpMcLoughlinExpansion

/**
 * @brief Computes the exponential function exp(x) using a Maclaurin series
 * expansion with precomputed tables for improved performance and accuracy.
 *
 * This function calculates exp(x) for a given double-precision input `x` using
 * a Maclaurin series expansion. The number of terms in the expansion is
 * specified by the template parameter `MACLAURIN_EXPANSION_REPEAT_NUMBER`. For
 * efficiency, the function uses precomputed tables and bit manipulation to
 * approximate the result, handling edge cases where `x` is outside the
 * representable range.
 *
 * @tparam MACLAURIN_EXPANSION_REPEAT_NUMBER Number of terms to use in the
 * Maclaurin series expansion. Must not exceed
 * `Base::Math::EXP_MACLAURIN_FACTOR_MAX_SIZE`.
 * @param x The input value for which to compute exp(x).
 * @return The computed value of exp(x) as a double.
 *
 * @note If `x` is greater than `Base::Math::EXP_INPUT_MAX`, the function
 * returns `Base::Math::EXP_OUTPUT_MAX`. If `x` is less than
 * `Base::Math::EXP_INPUT_MIN`, the function returns
 * `Base::Math::EXP_OUTPUT_MIN`. Otherwise, the function uses a combination of
 * table lookups, bit manipulation, and Maclaurin expansion.
 *
 * @see Base::Math::EXP_MACLAURIN_FACTOR
 * @see Base::Math::TABLE_FOR_EXP_DOUBLE
 * @see ExpMcLoughlinExpansion::Loop
 */
template <std::size_t MACLAURIN_EXPANSION_REPEAT_NUMBER>
inline double exp_double_maclaurin_expansion_with_table(const double &x) {
  static_assert(MACLAURIN_EXPANSION_REPEAT_NUMBER <=
                    Base::Math::EXP_MACLAURIN_FACTOR_MAX_SIZE,
                "The number of iterations is too large.");

  double result = static_cast<double>(1);

  if (x > static_cast<double>(Base::Math::EXP_INPUT_MAX)) {
    result = static_cast<double>(Base::Math::EXP_OUTPUT_MAX);
  } else if (x < static_cast<double>(Base::Math::EXP_INPUT_MIN)) {
    result = static_cast<double>(Base::Math::EXP_OUTPUT_MIN);
  } else {

    double y =
        Base::Math::EXP_MACLAURIN_FACTOR[MACLAURIN_EXPANSION_REPEAT_NUMBER - 1];

    double z = static_cast<double>(0);
    double r = static_cast<double>(0);
    int q = static_cast<int>(0);
    unsigned long long w = static_cast<unsigned long long>(0);

    z = x * (static_cast<double>(16) / static_cast<double>(Base::Math::LN_2));
    q = static_cast<int>(z) - (x < static_cast<double>(0));
    r = x -
        ((q << static_cast<int>(1)) + static_cast<int>(1)) *
            (static_cast<double>(Base::Math::LN_2) / static_cast<double>(32));
    w = static_cast<unsigned long long>(
            static_cast<int>(1023) + static_cast<int>(q >> static_cast<int>(4)))
            << static_cast<int>(52) ^
        static_cast<unsigned long long>(
            Base::Math::TABLE_FOR_EXP_DOUBLE[q & static_cast<int>(0xF)]);

    std::memcpy(&z, &w, static_cast<int>(8));

    ExpMcLoughlinExpansion::Loop<double,
                                 MACLAURIN_EXPANSION_REPEAT_NUMBER>::compute(y,
                                                                             r);

    result = y * z;
  }

  return result;
}

/**
 * @brief Computes the exponential of a floating-point value using a Maclaurin
 * series expansion with table-based optimizations.
 *
 * This function calculates exp(x) for a given float x using a Maclaurin series
 * expansion, with the number of terms specified by the template parameter
 * MACLAURIN_EXPANSION_REPEAT_NUMBER. It utilizes precomputed tables for
 * improved performance and handles input clamping for extreme values.
 *
 * @tparam MACLAURIN_EXPANSION_REPEAT_NUMBER Number of terms to use in the
 * Maclaurin series expansion. Must not exceed
 * Base::Math::EXP_MACLAURIN_FACTOR_MAX_SIZE.
 * @param x The input value for which to compute the exponential.
 * @return The computed exponential value as a float.
 *
 * @note For inputs greater than Base::Math::EXP_INPUT_MAX, the result is
 * clamped to Base::Math::EXP_OUTPUT_MAX. For inputs less than
 * Base::Math::EXP_INPUT_MIN, the result is clamped to
 * Base::Math::EXP_OUTPUT_MIN.
 * @note Relies on Base::Math::EXP_MACLAURIN_FACTOR,
 * Base::Math::TABLE_FOR_EXP_FLOAT, and other constants for computation.
 * @note Uses ExpMcLoughlinExpansion::Loop for the Maclaurin expansion loop.
 */
template <std::size_t MACLAURIN_EXPANSION_REPEAT_NUMBER>
inline float exp_float_maclaurin_expansion_with_table(const float &x) {
  static_assert(MACLAURIN_EXPANSION_REPEAT_NUMBER <=
                    Base::Math::EXP_MACLAURIN_FACTOR_MAX_SIZE,
                "The number of iterations is too large.");

  float result = static_cast<float>(1);

  if (x > static_cast<float>(Base::Math::EXP_INPUT_MAX)) {
    result = static_cast<float>(Base::Math::EXP_OUTPUT_MAX);
  } else if (x < static_cast<float>(Base::Math::EXP_INPUT_MIN)) {
    result = static_cast<float>(Base::Math::EXP_OUTPUT_MIN);
  } else {

    float y = static_cast<float>(
        Base::Math::EXP_MACLAURIN_FACTOR[MACLAURIN_EXPANSION_REPEAT_NUMBER -
                                         1]);

    float z = static_cast<float>(0);
    float r = static_cast<float>(0);
    int q = static_cast<int>(0);
    unsigned long w = static_cast<unsigned long>(0);

    z = x * (static_cast<float>(16) / static_cast<float>(Base::Math::LN_2));
    q = static_cast<int>(z) - (x < static_cast<float>(0));
    r = x - ((q << static_cast<int>(1)) + static_cast<int>(1)) *
                (static_cast<float>(Base::Math::LN_2) / static_cast<float>(32));
    w = static_cast<unsigned long>(static_cast<int>(127) +
                                   static_cast<int>(q >> static_cast<int>(4)))
            << static_cast<int>(23) ^
        static_cast<unsigned long>(
            Base::Math::TABLE_FOR_EXP_FLOAT[q & static_cast<int>(0xF)]);

    std::memcpy(&z, &w, static_cast<std::size_t>(4));

    ExpMcLoughlinExpansion::Loop<float,
                                 MACLAURIN_EXPANSION_REPEAT_NUMBER>::compute(y,
                                                                             r);

    result = y * z;
  }

  return result;
}

/**
 * @brief Computes the exponential of a floating-point value using the Maclaurin
 * series expansion with a precomputed table for optimization.
 *
 * This function template is enabled only for the `double` type. It delegates
 * the computation to `Base::Math::exp_double_maclaurin_expansion_with_table`,
 * which uses a Maclaurin series expansion repeated
 * `MACLAURIN_EXPANSION_REPEAT_NUMBER` times for improved accuracy and
 * performance.
 *
 * @tparam T The floating-point type (must be `double`).
 * @tparam MACLAURIN_EXPANSION_REPEAT_NUMBER The number of terms to use in the
 * Maclaurin expansion.
 * @param x The input value for which to compute the exponential.
 * @return The exponential of `x` computed using the Maclaurin series expansion.
 */
template <typename T, std::size_t MACLAURIN_EXPANSION_REPEAT_NUMBER>
typename std::enable_if<std::is_same<T, double>::value, T>::type
exp_FloatingPoint_maclaurin_expansion_with_table(const T &x) {

  return Base::Math::exp_double_maclaurin_expansion_with_table<
      MACLAURIN_EXPANSION_REPEAT_NUMBER>(x);
}

/**
 * @brief Computes the exponential of a floating-point value using the Maclaurin
 * series expansion with a lookup table.
 *
 * This function template is enabled only for the `float` type and delegates the
 * computation to `Base::Math::exp_float_maclaurin_expansion_with_table`,
 * parameterized by the specified number of Maclaurin expansion terms
 * (`MACLAURIN_EXPANSION_REPEAT_NUMBER`).
 *
 * @tparam T The floating-point type (must be `float`).
 * @tparam MACLAURIN_EXPANSION_REPEAT_NUMBER The number of terms to use in the
 * Maclaurin expansion.
 * @param x The input value for which to compute the exponential.
 * @return The exponential of `x` computed using the Maclaurin series expansion.
 */
template <typename T, std::size_t MACLAURIN_EXPANSION_REPEAT_NUMBER>
typename std::enable_if<std::is_same<T, float>::value, T>::type
exp_FloatingPoint_maclaurin_expansion_with_table(const T &x) {

  return Base::Math::exp_float_maclaurin_expansion_with_table<
      MACLAURIN_EXPANSION_REPEAT_NUMBER>(x);
}

/**
 * @brief Computes the exponential of a value using the Maclaurin series
 * expansion with a precomputed table.
 *
 * This function is a wrapper that calls the underlying implementation of the
 * Maclaurin series expansion for the exponential function, potentially
 * utilizing a lookup table for improved performance.
 *
 * @tparam T The floating-point type (e.g., float, double) for the computation.
 * @tparam MACLAURIN_EXPANSION_REPEAT_NUMBER The number of terms or repetitions
 * to use in the Maclaurin expansion.
 * @param x The input value for which to compute the exponential.
 * @return The computed exponential of x.
 */
template <typename T, std::size_t MACLAURIN_EXPANSION_REPEAT_NUMBER>
inline T exp_maclaurin_expansion_with_table(const T &x) {

  return exp_FloatingPoint_maclaurin_expansion_with_table<
      T, MACLAURIN_EXPANSION_REPEAT_NUMBER>(x);
}

/* exp */
/**
 * @brief Computes the exponential function e^x for the given input.
 *
 * This function provides multiple implementations of the exponential function,
 * depending on compile-time configuration macros:
 * - If `__BASE_MATH_USE_STD_MATH__` is defined, it uses the standard library's
 *   `std::exp` function, with input clamping to avoid overflow/underflow.
 * - If `__BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD__` is defined,
 *   it uses a Maclaurin series expansion with a lookup table for improved
 * accuracy.
 * - Otherwise, it uses a basic Maclaurin series expansion.
 *
 * Input values are clamped to predefined minimum and maximum values to ensure
 * numerical stability and prevent undefined behavior.
 *
 * @tparam T Numeric type of the input value.
 * @param x The exponent to raise e to.
 * @return The computed value of e^x, using the selected algorithm.
 */
template <typename T> inline T exp(const T &x) {

#ifdef __BASE_MATH_USE_STD_MATH__

  if (x > static_cast<float>(Base::Math::EXP_INPUT_MAX)) {
    return static_cast<float>(Base::Math::EXP_OUTPUT_MAX);
  } else if (x < static_cast<float>(Base::Math::EXP_INPUT_MIN)) {
    return static_cast<float>(Base::Math::EXP_OUTPUT_MIN);
  } else {

    return std::exp(x);
  }

#else // __BASE_MATH_USE_STD_MATH__

#ifdef __BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD__

  return exp_maclaurin_expansion_with_table<
      T, Base::Math::EXP_MACLAURIN_WITH_TABLE_REPEAT_NUMBER>(x);

#else // __BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD__

  return Base::Math::exp_maclaurin_expansion<T, Base::Math::EXP_REPEAT_NUMBER>(
      x);

#endif // __BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD__

#endif // __BASE_MATH_USE_STD_MATH__
}

/* exp2 */
namespace Exp2NewtonIteration {

template <typename T, std::size_t LOOP_MAX, std::size_t N> struct Loop {
  /**
   * @brief Computes and accumulates a term in a series expansion involving
   * logarithms.
   *
   * This static function multiplies the current term by (LN_2 * remainder) /
   * (LOOP_MAX - N), adds the result to the running total, and then recursively
   * calls itself with N - 1.
   *
   * @tparam T The numeric type used for computation.
   * @param result Reference to the accumulated result of the series.
   * @param term Reference to the current term in the series, which will be
   * updated.
   * @param remainder The remainder value used in the computation of the next
   * term.
   */
  static void compute(T &result, T &term, const T &remainder) {
    term *= static_cast<T>(Base::Math::LN_2) * remainder /
            static_cast<T>(LOOP_MAX - N);
    result += term;

    Loop<T, LOOP_MAX, N - 1>::compute(result, term, remainder);
  }
};

template <typename T, std::size_t LOOP_MAX> struct Loop<T, LOOP_MAX, 0> {
  /**
   * @brief Computes a value based on the provided result, term, and remainder.
   *
   * This is a placeholder implementation that performs no operation.
   * The function is intended to be specialized or overridden for specific
   * behavior.
   *
   * @tparam T The type of the input and output parameters.
   * @param result Reference to the result variable to be computed.
   * @param term Reference to the term variable used in computation.
   * @param remainder Constant reference to the remainder used in computation.
   */
  static void compute(T &result, T &term, const T &remainder) {
    /* Do Nothing. */
    static_cast<void>(result);
    static_cast<void>(term);
    static_cast<void>(remainder);
  }
};

} // namespace Exp2NewtonIteration

/**
 * @brief Computes the exponential function 2^x using a Maclaurin series
 * expansion.
 *
 * This function calculates 2 raised to the power of x (2^x) for a given input
 * x, using a Maclaurin series expansion for the fractional part of x and
 * integer exponentiation for the integer part. The computation is limited by
 * predefined input and output bounds.
 *
 * @tparam T The floating-point type for computation (e.g., float, double).
 * @tparam LOOP_NUMBER The number of terms to use in the Maclaurin series
 * expansion.
 * @param x The exponent value for which 2^x is to be computed.
 * @return The computed value of 2^x, clamped to output bounds if x is out of
 * range.
 */
template <typename T, std::size_t LOOP_NUMBER>
inline T exp2_maclaurin_expansion(const T &x) {

  T result = static_cast<T>(1);

  if (x > static_cast<T>(Base::Math::EXP2_INPUT_MAX)) {
    result = static_cast<T>(Base::Math::EXP2_OUTPUT_MAX);
  } else if (x < static_cast<T>(Base::Math::EXP2_INPUT_MIN)) {
    result = static_cast<T>(Base::Math::EXP2_OUTPUT_MIN);
  } else {

    int n = static_cast<int>(x);
    T remainder = x - static_cast<T>(n);

    T term = static_cast<T>(1);

    Exp2NewtonIteration::Loop<T, LOOP_NUMBER, LOOP_NUMBER - 1>::compute(
        result, term, remainder);

    result = fast_ldexp(result, n);
  }

  return result;
}

/**
 * @brief Computes the base-2 exponential of the input value.
 *
 * This function calculates 2 raised to the power of x (i.e., 2^x) for the given
 * input. Depending on the compilation flag `__BASE_MATH_USE_STD_MATH__`, it
 * either uses the standard library's `std::exp2` function or a custom Maclaurin
 * series expansion.
 *
 * If the input x exceeds predefined maximum or minimum bounds (`EXP2_INPUT_MAX`
 * or `EXP2_INPUT_MIN`), the function returns corresponding clamped output
 * values
 * (`EXP2_OUTPUT_MAX` or `EXP2_OUTPUT_MIN`) to prevent overflow or underflow.
 *
 * @tparam T Numeric type of the input and output (e.g., float, double).
 * @param x The exponent to which 2 is raised.
 * @return The computed value of 2^x, possibly clamped to avoid
 * overflow/underflow.
 */
template <typename T> inline T exp2(const T &x) {

#ifdef __BASE_MATH_USE_STD_MATH__

  if (x > static_cast<T>(Base::Math::EXP2_INPUT_MAX)) {
    return static_cast<T>(Base::Math::EXP2_OUTPUT_MAX);
  } else if (x < static_cast<T>(Base::Math::EXP2_INPUT_MIN)) {
    return static_cast<T>(Base::Math::EXP2_OUTPUT_MIN);
  } else {

    return std::exp2(x);
  }
#else // __BASE_MATH_USE_STD_MATH__
  return Base::Math::exp2_maclaurin_expansion<T,
                                              Base::Math::EXP2_REPEAT_NUMBER>(
      x);
#endif
}

/* log maclaurin expansion with table */

/**
 * @brief Approximates the natural logarithm (ln) of a double-precision
 * floating-point number using a Maclaurin series expansion and a precomputed
 * lookup table for improved performance.
 *
 * This function computes an approximation of ln(x) for a given positive double
 * x. It handles special cases for non-positive inputs by returning a predefined
 * minimum output. For positive inputs, it decomposes the floating-point number
 * into exponent and mantissa, applies a Maclaurin series expansion for the
 * mantissa, and uses a lookup table for part of the result.
 *
 * @param x The input value (must be greater than 0 for valid results).
 * @return The approximate value of ln(x). For x <= 0, returns
 * Base::Math::LOG_OUTPUT_MIN.
 *
 * @note
 * - The function assumes the existence of Base::Math::LOG_OUTPUT_MIN,
 * Base::Math::LN_2, and Base::Math::TABLE_FOR_LOG_DOUBLE[].
 * - The approximation is designed for performance and may not be as accurate as
 * std::log.
 * - The function uses bit manipulation and reinterpretation of the double's
 * binary representation.
 */
inline double log_double_maclaurin_expansion_with_table(const double &x) {

  double result = static_cast<double>(0);

  if (x <= static_cast<double>(0)) {
    result = static_cast<double>(Base::Math::LOG_OUTPUT_MIN);
  } else {

    unsigned long long w = static_cast<unsigned long long>(0);
    unsigned long long mantissa16 = static_cast<unsigned long long>(0);
    int q = static_cast<int>(0);
    double y = static_cast<double>(0);
    double h = static_cast<double>(0);
    double z = static_cast<double>(0);

    std::memcpy(&w, &x, static_cast<std::size_t>(8));

    q = ((static_cast<int>(w >> static_cast<int>(47)) &
          static_cast<int>(0x1F)) +
         static_cast<int>(1)) >>
        static_cast<int>(1);
    mantissa16 = (w & static_cast<unsigned long long>(0xFFFFFFFFFFFFF)) ^
                 static_cast<unsigned long long>(
                     0x4030000000000000); // mantissa*16  16<=mantissa16<32

    std::memcpy(&h, &mantissa16, static_cast<std::size_t>(8));

    z = static_cast<double>(q + static_cast<int>(16));
    h = (h - z) / (h + z);
    z = h * h;

    y = static_cast<double>(1.5) * z + static_cast<double>(2);
    y = y * h;

    result = static_cast<double>(static_cast<int>(w >> static_cast<int>(52)) -
                                 static_cast<int>(1023)) *
                 static_cast<double>(Base::Math::LN_2) +
             static_cast<double>(Base::Math::TABLE_FOR_LOG_DOUBLE[q]) + y;
  }

  return result;
}

/**
 * @brief Approximates the natural logarithm (ln) of a floating-point number
 * using a Maclaurin series expansion and a lookup table.
 *
 * This function computes an efficient approximation of the natural logarithm
 * for a given positive float value `x`. It uses bit manipulation to extract the
 * exponent and mantissa, applies a Maclaurin series expansion for improved
 * accuracy, and utilizes a precomputed lookup table for further optimization.
 *
 * @param x The input value for which to compute the natural logarithm. Must be
 * greater than 0.
 * @return The approximate value of ln(x). If x <= 0, returns a predefined
 * minimum output value (Base::Math::LOG_OUTPUT_MIN).
 *
 * @note
 * - The function assumes IEEE 754 floating-point representation for bit
 * manipulation.
 * - The accuracy depends on the quality of the Maclaurin expansion and the
 * lookup table (Base::Math::TABLE_FOR_LOG_DOUBLE).
 * - For x <= 0, the function returns a minimum output value to avoid domain
 * errors.
 *
 * @see Base::Math::LOG_OUTPUT_MIN
 * @see Base::Math::LN_2
 * @see Base::Math::TABLE_FOR_LOG_DOUBLE
 */
inline float log_float_maclaurin_expansion_with_table(const float &x) {

  float result = static_cast<float>(0);

  if (x <= static_cast<float>(0)) {
    result = static_cast<float>(Base::Math::LOG_OUTPUT_MIN);
  } else {

    unsigned long w = static_cast<unsigned long>(0);
    unsigned long mantissa16 = static_cast<unsigned long>(0);
    int q = static_cast<int>(0);
    float y = static_cast<float>(0);
    float h = static_cast<float>(0);
    float z = static_cast<float>(0);

    std::memcpy(&w, &x, static_cast<std::size_t>(4));

    q = ((static_cast<int>(w >> static_cast<int>(18)) &
          static_cast<int>(0x1F)) +
         static_cast<int>(1)) >>
        static_cast<int>(1);
    mantissa16 =
        (w & static_cast<unsigned long>(0x7FFFFF)) ^
        static_cast<unsigned long>(0x41800000); // mantissa*16 16<=mantissa16<32

    std::memcpy(&h, &mantissa16, static_cast<std::size_t>(4));

    z = static_cast<float>(q + static_cast<int>(16));
    h = (h - z) / (h + z);
    z = h * h;

    y = static_cast<float>(1.5) * z + static_cast<float>(2);
    y = y * h;

    result = static_cast<float>(static_cast<int>(w >> static_cast<int>(23)) -
                                static_cast<int>(127)) *
                 static_cast<float>(Base::Math::LN_2) +
             static_cast<float>(Base::Math::TABLE_FOR_LOG_DOUBLE[q]) + y;
  }

  return result;
}

/**
 * @brief Computes the natural logarithm of a floating-point value using a
 * Maclaurin series expansion with a lookup table for optimization.
 *
 * This function is enabled only for the type `double`. It delegates the
 * computation to `Base::Math::log_double_maclaurin_expansion_with_table`, which
 * performs the actual calculation.
 *
 * @tparam T The floating-point type (must be `double` for this specialization).
 * @param x The input value for which to compute the natural logarithm.
 * @return The natural logarithm of `x` as a value of type `T`.
 */
template <typename T>
typename std::enable_if<std::is_same<T, double>::value, T>::type
log_FloatingPoint_maclaurin_expansion_with_table(const T &x) {

  return Base::Math::log_double_maclaurin_expansion_with_table(x);
}

/**
 * @brief Computes the natural logarithm of a floating-point value using the
 * Maclaurin series expansion with a lookup table.
 *
 * This function is enabled only for the float type. It delegates the
 * computation to Base::Math::log_float_maclaurin_expansion_with_table, which
 * utilizes a Maclaurin series expansion and a precomputed table for improved
 * performance and accuracy.
 *
 * @tparam T The floating-point type (only enabled for float).
 * @param x The input value for which to compute the natural logarithm.
 * @return The computed natural logarithm of x as a float.
 */
template <typename T>
typename std::enable_if<std::is_same<T, float>::value, T>::type
log_FloatingPoint_maclaurin_expansion_with_table(const T &x) {

  return Base::Math::log_float_maclaurin_expansion_with_table(x);
}

/**
 * @brief Computes the natural logarithm of a value using the Maclaurin series
 * expansion with a lookup table for improved performance.
 *
 * This function is a wrapper that calls the underlying implementation in
 * Base::Math::log_FloatingPoint_maclaurin_expansion_with_table<T>. It is
 * intended for use with floating-point types and leverages a precomputed table
 * to accelerate the Maclaurin series expansion for logarithms.
 *
 * @tparam T The floating-point type of the input value.
 * @param x The value for which to compute the natural logarithm.
 * @return The natural logarithm of x, computed using the Maclaurin series with
 * table optimization.
 */
template <typename T> inline T log_maclaurin_expansion_with_table(const T &x) {

  return Base::Math::log_FloatingPoint_maclaurin_expansion_with_table<T>(x);
}

/* log newton method */

/**
 * @brief Performs a single iteration of Newton's method for logarithm
 * calculation.
 *
 * This struct template recursively applies Newton's method to refine the guess
 * for the logarithm of a value. On each iteration, it updates the exponential
 * guess and the current guess using the Newton-Raphson formula: guess = guess -
 * (exp(guess) - scaled_x) / exp(guess)
 *
 * @tparam T The floating-point type used for calculations.
 * @tparam N The number of remaining iterations to perform (recursively
 * decremented).
 * @param exp_guess Reference to the current exponential guess (exp(guess)).
 * @param guess Reference to the current guess for the logarithm.
 * @param scaled_x The scaled input value for which the logarithm is being
 * computed.
 */
template <typename T, std::size_t N> struct LogNewtonIterationLoop {
  static void compute(T &exp_guess, T &guess, const T &scaled_x) {
    exp_guess = Base::Math::exp(guess);
    guess = guess - (exp_guess - scaled_x) / exp_guess;

    LogNewtonIterationLoop<T, N - 1>::compute(exp_guess, guess, scaled_x);
  }
};

template <typename T> struct LogNewtonIterationLoop<T, 0> {
  /**
   * @brief Computes values related to exponential and logarithmic operations.
   *
   * This static function is intended to perform computations involving
   * exponential guesses and scaled input values. In its current implementation,
   * it does not perform any operations and simply suppresses unused parameter
   * warnings.
   *
   * @tparam T The data type of the input and output parameters.
   * @param exp_guess Reference to the variable representing the exponential
   * guess (output or intermediate value).
   * @param guess Reference to the variable representing the current guess
   * (output or intermediate value).
   * @param scaled_x Constant reference to the scaled input value.
   */
  static void compute(T &exp_guess, T &guess, const T &scaled_x) {
    /* Do Nothing. */
    static_cast<void>(exp_guess);
    static_cast<void>(guess);
    static_cast<void>(scaled_x);
  }
};

/**
 * @brief Computes the natural logarithm of a value using the Newton-Raphson
 * method.
 *
 * This function estimates the natural logarithm (ln) of the input value `x`
 * using an iterative Newton-Raphson approach. The number of iterations is
 * controlled by the template parameter `LOOP_NUMBER`. The input is scaled to
 * improve convergence, and the result is adjusted accordingly. Special handling
 * is provided for non-positive inputs, returning a minimum output value.
 *
 * @tparam T Numeric type of the input and output (e.g., float, double).
 * @tparam LOOP_NUMBER Number of Newton-Raphson iterations to perform.
 * @param x The value for which to compute the natural logarithm.
 * @return The estimated natural logarithm of `x`.
 *
 * @note The function relies on constants and helper classes/functions from the
 *       `Base::Math` namespace, such as `LOG_OUTPUT_MIN`, `LOG_SCALE_FACTOR`,
 *       `LOG_SCALING_LOOP_MAX`, and `LOG_OF_LOG_SCALE_FACTOR`.
 * @note For non-positive values of `x`, the function returns
 * `Base::Math::LOG_OUTPUT_MIN`.
 */
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

/* log */

/**
 * @brief Computes the natural logarithm (base e) of the input value.
 *
 * This function provides multiple implementations for calculating the natural
 * logarithm, depending on compile-time configuration macros:
 * - If `__BASE_MATH_USE_STD_MATH__` is defined, it uses the standard library's
 * `std::log` function, with a check to ensure the input is positive. If the
 * input is less than or equal to zero, it returns a predefined minimum output
 * value (`Base::Math::LOG_OUTPUT_MIN`).
 * - If `__BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD__` is defined
 * (and standard math is not used), it uses a Maclaurin expansion with a lookup
 * table.
 * - Otherwise, it uses Newton's method for logarithm calculation, with a
 * configurable number of iterations (`Base::Math::LOG_REPEAT_NUMBER`).
 *
 * @tparam T Numeric type of the input value.
 * @param x The value for which to compute the natural logarithm.
 * @return The natural logarithm of `x`, or a minimum output value if `x` is not
 * positive.
 */
template <typename T> inline T log(const T &x) {

#ifdef __BASE_MATH_USE_STD_MATH__

  if (x <= static_cast<T>(0)) {
    return static_cast<T>(Base::Math::LOG_OUTPUT_MIN);
  } else {

    return std::log(x);
  }
#else // __BASE_MATH_USE_STD_MATH__

#ifdef __BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD__

  return log_maclaurin_expansion_with_table<T>(x);

#else // __BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD__

  return Base::Math::log_newton_method<T, Base::Math::LOG_REPEAT_NUMBER>(x);

#endif // __BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD__

#endif // __BASE_MATH_USE_STD_MATH__
}

/* log2 newton method */
template <typename T, std::size_t LOOP_NUMBER>
inline T log2_newton_method(const T &x) {

  return Base::Math::log_newton_method<T, LOOP_NUMBER>(x) /
         static_cast<T>(Base::Math::LN_2);
}

/* log2 maclaurin expansion with table */

/**
 * @brief Computes the base-2 logarithm of a value using a Maclaurin series
 * expansion with a lookup table.
 *
 * This function calculates log2(x) by first computing the natural logarithm of
 * x using a Maclaurin series expansion (with possible table-based
 * optimizations), and then dividing the result by the natural logarithm of 2 to
 * convert it to base-2.
 *
 * @tparam T Numeric type of the input value (e.g., float, double).
 * @param x The value for which to compute the base-2 logarithm.
 * @return The base-2 logarithm of x.
 */
template <typename T> inline T log2_maclaurin_expansion_with_table(const T &x) {

  return Base::Math::log_maclaurin_expansion_with_table<T>(x) /
         static_cast<T>(Base::Math::LN_2);
}

/* log2 */

/**
 * @brief Computes the base-2 logarithm of the input value.
 *
 * This function calculates log2(x) using different algorithms depending on
 * compile-time macros:
 * - If `__BASE_MATH_USE_STD_MATH__` is defined, it uses the standard library's
 * `std::log2` function. If the input `x` is less than or equal to zero, it
 * returns a minimum output value defined by `Base::Math::LOG_OUTPUT_MIN`
 * divided by `Base::Math::LN_2`.
 * - If `__BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD__` is defined
 * (and standard math is not used), it uses a Maclaurin expansion with a lookup
 * table.
 * - Otherwise, it uses Newton's method with a repeat count defined by
 * `Base::Math::LOG_REPEAT_NUMBER`.
 *
 * @tparam T Numeric type of the input value.
 * @param x The value for which to compute the base-2 logarithm.
 * @return The base-2 logarithm of `x`.
 */
template <typename T> inline T log2(const T &x) {

#ifdef __BASE_MATH_USE_STD_MATH__

  if (x <= static_cast<T>(0)) {
    return static_cast<T>(Base::Math::LOG_OUTPUT_MIN) /
           static_cast<T>(Base::Math::LN_2);
  } else {

    return std::log2(x);
  }
#else // __BASE_MATH_USE_STD_MATH__

#ifdef __BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD__

  return Base::Math::log2_maclaurin_expansion_with_table<T>(x);

#else // __BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD__

  return Base::Math::log2_newton_method<T, Base::Math::LOG_REPEAT_NUMBER>(x);

#endif // __BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD__

#endif // __BASE_MATH_USE_STD_MATH__
}

/* log10 newton method  */

/**
 * @brief Computes the base-10 logarithm of a value using Newton's method.
 *
 * This function calculates the logarithm of the input value `x` to base 10
 * by first computing the natural logarithm using Newton's method and then
 * dividing by the natural logarithm of 10.
 *
 * @tparam T The floating-point type of the input and output.
 * @tparam LOOP_NUMBER The number of iterations for Newton's method.
 * @param x The value for which to compute the base-10 logarithm.
 * @return The base-10 logarithm of `x`.
 */
template <typename T, std::size_t LOOP_NUMBER>
inline T log10_newton_method(const T &x) {

  return Base::Math::log_newton_method<T, LOOP_NUMBER>(x) /
         static_cast<T>(Base::Math::LN_10);
}

/* log10 maclaurin expansion with table */

/**
 * @brief Computes the base-10 logarithm of a value using a Maclaurin series
 * expansion with a lookup table.
 *
 * This function calculates log10(x) by first computing the natural logarithm of
 * x using a Maclaurin series expansion (with possible table-based
 * optimizations), and then dividing the result by the natural logarithm of 10
 * to convert it to base-10.
 *
 * @tparam T Numeric type of the input and output (e.g., float, double).
 * @param x The value for which to compute the base-10 logarithm.
 * @return The base-10 logarithm of x.
 */
template <typename T>
inline T log10_maclaurin_expansion_with_table(const T &x) {

  return Base::Math::log_maclaurin_expansion_with_table<T>(x) /
         static_cast<T>(Base::Math::LN_10);
}

/* log10 */

/**
 * @brief Computes the base-10 logarithm of the input value.
 *
 * This function calculates the logarithm to base 10 for the given input `x`.
 * The implementation can use different algorithms depending on compile-time
 * macros:
 * - If `__BASE_MATH_USE_STD_MATH__` is defined, it uses the standard library's
 *   `std::log10` function. For non-positive values of `x`, it returns a minimum
 *   output value defined by `Base::Math::LOG_OUTPUT_MIN` divided by
 *   `Base::Math::LN_10`.
 * - If `__BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD__` is defined,
 *   it uses a Maclaurin expansion with a lookup table.
 * - Otherwise, it uses Newton's method with a specified number of iterations.
 *
 * @tparam T Numeric type of the input value.
 * @param x The value for which to compute the base-10 logarithm.
 * @return The base-10 logarithm of `x`.
 */
template <typename T> inline T log10(const T &x) {

#ifdef __BASE_MATH_USE_STD_MATH__

  if (x <= static_cast<T>(0)) {
    return static_cast<T>(Base::Math::LOG_OUTPUT_MIN) /
           static_cast<T>(Base::Math::LN_10);
  } else {

    return std::log10(x);
  }
#else // __BASE_MATH_USE_STD_MATH__

#ifdef __BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD__

  return Base::Math::log10_maclaurin_expansion_with_table<T>(x);

#else // __BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD__

  return Base::Math::log10_newton_method<T, Base::Math::LOG_REPEAT_NUMBER>(x);

#endif // __BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD__

#endif // __BASE_MATH_USE_STD_MATH__
}

/* pow maclaurin expansion */

/**
 * @brief Computes the power function x^y using Maclaurin expansion for
 * exponentiation and Newton's method for logarithm.
 *
 * This function calculates x raised to the power y by first computing the
 * natural logarithm of x using the Newton-Raphson method with a specified
 * number of iterations (LOG_LOOP_NUMBER), then multiplying the result by y, and
 * finally applying the Maclaurin series expansion for the exponential function
 * with a specified number of terms (EXP_LOOP_NUMBER).
 *
 * @tparam T                The numeric type (e.g., float, double).
 * @tparam EXP_LOOP_NUMBER  Number of terms to use in the Maclaurin expansion
 * for exponentiation.
 * @tparam LOG_LOOP_NUMBER  Number of iterations to use in the Newton-Raphson
 * method for logarithm.
 * @param x                 The base value.
 * @param y                 The exponent value.
 * @return T                The computed value of x raised to the power y.
 */
template <typename T, std::size_t EXP_LOOP_NUMBER, std::size_t LOG_LOOP_NUMBER>
inline T pow_base_math(const T &x, const T &y) {
  return Base::Math::exp_maclaurin_expansion<T, EXP_LOOP_NUMBER>(
      y * Base::Math::log_newton_method<T, LOG_LOOP_NUMBER>(x));
}

/**
 * @brief Computes x raised to the power y (x^y) using Maclaurin series
 * expansions for logarithm and exponential functions, with optional table-based
 * optimization.
 *
 * This function calculates x^y by first computing the natural logarithm of x
 * using a Maclaurin series expansion (potentially with a lookup table for
 * optimization), then multiplying the result by y, and finally applying the
 * exponential function (also via Maclaurin expansion with table) to obtain the
 * final result.
 *
 * @tparam T The numeric type (e.g., float, double) for the computation.
 * @tparam EXP_LOOP_NUMBER The number of terms to use in the Maclaurin series
 * expansion for the exponential function.
 * @param x The base value.
 * @param y The exponent value.
 * @return T The computed value of x raised to the power y.
 *
 * @note This method is suitable for cases where high performance and reasonable
 * accuracy are required, and is especially useful in environments where
 * standard library functions may not be available or are too slow.
 */
template <typename T, std::size_t EXP_LOOP_NUMBER>
inline T pow_maclaurin_expansion_with_table(const T &x, const T &y) {
  return Base::Math::exp_maclaurin_expansion_with_table<T, EXP_LOOP_NUMBER>(
      y * Base::Math::log_maclaurin_expansion_with_table<T>(x));
}

/* pow */

/**
 * @brief Computes the power of a number (x raised to the power y).
 *
 * This function calculates x^y for the given values of x and y.
 * Depending on the compilation flag __BASE_MATH_USE_STD_MATH__, it either:
 * - Uses the standard library's std::pow function, or
 * - Uses a custom implementation based on the Maclaurin expansion with a lookup
 * table.
 *
 * @tparam T The numeric type of the input values.
 * @param x The base value.
 * @param y The exponent value.
 * @return The result of x raised to the power y.
 */
template <typename T> inline T pow(const T &x, const T &y) {

#ifdef __BASE_MATH_USE_STD_MATH__
  return std::pow(x, y);
#else // __BASE_MATH_USE_STD_MATH__

  return Base::Math::pow_maclaurin_expansion_with_table<
      T, Base::Math::EXP_MACLAURIN_WITH_TABLE_REPEAT_NUMBER>(x, y);

#endif // __BASE_MATH_USE_STD_MATH__
}

} // namespace Math
} // namespace Base

#endif // __BASE_MATH_EXPONENTIAL_LOGARITHMIC_HPP__
