#ifndef BASE_MATH_MACROS_HPP
#define BASE_MATH_MACROS_HPP

namespace Base {
namespace Math {

/* Comment out if you want to use the standard math library. */
// #define BASE_MATH_USE_STD_MATH

#ifndef BASE_MATH_USE_STD_MATH

/* Comment out if you want to use the algorithms fast but dependent on IEEE 754
 * standard. */
// #define BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD

/* Comment out if you want to use the fast but rough approximations. */
// #define BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS

#else  // BASE_MATH_USE_STD_MATH
#endif // BASE_MATH_USE_STD_MATH

} // namespace Math
} // namespace Base

#endif // BASE_MATH_MACROS_HPP
