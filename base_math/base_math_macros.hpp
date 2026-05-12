#ifndef BASE_MATH_MACROS_HPP_
#define BASE_MATH_MACROS_HPP_

namespace Base {
namespace Math {

/* Comment out if you want to use the standard math library. */
// #define BASE_MATH_USE_STD_MATH_

#ifndef BASE_MATH_USE_STD_MATH_

/* Comment out if you want to use the algorithms fast but dependent on IEEE 754
 * standard. */
// #define BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD_

/* Comment out if you want to use the fast but rough approximations. */
// #define BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS_

#else  // BASE_MATH_USE_STD_MATH_
#endif // BASE_MATH_USE_STD_MATH_

} // namespace Math
} // namespace Base

#endif // BASE_MATH_MACROS_HPP_
