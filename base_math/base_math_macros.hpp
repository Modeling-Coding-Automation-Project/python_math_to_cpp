#ifndef __BASE_MATH_MACROS_HPP__
#define __BASE_MATH_MACROS_HPP__

namespace Base {
namespace Math {

/* Comment out if you want to use the standard math library. */
// #define __BASE_MATH_USE_STD_MATH__

#ifndef __BASE_MATH_USE_STD_MATH__

/* Comment out if you want to use the algorithms fast but dependent on IEEE 754
 * standard. */
// #define __BASE_MATH_USE_ALGORITHM_DEPENDENT_ON_IEEE_754_STANDARD__

/* Comment out if you want to use the fast but rough approximations. */
// #define __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__

#else  // __BASE_MATH_USE_STD_MATH__
#endif // __BASE_MATH_USE_STD_MATH__

} // namespace Math
} // namespace Base

#endif // __BASE_MATH_MACROS_HPP__
