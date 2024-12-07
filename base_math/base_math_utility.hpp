#ifndef BASE_MATH_UTILITY_HPP
#define BASE_MATH_UTILITY_HPP

namespace Base {
namespace Math {

template <typename T> inline T avoid_zero_divide(T in, const T &division_min) {
  if (in < division_min) {
    if (in >= 0) {
      return division_min;
    } else if (in > -division_min) {
      return -division_min;
    }
  }

  return in;
}

template <typename T> inline bool near_zero(T in, T division_min) {
  bool flag = false;
  if (in < division_min) {
    if (in >= 0) {
      flag = true;
    } else if (in > -division_min) {
      flag = true;
    }
  }

  return flag;
}

} // namespace Math
} // namespace Base

#endif // BASE_MATH_UTILITY_HPP
