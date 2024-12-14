#ifndef BASE_MATH_UTILITY_HPP
#define BASE_MATH_UTILITY_HPP

namespace Base {
namespace Math {

template <typename T> inline T ldexp(const T &x, const int &n) {

  T result = x;

  for (int i = 0; i < n; ++i) {
    result *= static_cast<T>(2);
  }
  for (int i = 0; i < -n; ++i) {
    result *= static_cast<T>(0.5);
  }

  return result;
}

} // namespace Math
} // namespace Base

#endif // BASE_MATH_UTILITY_HPP
