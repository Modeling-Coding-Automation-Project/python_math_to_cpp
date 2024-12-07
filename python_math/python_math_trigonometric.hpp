#ifndef PYTHON_MATH_TRIGONOMETRIC_HPP
#define PYTHON_MATH_TRIGONOMETRIC_HPP

#include "base_math.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace PythonMath {

/* sin */
template <typename T> inline T sin(const T &x) { return Base::Math::sin(x); }

template <typename T> inline std::vector<T> sin(const std::vector<T> &vector) {
  std::vector<T> result;
  result.reserve(vector.size());
  for (const auto &element : vector) {
    result.push_back(Base::Math::sin(element));
  }
  return result;
}

template <typename T, std::size_t N>
inline std::array<T, N> sin(const std::array<T, N> &array) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::sin(array[i]);
  }
  return result;
}

/* cos */
template <typename T> inline T cos(const T &x) { return Base::Math::cos(x); }

template <typename T> inline std::vector<T> cos(const std::vector<T> &vector) {
  std::vector<T> result;
  result.reserve(vector.size());
  for (const auto &element : vector) {
    result.push_back(Base::Math::cos(element));
  }
  return result;
}

template <typename T, std::size_t N>
inline std::array<T, N> cos(const std::array<T, N> &array) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::cos(array[i]);
  }
  return result;
}

/* tan */
template <typename T> inline T tan(const T &x) { return Base::Math::tan(x); }

template <typename T> inline std::vector<T> tan(const std::vector<T> &vector) {
  std::vector<T> result;
  result.reserve(vector.size());
  for (const auto &element : vector) {
    result.push_back(Base::Math::tan(element));
  }
  return result;
}

template <typename T, std::size_t N>
inline std::array<T, N> tan(const std::array<T, N> &array) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::tan(array[i]);
  }
  return result;
}

/* atan */
template <typename T> inline T atan(const T &x) { return Base::Math::atan(x); }

template <typename T> inline std::vector<T> atan(const std::vector<T> &vector) {
  std::vector<T> result;
  result.reserve(vector.size());
  for (const auto &element : vector) {
    result.push_back(Base::Math::atan(element));
  }
  return result;
}

template <typename T, std::size_t N>
inline std::array<T, N> atan(const std::array<T, N> &array) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::atan(array[i]);
  }
  return result;
}

/* atan2 */
template <typename T> inline T atan2(const T &y, const T &x) {
  return Base::Math::atan2(y, x);
}

template <typename T>
inline std::vector<T> atan2(const std::vector<T> &vector_y, const T &x) {
  std::vector<T> result;
  result.reserve(vector_y.size());
  for (const auto &element : vector_y) {
    result.push_back(Base::Math::atan2(element, x));
  }
  return result;
}

template <typename T>
inline std::vector<T> atan2(const T &y, const std::vector<T> &vector_x) {
  std::vector<T> result;
  result.reserve(vector_x.size());
  for (const auto &element : vector_x) {
    result.push_back(Base::Math::atan2(y, element));
  }
  return result;
}

template <typename T>
inline std::vector<T> atan2(const std::vector<T> &vector_y,
                            const std::vector<T> &vector_x) {
  std::vector<T> result;
  result.reserve(vector_y.size());
  for (std::size_t i = 0; i < vector_y.size(); ++i) {
    result.push_back(Base::Math::atan2(vector_y[i], vector_x[i]));
  }
  return result;
}

template <typename T, std::size_t N>
inline std::array<T, N> atan2(const std::array<T, N> &array_y, const T &x) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::atan2(array_y[i], x);
  }
  return result;
}

template <typename T, std::size_t N>
inline std::array<T, N> atan2(const T &y, const std::array<T, N> &array_x) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::atan2(y, array_x[i]);
  }
  return result;
}

template <typename T, std::size_t N>
inline std::array<T, N> atan2(const std::array<T, N> &array_y,
                              const std::array<T, N> &array_x) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::atan2(array_y[i], array_x[i]);
  }
  return result;
}

/* asin */
template <typename T> inline T asin(const T &x) { return Base::Math::asin(x); }

template <typename T> inline std::vector<T> asin(const std::vector<T> &vector) {
  std::vector<T> result;
  result.reserve(vector.size());
  for (const auto &element : vector) {
    result.push_back(Base::Math::asin(element));
  }
  return result;
}

template <typename T, std::size_t N>
inline std::array<T, N> asin(const std::array<T, N> &array) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::asin(array[i]);
  }
  return result;
}

/* acos */
template <typename T> inline T acos(const T &x) { return Base::Math::acos(x); }

template <typename T> inline std::vector<T> acos(const std::vector<T> &vector) {
  std::vector<T> result;
  result.reserve(vector.size());
  for (const auto &element : vector) {
    result.push_back(Base::Math::acos(element));
  }
  return result;
}

template <typename T, std::size_t N>
inline std::array<T, N> acos(const std::array<T, N> &array) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::acos(array[i]);
  }
  return result;
}

/* sinh */
template <typename T> inline T sinh(const T &x) { return Base::Math::sinh(x); }

template <typename T> inline std::vector<T> sinh(const std::vector<T> &vector) {
  std::vector<T> result;
  result.reserve(vector.size());
  for (const auto &element : vector) {
    result.push_back(Base::Math::sinh(element));
  }
  return result;
}

template <typename T, std::size_t N>
inline std::array<T, N> sinh(const std::array<T, N> &array) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::sinh(array[i]);
  }
  return result;
}

/* cosh */
template <typename T> inline T cosh(const T &x) { return Base::Math::cosh(x); }

template <typename T> inline std::vector<T> cosh(const std::vector<T> &vector) {
  std::vector<T> result;
  result.reserve(vector.size());
  for (const auto &element : vector) {
    result.push_back(Base::Math::cosh(element));
  }
  return result;
}

template <typename T, std::size_t N>
inline std::array<T, N> cosh(const std::array<T, N> &array) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::cosh(array[i]);
  }
  return result;
}

/* tanh */
template <typename T> inline T tanh(const T &x) { return Base::Math::tanh(x); }

template <typename T> inline std::vector<T> tanh(const std::vector<T> &vector) {
  std::vector<T> result;
  result.reserve(vector.size());
  for (const auto &element : vector) {
    result.push_back(Base::Math::tanh(element));
  }
  return result;
}

template <typename T, std::size_t N>
inline std::array<T, N> tanh(const std::array<T, N> &array) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::tanh(array[i]);
  }
  return result;
}

} // namespace PythonMath

#endif // PYTHON_MATH_TRIGONOMETRIC_HPP
