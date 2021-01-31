// Copyright Matt Borland, 2021
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Standardized C++20 concepts (with backwards compatiable macros)
// https://www.stroustrup.com/good_concepts.pdf

#ifndef BOOST_MATH_CONCEPTS
#define BOOST_MATH_CONCEPTS

#if __cplusplus > 202000L || _MSVC_LANG > 202000L && defined __has_include && __has_include(<concepts>)
#include <concepts>

// RealType that excludes multiprecision types
template<typename T>
concept SimpleRealType = std::is_floating_point_v<T>;

template<typename T>
concept Functor = std::is_function_v<T>;

#else
#define SimpleRealType typename
#define Functor typename
#endif

#endif // BOOST_MATH_CONCEPTS
