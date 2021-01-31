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

#ifdef __has_include && defined __has_include(<concepts>)
#include <concepts>

template<typename T>
concept RealType = std::is_floating_point_v<T>;

#else
#define RealType typename
#endif

#endif // BOOST_MATH_CONCEPTS
