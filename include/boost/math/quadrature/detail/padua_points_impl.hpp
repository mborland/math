// Copyright Matt Borland, 2021
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// References:
// Summary: https://www.math.unipd.it/~marcov/pdf/poster5ecm.pdf
// https://www.math.unipd.it/~marcov/pdf/fastpadua.pdf
// https://www.math.unipd.it/~demarchi/papers/padua2d.pdf

#ifndef BOOST_MATH_QUADRATURE_DETAIL_PADUA_POINTS_IMPL
#define BOOST_MATH_QUADRATURE_DETAIL_PADUA_POINTS_IMPL

#include <cstddef>

namespace boost { namespace math { namespace quadrature { namespace detail {

// Constrained by RealType concept at user interface
template<typename Real>
class padua_points_impl
{
private:
    std::size_t levels_;

public:
    explicit padua_points_impl(std::size_t levels) : levels_ {levels} {};
};

}}}} // namespaces

#endif // BOOST_MATH_QUADRATURE_DETAIL_PADUA_POINTS_IMPL
