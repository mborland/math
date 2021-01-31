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

#ifndef BOOST_MATH_QUADRATURE_PADUA_POINTS
#define BOOST_MATH_QUADRATURE_PADUA_POINTS

#include <vector>
#include <memory>
#include <type_traits>
#include <cstddef>
#include <boost/math/quadrature/detail/padua_points_impl.hpp>
#include <boost/math/concepts/concepts.hpp>

namespace boost { namespace math { namespace quadrature {

template<RealType Real>
class padua_points
{
private:
    std::unique_ptr<detail::padua_points_impl<Real>> impl_;

public:
    explicit padua_points(std::size_t levels) : impl_(std::make_unique<detail::padua_points_impl<Real>>(levels)) {};
};

}}} // Namespaces
#endif // BOOST_MATH_QUADRATURE_PADUA_POINTS
