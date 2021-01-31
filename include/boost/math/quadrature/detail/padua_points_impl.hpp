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

#include <vector>
#include <cstddef>
#include <cmath>
#include <boost/math/constants/constants.hpp>
#include <boost/math/concepts/concepts.hpp>

namespace boost { namespace math { namespace quadrature { namespace detail {

template<RealType Real>
class padua_points_impl
{
private:
    std::size_t levels_;
    std::size_t order_;
    std::vector<Real> points_;

    void calculate_order();
    void calculate_points();

public:
    explicit padua_points_impl(std::size_t levels) : levels_ {levels} 
    {
        calculate_order();
        calculate_points();
    };
};

template<RealType Real>
void padua_points_impl<Real>::calculate_order()
{
    for(std::size_t i = 0; i <= levels_; ++i)
    {
        order_ += levels_/2 + 1;
        if(levels_%2 == 1 && i%2 == 1)
        {
            ++order_;
        }
    }
}

template<RealType Real>
void padua_points_impl<Real>::calculate_points()
{
    const Real pi = boost::math::constants::pi<Real>();
    const std::size_t num_points = ((levels_+1)*(levels_+2))/2;
    points_.reserve(2*num_points);

    if(levels_ == 0)
    {
        points_[0] = Real(0);
        points_[1] = Real(0);
        return;
    }

    for(std::size_t i = 0, k = 0; i <= levels_; ++i, ++k)
    {
        std::size_t j_hi = levels_/2 + 1;
        if(levels_%2 == 1 && i%2 == 1)
        {
            ++j_hi;
        }

        for(std::size_t j = 0; j <= j_hi; ++j)
        {
            if(i*2 == levels_)
            {
                points_[0+k*2] = Real(0);
            }
            else
            {
                Real angle_1 = i*pi/levels_;
                points_[0+k*2] = std::cos(angle_1);
            }

            if(i%2 == 0)
            {
                if(2*(2*j-1) == levels_+1)
                {
                    points_[1+k*2] = Real(0);
                }
                else
                {
                    Real angle_2 = (2*j-1)*pi/(levels_+1);
                    points_[1+k*2] = std::cos(angle_2);
                }
            }
            else
            {
                if(2*(2*j-2) == levels_+1)
                {
                    points_[1+k*2] = Real(0);
                }
                else
                {
                    Real angle_2 = (2*j-2)*pi/(levels_+1);
                    points_[1+k*2] = std::cos(angle_2); 
                }
            }
        }
    }
}

}}}} // namespaces

#endif // BOOST_MATH_QUADRATURE_DETAIL_PADUA_POINTS_IMPL
