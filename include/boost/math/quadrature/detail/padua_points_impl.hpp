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
#include <iostream>
#include <cstddef>
#include <cmath>
#include <boost/math/constants/constants.hpp>
#include <boost/math/concepts/concepts.hpp>

namespace boost { namespace math { namespace quadrature { namespace detail {

template<SimpleRealType Real>
class padua_points_impl
{
private:
    std::size_t levels_;
    std::size_t order_;
    std::vector<Real> points_;

    void calculate_order();
    void calculate_points();
    void print_points();

public:
    explicit padua_points_impl(std::size_t levels) : levels_ {levels} 
    {
        std::cout << "Levels: " << levels_ << '\n';
        calculate_order();
        std::cout << "Order: " << order_ << '\n';
        calculate_points();
        print_points();
    };
};

template<SimpleRealType Real>
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

template<SimpleRealType Real>
void padua_points_impl<Real>::calculate_points()
{
    const Real pi = boost::math::constants::pi<Real>();
    const std::size_t num_points = ((levels_+1)*(levels_+2));
    points_.resize(num_points);

    if(levels_ == 0)
    {
        points_[0] = Real(0);
        points_[1] = Real(0);
        return;
    }

    std::size_t index = 0;
    for(std::size_t i = 0; i <= levels_; ++i)
    {
        std::size_t j_max = (levels_/2) + 1;
        if(levels_ % 2 == 1 && i % 2 == 1)
        {
            ++j_max;
        }

        for(std::size_t j = 1; j <= j_max; ++j)
        {
            if(i * 2 == levels_)
            {
                points_[index*2] = Real(0);
            }
            else
            {
                points_[index*2] = std::cos(i*pi/levels_);
            }

            if(i % 2 == 0)
            {
                if(2*(2*j-1) == levels_+1)
                {
                    points_[1+index*2] = Real(0);
                }
                else
                {
                    points_[1+index*2] = std::cos((2*j-1)*pi/(levels_+1));
                }
            }

            else
            {
                if(2*(2*j-2) == levels_+1)
                {
                    points_[1+index*2] = Real(0);
                }
                else
                {
                    points_[1+index*2] = std::cos((2*j-2)*pi/(levels_+1));
                }
            }
            ++index;
        }
    }
}

template<SimpleRealType Real>
void padua_points_impl<Real>::print_points()
{
    for(std::size_t i = 0; i < points_.size(); ++i)
    {
        std::cout << "X: " << i << " Y: " << points_[i] << std::endl;
    }
    std::cout << std::endl;
}
}}}} // namespaces

#endif // BOOST_MATH_QUADRATURE_DETAIL_PADUA_POINTS_IMPL
