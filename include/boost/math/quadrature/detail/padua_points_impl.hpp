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

namespace boost { namespace math { namespace quadrature { namespace detail {

template<typename Real>
class padua_points_impl
{
private:
    std::size_t levels_;
    std::size_t num_points_;
    std::size_t order_;
    std::vector<Real> points_;
    std::vector<Real> weights_;

    void calculate_order();
    void calculate_points();
    void print_points();
    void calculate_weights();

public:
    explicit padua_points_impl(std::size_t levels) : levels_ {levels}, num_points_ {(levels+1)*(levels+2)} 
    {
        std::cout << "Levels: " << levels_ << '\n';
        
        calculate_order();
        std::cout << "Order: " << order_ << '\n';
        
        calculate_points();
        print_points();
        
        calculate_weights();
    };
};

template<typename Real>
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

template<typename Real>
void padua_points_impl<Real>::calculate_points()
{
    using std::cos;
    constexpr Real pi = boost::math::constants::pi<Real>();
    
    points_.resize(num_points_);

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
                points_[index*2] = cos(i*pi/levels_);
            }

            if(i % 2 == 0)
            {
                if(2*(2*j-1) == levels_+1)
                {
                    points_[1+index*2] = Real(0);
                }
                else
                {
                    points_[1+index*2] = cos((2*j-1)*pi/(levels_+1));
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
                    points_[1+index*2] = cos((2*j-2)*pi/(levels_+1));
                }
            }
            ++index;
        }
    }
}

template<typename Real>
void padua_points_impl<Real>::print_points()
{
    for(std::size_t i = 0; i < points_.size(); ++i)
    {
        std::cout << "X: " << i << " Y: " << points_[i] << std::endl;
    }
    std::cout << std::endl;
}

template<typename Real>
void padua_points_impl<Real>::calculate_weights()
{
    using std::cos;
    constexpr Real pi = boost::math::constants::pi<Real>();
    constexpr Real sqrt2 = boost::math::constants::root_two<Real>();

    weights_.resize(num_points_/2);

    if(levels_ == 0)
    {
        weights_[0] = Real(4);
        return;
    }

    // Levels to represent subgrids
    const std::size_t levels_1 = (levels_+1)/2;
    const std::size_t levels_2 = (levels_+2)/2;
    const std::size_t levels_3 = (levels_+3)/2;

    // Calculate even and odd Chebyshev polynomials on subgrids 1 and 2
    std::vector<Real> chebyshev_even_1(levels_2 * levels_2);
    std::vector<Real> chebyshev_odd_1(levels_2 * levels_1);
    std::vector<Real> chebyshev_even_2(levels_2 * levels_3);
    std::vector<Real> chebyshev_odd_2(levels_2 * levels_2);

    for(std::size_t i = 0; i < levels_2; ++i)
    {
        chebyshev_even_1[i*levels_2] = Real(1);
        for(std::size_t j = 1; j < levels_2; ++j)
        {
            chebyshev_even_1[i*levels_2+j] = cos(4*i*j/levels_) * sqrt2;
        }
    }

    for(std::size_t i = 0; i < levels_1; ++i)
    {
        chebyshev_odd_1[i*levels_2] = Real(1);
        for(std::size_t j = 1; j < levels_2; ++j)
        {
            chebyshev_odd_1[i*levels_2+j] = cos(2*j*(2*i+1)/levels_) * sqrt2;
        }
    }

    for(std::size_t i = 0; i < levels_3; ++i)
    {
        chebyshev_even_2[i*levels_2] = Real(1);
        for(std::size_t j = 1; j < levels_2; ++j)
        {
            chebyshev_even_2[i*levels_2+j] = cos(4*i*j/(levels_+1)) * sqrt2;
        }
    }

    for(std::size_t i = 0; i < levels_2; ++i)
    {
        chebyshev_odd_2[i*levels_2] = Real(1);
        for(std::size_t j = 1; j < levels_2; ++j)
        {
            chebyshev_odd_2[i*levels_2+j] = cos(2*j*(2*i+1)/(levels_+1)) * sqrt2;
        }
    }

    // Calculate the moments matrix
    std::vector<Real> moments(levels_2 * levels_2);

    for(std::size_t i = 0; i < levels_2; ++i)
    {
        Real moment_i = 2*sqrt2/(1-(4*i*i));
        for(std::size_t j = 0; j < levels_2 - i; ++j)
        {
            Real moment_j = 2*sqrt2/(1-(4*j*j));
            moments[i*levels_2+j] = moment_i * moment_j;
        }
    }

    for(std::size_t i = 0; i < levels_2; ++i)
    {
        moments[i*levels_2] /= sqrt2;
        moments[i] /= sqrt2;
    }

    if((levels_ % 2) == 0)
    {
        moments[levels_2-1] /= 2;
    }

    // Calculate matrix products
    std::vector<Real> odd_even_product(levels_2*levels_2, Real(0));

    for(std::size_t i = 0; i < levels_2; ++i)
    {
        for(std::size_t j = 0; j < levels_2 - i; ++j)
        {
            for(std::size_t i2 = 0; i2 < levels_2; ++i2)
            {
                for(std::size_t j2 = 0; j2 < levels_2; ++j2)
                {
                    odd_even_product[j2+i2*levels_2] += chebyshev_odd_2[j2*levels_2+j] * moments[j*levels_2+i] * chebyshev_even_1[i2*levels_2+i];
                }
            }
        }
    }

    std::vector<Real> even_odd_product(levels_3*levels_1, Real(0));

    for(std::size_t i = 0; i < levels_2; ++i)
    {
        for(std::size_t j = 0; j < levels_2 - i; ++j)
        {
            for(std::size_t i2 = 0; i2 < levels_1; ++i2)
            {
                for(std::size_t j2 = 0; j < levels_3; ++j2)
                {
                    even_odd_product[i2*levels_3+j2] += chebyshev_even_2[j2*levels_2+j] * moments[j*levels_2+i] * chebyshev_odd_1[i2*levels_2+i];
                }
            }
        }
    }

    // Calculate interpolation weight matricies
    std::vector<Real> weights_1(levels_2*levels_2);

    for(std::size_t i = 0; i < levels_2; ++i)
    {
        for(std::size_t j = 0; j < levels_2; ++j)
        {
            weights_1[i*levels_2 + j] = Real(2)/(levels_*(levels_+1));
        }
        weights_1[i] /= Real(2);
    }

    if(levels_ % 2 == 0)
    {
        for(std::size_t i = 0; i < levels_2; ++i)
        {
            weights_1[i+levels_2*levels_2-1] /= Real(2);
        }

        for(std::size_t i = 0; i < levels_2; ++i)
        {
            weights_1[i*levels_2+levels_2-1] /= Real(2);
        }
    }

    std::vector<Real> weights_2(levels_1*levels_3);
}

}}}} // namespaces

#endif // BOOST_MATH_QUADRATURE_DETAIL_PADUA_POINTS_IMPL
