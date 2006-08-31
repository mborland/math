// Copyright John Maddock 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_DISTRIBUTIONS_FISHER_F_HPP
#define BOOST_MATH_DISTRIBUTIONS_FISHER_F_HPP

#include <boost/math/special_functions/beta.hpp> // for incomplete beta.
#include <boost/math/distributions/complement.hpp> // complements
#include <boost/math/distributions/detail/common_error_handling.hpp> // error checks
#include <boost/math/special_functions/fpclassify.hpp>

namespace boost{ namespace math{ 
   
template <class RealType>
class fisher_f_distribution
{
public:
   typedef RealType value_type;

   fisher_f_distribution(const RealType& i, const RealType& j) : m_df1(i), m_df2(j)
   {
      RealType result;
      detail::check_df(
         BOOST_CURRENT_FUNCTION, m_df1, &result);
      detail::check_df(
         BOOST_CURRENT_FUNCTION, m_df2, &result);
   } // fisher_f_distribution

   RealType degrees_of_freedom1()const
   {
      return m_df1;
   }
   RealType degrees_of_freedom2()const
   {
      return m_df2;
   }

private:
   //
   // Data members:
   //
   RealType m_df1;  // degrees of freedom are a real number.
   RealType m_df2;  // degrees of freedom are a real number.
};

typedef fisher_f_distribution<double> fisher_f;

template <class RealType>
RealType pdf(const fisher_f_distribution<RealType>& dist, const RealType& x)
{
   using namespace std;  // for ADL of std functions
   RealType df1 = dist.degrees_of_freedom1();
   RealType df2 = dist.degrees_of_freedom2();
   // Error check:
   RealType error_result;
   if(false == detail::check_df(
         BOOST_CURRENT_FUNCTION, df1, &error_result)
         && detail::check_df(
         BOOST_CURRENT_FUNCTION, df2, &error_result))
      return error_result;

   if((x < 0) || !(boost::math::isfinite)(x))
   {
      return tools::domain_error<RealType>(
         BOOST_CURRENT_FUNCTION, "Random variable parameter was %1%, but must be > 0 !", x);
   }

   if(x == 0)
   {
      // special cases:
      if(df1 < 2)
         return tools::overflow_error<RealType>(
            BOOST_CURRENT_FUNCTION, 0);
      else if(df1 == 2)
         return 1;
      else
         return 0;
   }

   //
   // You reach this formula by direct differentiation of the
   // cdf expressed in terms of the incomplete beta.
   //
   // There are two versions so we don't pass a value of z
   // that is very close to 1 to ibeta_derivative: for some values
   // of df1 and df2, all the change takes place in this area.
   //
   RealType v1x = df1 * x;
   RealType result;
   if(v1x > df2)
   {
      result = (df2 * df1) / ((df2 + v1x) * (df2 + v1x));
      result *= ibeta_derivative(df2 / 2, df1 / 2, df2 / (df2 + v1x));
   }
   else
   {
      result = df2 + df1 * x;
      result = (result * df1 - x * df1 * df1) / (result * result);
      result *= ibeta_derivative(df1 / 2, df2 / 2, v1x / (df2 + v1x));
   }
   return result;
} // pdf

template <class RealType>
RealType cdf(const fisher_f_distribution<RealType>& dist, const RealType& x)
{
   RealType df1 = dist.degrees_of_freedom1();
   RealType df2 = dist.degrees_of_freedom2();
   // Error check:
   RealType error_result;
   if(false == detail::check_df(
         BOOST_CURRENT_FUNCTION, df1, &error_result)
         && detail::check_df(
         BOOST_CURRENT_FUNCTION, df2, &error_result))
      return error_result;

   if((x < 0) || !(boost::math::isfinite)(x))
   {
      return tools::domain_error<RealType>(
         BOOST_CURRENT_FUNCTION, "Random Variable parameter was %1%, but must be > 0 !", x);
   }

   RealType v1x = df1 * x;
   //
   // There are two equivalent formulas used here, the aim is
   // to prevent the final argument to the incomplete beta
   // from being too close to 1: for some values of df1 and df2
   // the rate of change can be arbitrarily large in this area,
   // whilst the value we're passing will have lost information
   // content as a result of being 0.999999something.  Better
   // to switch things around so we're passing 1-z instead.
   //
   return v1x > df2 
      ? boost::math::ibetac(df2 / 2, df1 / 2, df2 / (df2 + v1x))
      : boost::math::ibeta(df1 / 2, df2 / 2, v1x / (df2 + v1x));
} // cdf

template <class RealType>
RealType quantile(const fisher_f_distribution<RealType>& dist, const RealType& p)
{
   RealType df1 = dist.degrees_of_freedom1();
   RealType df2 = dist.degrees_of_freedom2();
   // Error check:
   RealType error_result;
   if(false == detail::check_df(
            BOOST_CURRENT_FUNCTION, df1, &error_result)
         && detail::check_df(
            BOOST_CURRENT_FUNCTION, df2, &error_result)
         && detail::check_probability(
            BOOST_CURRENT_FUNCTION, p, &error_result))
      return error_result;

   RealType x, y;

   x = boost::math::ibeta_inv(df1 / 2, df2 / 2, p, &y);

   return df2 * x / (df1 * y);
} // quantile

template <class RealType>
RealType cdf(const complemented2_type<fisher_f_distribution<RealType>, RealType>& c)
{
   RealType df1 = c.dist.degrees_of_freedom1();
   RealType df2 = c.dist.degrees_of_freedom2();
   RealType x = c.param;
   // Error check:
   RealType error_result;
   if(false == detail::check_df(
         BOOST_CURRENT_FUNCTION, df1, &error_result)
         && detail::check_df(
         BOOST_CURRENT_FUNCTION, df2, &error_result))
      return error_result;

   if((x < 0) || !(boost::math::isfinite)(x))
   {
      return tools::domain_error<RealType>(
         BOOST_CURRENT_FUNCTION, "Random Variable parameter was %1%, but must be > 0 !", x);
   }

   RealType v1x = df1 * x;
   //
   // There are two equivalent formulas used here, the aim is
   // to prevent the final argument to the incomplete beta
   // from being too close to 1: for some values of df1 and df2
   // the rate of change can be arbitrarily large in this area,
   // whilst the value we're passing will have lost information
   // content as a result of being 0.999999something.  Better
   // to switch things around so we're passing 1-z instead.
   //
   return v1x > df2 
      ? boost::math::ibeta(df2 / 2, df1 / 2, df2 / (df2 + v1x))
      : boost::math::ibetac(df1 / 2, df2 / 2, v1x / (df2 + v1x));
}

template <class RealType>
RealType quantile(const complemented2_type<fisher_f_distribution<RealType>, RealType>& c)
{
   RealType df1 = c.dist.degrees_of_freedom1();
   RealType df2 = c.dist.degrees_of_freedom2();
   RealType p = c.param;
   // Error check:
   RealType error_result;
   if(false == detail::check_df(
            BOOST_CURRENT_FUNCTION, df1, &error_result)
         && detail::check_df(
            BOOST_CURRENT_FUNCTION, df2, &error_result)
         && detail::check_probability(
            BOOST_CURRENT_FUNCTION, p, &error_result))
      return error_result;

   RealType x, y;

   x = boost::math::ibetac_inv(df1 / 2, df2 / 2, p, &y);

   return df2 * x / (df1 * y);
}

template <class RealType>
inline RealType mean(const fisher_f_distribution<RealType>& dist)
{ // Mean of F distribution = v.
   RealType df1 = dist.degrees_of_freedom1();
   RealType df2 = dist.degrees_of_freedom2();
   // Error check:
   RealType error_result;
   if(false == detail::check_df(
            BOOST_CURRENT_FUNCTION, df1, &error_result)
         && detail::check_df(
            BOOST_CURRENT_FUNCTION, df2, &error_result))
      return error_result;
   if(df2 <= 2)
   {
      return tools::domain_error<RealType>(
         BOOST_CURRENT_FUNCTION, "Second degree of freedom was %1% but must be > 2 in order for the distribution to have a mean.", df2);
   }
   return df2 / (df2 - 2);
} // mean

template <class RealType>
inline RealType variance(const fisher_f_distribution<RealType>& dist)
{ // Variance of F distribution.
   RealType df1 = dist.degrees_of_freedom1();
   RealType df2 = dist.degrees_of_freedom2();
   // Error check:
   RealType error_result;
   if(false == detail::check_df(
            BOOST_CURRENT_FUNCTION, df1, &error_result)
         && detail::check_df(
            BOOST_CURRENT_FUNCTION, df2, &error_result))
      return error_result;
   if(df2 <= 4)
   {
      return tools::domain_error<RealType>(
         BOOST_CURRENT_FUNCTION, "Second degree of freedom was %1% but must be > 4 in order for the distribution to have a valid variance.", df2);
   }
   return 2 * df2 * df2 * (df1 + df2 - 2) / (df1 * (df2 - 2) * (df2 - 2) * (df2 - 4));
} // variance

} // namespace math
} // namespace boost

// This include must be at the end, *after* the accessors
// for this distribution have been defined, in order to
// keep compilers that support two-phase lookup happy.
#include <boost/math/distributions/detail/derived_accessors.hpp>

#endif // BOOST_MATH_DISTRIBUTIONS_FISHER_F_HPP
