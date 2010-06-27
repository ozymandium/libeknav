#ifndef AHRS_RANDOM_QUATERNION_HPP
#define AHRS_RANDOM_QUATERNION_HPP

#include "random_vector.hpp"
#include "quaternions.hpp"

template<typename FloatT>
Quaternion<FloatT>
random_quaternion(void)
{
	using namespace Eigen;
	static auto engine = boost::minstd_rand0();
	static RandomVector<FloatT, 3, boost::minstd_rand0> entropy(
			Matrix<FloatT, 3, 1>::Zero(),
			(Matrix<FloatT, 3, 1>() << M_PI, M_PI, M_PI).finished().asDiagonal(),
			&engine);

	return exp<FloatT>(entropy());
}


#endif // !defined AHRS_RANDOM_QUATERNION_HPP
