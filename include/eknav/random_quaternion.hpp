#ifndef LIBEKNAV_RANDOM_QUATERNION_HPP
#define LIBEKNAV_RANDOM_QUATERNION_HPP
/*
 * random_quaternion.hpp
 *
 *      Author: Jonathan Brandmeyer
 *
 *          This file is part of libeknav.
 *
 *  Libeknav is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, version 3.
 *
 *  Libeknav is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.

 *  You should have received a copy of the GNU General Public License
 *  along with libeknav.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "eknav/random_vector.hpp"
#include "eknav/quaternions.hpp"

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
