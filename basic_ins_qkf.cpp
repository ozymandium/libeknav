/*
 * basic_ins_qkf.cpp
 *
 *  Created on: Aug 11, 2009
 *      Author: Jonathan Brandmeyer
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
 *
 */

#include "ins_qkf.hpp"
#include "assertions.hpp"

using namespace Eigen;

basic_ins_qkf::basic_ins_qkf(const Vector3d& estimate,
		double pos_error, double bias_error, double v_error,
		const Vector3d& gyro_white_noise,
		const Vector3d& gyro_stability_noise,
		const Vector3d& accel_white_noise,
		const Vector3d& vel_estimate)
	: gyro_stability_noise(gyro_stability_noise)
	, gyro_white_noise(gyro_white_noise)
	, accel_white_noise(accel_white_noise)
{
	avg_state.position = estimate;
	avg_state.gyro_bias = Vector3d::Zero();
	avg_state.orientation = Quaterniond::Identity();
	avg_state.velocity = vel_estimate;

	cov << Matrix3d::Identity()*bias_error*bias_error, Matrix<double, 3, 9>::Zero(),
		Matrix3d::Zero(), Matrix3d::Identity()*M_PI*M_PI*0.5, Matrix<double, 3, 6>::Zero(),
		Matrix<double, 3, 6>::Zero(), Matrix3d::Identity()*pos_error*pos_error, Matrix3d::Zero(),
		Matrix<double, 3, 9>::Zero(), Matrix3d::Identity()*v_error*v_error;
	assert(invariants_met());
}

Quaterniond
basic_ins_qkf::state::apply_kalman_vec_update(const Matrix<double, 12, 1> update)
{
	// std::cout << "***update available***\n"
	// 		<< "\tstate: "; print(std::cout);
	// std::cout << "\n\tupdate: " << update.transpose() << "\n";
	gyro_bias += update.segment<3>(0);
	Quaterniond posterior_update = exp<double>(update.segment<3>(3));
	orientation = incremental_normalized(orientation * posterior_update);
	position += update.segment<3>(6);
	velocity += update.segment<3>(9);
	return posterior_update;
}


