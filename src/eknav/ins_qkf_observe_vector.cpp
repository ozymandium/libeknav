/*
 * ins_qkf_observe_vector.cpp
 *
 *  Created on: Sep 2, 2009
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
 */

#include "eknav/ins_qkf.hpp"
#include "eknav/assertions.hpp"
#include "eknav/posix/timer.hpp"

#ifdef TIME_OPS
#include <iostream>
#endif

using namespace Eigen;

void
basic_ins_qkf::obs_gyro_bias(const Vector3d& bias, const Vector3d& bias_error)
{
	Matrix<double, 3, 3> s = cov.block<3, 3>(0, 0);
	s.diagonal() += bias_error;

	Matrix<double, 12, 3> kalman_gain = cov.block<12, 3>(0, 0) * s.inverse();
	cov -= kalman_gain * cov.block<3, 12>(0, 0);
	Vector3d innovation = bias - avg_state.gyro_bias;

	// Apply the Kalman gain to obtain the posterior state and error estimates.
	avg_state.apply_kalman_vec_update(kalman_gain * innovation);
	assert(invariants_met());
}

void
basic_ins_qkf::obs_vector(const Vector3d& ref,
		const Vector3d& obs,
		double error)
{
#ifdef TIME_OPS
	timer clock;
	clock.start();
#endif
#define DEBUG_VECTOR_OBS 0

	Vector3d obs_ref = avg_state.orientation.conjugate()*obs;
	Vector3d v_residual = log<double>(Quaterniond().setFromTwoVectors(ref, obs_ref));

	Matrix<double, 3, 2> h_trans;
	h_trans.col(0) = ref.cross(
		(abs(ref.dot(obs_ref)) < 0.9994) ? obs_ref :
			(abs(ref.dot(Vector3d::UnitX())) < 0.707)
				? Vector3d::UnitX() : Vector3d::UnitY()).normalized();
	h_trans.col(1) = -ref.cross(h_trans.col(0));
	assert(!hasNaN(h_trans));
	assert(h_trans.isUnitary());

#ifdef RANK_ONE_UPDATES
	// Running a rank-one update here is a strict win.
	Matrix<double, 12, 1> update = Matrix<double, 12, 1>::Zero();
	for (int i = 0; i < 2; ++i) {
		double obs_error = error;
		double obs_cov = (h_trans.col(i).transpose() * cov.block<3, 3>(3, 3) * h_trans.col(i))[0];
		Matrix<double, 12, 1> gain = cov.block<12, 3>(0, 3) * h_trans.col(i) / (obs_error + obs_cov);
		update += gain * h_trans.col(i).transpose() * v_residual;
		// TODO: Get Eigen to treat cov as self-adjoint
		cov -= gain * h_trans.col(i).transpose() * cov.block<3, 12>(3, 0);
	}
#else
	// block-wise form.  This is much less efficient.
	Vector2d innovation = h_trans.transpose() * v_residual;
	Matrix<double, 12, 2> kalman_gain = cov.block<12, 3>(0, 3) * h_trans
			* (h_trans.transpose() * cov.block<3, 3>(3, 3) * h_trans
				+ (Vector2d() << error, error).finished().asDiagonal()).inverse();
	// TODO: Get Eigen to treat cov as self-adjoint
	cov -= kalman_gain * h_trans.transpose() * cov.block<3, 12>(3, 0);
	Matrix<double, 12, 1> update = (kalman_gain * innovation);
#endif


#if DEBUG_VECTOR_OBS
	// std::cout << "projected update: " << (obs_projection * update.segment<3>(3)).transpose() << "\n";
	std::cout << "deprojected update: " << update.segment<3>(3).transpose() << "\n";
#endif
	avg_state.apply_kalman_vec_update(update);

	assert(invariants_met());
#ifdef TIME_OPS
	double time = clock.stop() * 1e6;
	std::cout << "observe_vector(): " << time << "\n";
#endif

}
