/*	Author: Jonathan Brandmeyer
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
#include "ins_qkf.hpp"
#include "assertions.hpp"
#include "timer.hpp"

#ifdef TIME_OPS
#include <iostream>
#endif

using namespace Eigen;

void
basic_ins_qkf::predict_ned(const Vector3d& gyro_meas,
		const Vector3d& accel_meas,
		double dt)
{
#ifdef TIME_OPS
	timer clock;
	clock.start();
#endif
	Vector3d accel_body = avg_state.orientation.conjugate()*accel_meas;
	// Vector3d accel_dir = accel_body.normalized();
	Matrix<double, 3, 3> accel_cov = cross(-accel_body);

	// 1500x realtime, without vectorization, on 2.2 GHz Athlon X2
	// See basic_ins_qkf::predict() for the full blockwise form
	const Matrix<double, 12, 12> pcov = cov;
	const Matrix3d dtR = dt * avg_state.orientation.conjugate().toRotationMatrix();
	const Matrix3d dtQ = accel_cov * dt;

	cov.block<3, 3>(0, 3) -= pcov.block<3,3>(0, 0)*dtR.transpose();
	cov.block<3, 3>(0, 6) += dt * pcov.block<3, 3>(0, 9);
	cov.block<3, 3>(0, 9) -= pcov.block<3, 3>(0, 3) * dtQ.transpose();
	cov.block<3, 3>(3, 3).part<Eigen::SelfAdjoint>() += dtR*pcov.block<3, 3>(0, 0)*dtR.transpose()
			- dtR*pcov.block<3, 3>(0, 3) - pcov.block<3, 3>(3, 0)*dtR.transpose();
	cov.block<3, 3>(3, 6) += -dtR * (pcov.block<3, 3>(0, 6) + dt*pcov.block<3, 3>(0, 9))
			+ dt*pcov.block<3, 3>(3, 9);
	cov.block<3, 3>(3, 9) += -dtR*( -pcov.block<3, 3>(0, 3)*dtQ.transpose() + pcov.block<3, 3>(0, 9))
			- pcov.block<3, 3>(3, 3)*dtQ.transpose();
	cov.block<3, 3>(6, 6).part<Eigen::SelfAdjoint>() += dt*pcov.block<3, 3>(6, 9) + dt*dt*pcov.block<3, 3>(9, 9)
			+ dt*pcov.block<3, 3>(9, 6);
	cov.block<3, 3>(6, 9) += -pcov.block<3, 3>(6, 3)*dtQ.transpose() + dt*pcov.block<3, 3>(9, 9)
			- dt*pcov.block<3, 3>(9, 3)*dtQ.transpose();
	cov.block<3, 3>(9, 9).part<Eigen::SelfAdjoint>() += dtQ*pcov.block<3, 3>(3, 3)*dtQ.transpose()
			- dtQ*pcov.block<3, 3>(3, 9) - pcov.block<3, 3>(9, 3)*dtQ.transpose();

	// Update symmetric cross-covariance terms
	cov.block<3, 3>(3, 0) = cov.block<3, 3>(0, 3).transpose();
	cov.block<3, 3>(6, 0) = cov.block<3, 3>(0, 6).transpose();
	cov.block<3, 3>(6, 3) = cov.block<3, 3>(3, 6).transpose();
	cov.block<3, 3>(9, 0) = cov.block<3, 3>(0, 9).transpose();
	cov.block<3, 3>(9, 3) = cov.block<3, 3>(3, 9).transpose();
	cov.block<3, 3>(9, 6) = cov.block<3, 3>(6, 9).transpose();

	// Add state transition noise
	cov.block<3, 3>(0, 0) += gyro_stability_noise.asDiagonal() * dt;
	cov.block<3, 3>(3, 3) += gyro_white_noise.asDiagonal() * dt;
	cov.block<3, 3>(6, 6) += accel_white_noise.asDiagonal() * 0.5*dt*dt;
	cov.block<3, 3>(9, 9) += accel_white_noise.asDiagonal() * dt;

	Quaterniond orientation = exp<double>((gyro_meas - avg_state.gyro_bias) * dt)
			* avg_state.orientation;
	// The contribution of the acceleration due to gravity on the accelerometer.
	// Note that in NED, the z axis is down, and the accelerometer observes
	// an extra acceleration due to gravity in the "up" direction.
	Vector3d accel_gravity = -Vector3d::UnitZ() * 9.81;
	Vector3d accel = accel_body - accel_gravity;
	Vector3d position = avg_state.position + avg_state.velocity * dt + 0.5*accel*dt*dt;
	Vector3d velocity = avg_state.velocity + accel*dt;

	avg_state.position = position;
	avg_state.velocity = velocity;
	// Note: Renormalization occurs during all measurement updates.
	avg_state.orientation = orientation;

	assert(invariants_met());
#ifdef TIME_OPS
	double time = clock.stop()*1e6;
	std::cout << "unscented predict time: " << time << "\n";
#endif
}

void
basic_ins_qkf::obs_gps_vtg_report(const Vector2d vel, const double v_error)
{
	Vector2d residual = vel - avg_state.velocity.start<2>();
#if 0
	// Compute using a single block update
	Matrix<double, 2, 2> innovation_cov = cov.block<2, 2>(9, 9)
			+ (Vector2d() << v_error, v_error).finished().asDiagonal();

	Matrix<double, 12, 2> kalman_gain = cov.block<12, 2>(0, 9)
		* innovation_cov.part<Eigen::SelfAdjoint>().inverse();

	cov.part<Eigen::SelfAdjoint>() -= kalman_gain * cov.block<2, 12>(9, 0);
	Matrix<double, 12, 1> update = kalman_gain * residual;

#else
	// Compute using rank-one updates
	Matrix<double, 12, 1> update = Matrix<double, 12, 1>::Zero();
	for (int i = 0; i < 2; ++i) {
		double innovation_cov_inv = 1.0/(cov(9+i, 9+i) + v_error);
		Matrix<double, 12, 1> gain = cov.block<12, 1>(0, 9+i) / innovation_cov_inv;
		update += gain * (residual[i] - update[9+i]);
		cov.part<Eigen::SelfAdjoint>() -= gain * cov.block<1, 12>(9+i, 0);
	}
#endif

	avg_state.apply_kalman_vec_update(update);
	assert(invariants_met());
}

