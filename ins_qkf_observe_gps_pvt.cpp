/*
 * ins_qkf_observe_gps_pvt.cpp
 *
 *  Created on: Sep 2, 2009
 *      Author: Jonathan Brandmeyer
 */

#include "ins_qkf.hpp"
#include "assertions.hpp"
#include <Eigen/LU>
#include "timer.hpp"

using namespace Eigen;

#define RANK_ONE_UPDATES

void
basic_ins_qkf::obs_gps_v_report(const Vector3d& vel, const Vector3d& v_error)
{
	Vector3d residual = vel - avg_state.velocity;
	Matrix<double, 3, 3> innovation_cov = cov.block<3, 3>(9, 9);
	innovation_cov += v_error.asDiagonal();

#if 0
	// WARNING: Running rank-one updates in velocity tends to destabilize
	// the filter during initialization.
	Matrix<double, 12, 1> update = Matrix<double, 12, 1>::Zero();
	for (int i = 0; i < 3; ++i) {
		Matrix<double, 12, 1> gain = cov.block<12, 1>(0, 9+i) / innovation_cov(i, i);
		update += gain * (residual[i] - update[9+i]);
		cov -= gain * cov.block<1, 12>(9+i, 0);
	}
#else
	Matrix<double, 3, 12> kalman_gain_t;
	innovation_cov.qr().solve(cov.block<3, 12>(9, 0), &kalman_gain_t);
	cov.part<Eigen::SelfAdjoint>() -= cov.block<12, 3>(0, 9) * kalman_gain_t; // .transpose() * cov.block<3, 12>(9, 0);
	Matrix<double, 12, 1> update = kalman_gain_t.transpose() * residual;
#endif

	Quaterniond rotor = avg_state.apply_kalman_vec_update(update);
	counter_rotate_cov(rotor);
	assert(is_real());
}

void
basic_ins_qkf::obs_gps_pv_report(
		const Vector3d& pos, const Vector3d& vel,
		const Vector3d& p_error, const Vector3d v_error)
{
#ifdef TIME_OPS
	timer clock;
	clock.start();
#endif
#if 1
	obs_gps_p_report(pos, p_error);
	obs_gps_v_report(vel, v_error);
#else

	// The observation model is strictly linear here, so use the linear
	// form of the kalman gain and update
	Matrix<double, 6, 12> obs_matrix;
	obs_matrix << Matrix<double, 6, 6>::Zero(), Matrix<double, 6, 6>::Identity();

	Matrix<double, 6, 1> residual;
	residual.segment<3>(0) = (pos - avg_state.position);
	residual.segment<3>(3) = vel - avg_state.velocity;

	Matrix<double, 6, 6> innovation_cov = cov.corner<6, 6>(Eigen::BottomRight);
	//innovation_cov = obs_matrix * cov * obs_matrix.transpose();
	innovation_cov.corner<3, 3>(Eigen::TopLeft) += p_error.asDiagonal();
	innovation_cov.corner<3, 3>(Eigen::BottomRight) += v_error.asDiagonal();
#if 1
	// Perform matrix inverse by LU decomposition instead of cofactor expansion.
	// K = P*transpose(H)*inverse(S)
	// K = P*transpose(transpose(transpose(H)*inverse(S)))
	// K = P*transpose(transpose(inverse(S))*H)
	// K = P*transpose(inverse(transpose(S))*H)
	// S == transpose(S)
	// K = P*transpose(inverse(S)*H)
	// obs_matrx <- inverse(S)*H
	Matrix<double, 6, 12> inv_s_h;
	innovation_cov.qr().solve(obs_matrix, &inv_s_h);
	Matrix<double, 12, 6> kalman_gain = cov * inv_s_h.transpose();
#else
	Matrix<double, 12, 6> kalman_gain = cov.block<12, 6>(0, 6) * innovation_cov.inverse();
#endif
	Quaterniond rotor = avg_state.apply_kalman_vec_update(kalman_gain * residual);
	cov.part<Eigen::SelfAdjoint>() -= kalman_gain * obs_matrix * cov;
	counter_rotate_cov(rotor);
	assert(is_real());
#endif

#ifdef TIME_OPS
	double time = clock.stop() * 1e6;
	std::cout << "obs_gps_pvt time: " << time << "\n";
#endif
}

