/*
 * pr_ins_qkf.cpp
 *
 *  Created on: May 1, 2011
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

#include <pr_ins_qkf.hpp>
#include <assertions.hpp>

using namespace Eigen;

pseudorange_ins_qkf::pseudorange_ins_qkf()
{
	avg_state.position = Vector3d::Zero();
	avg_state.gyro_bias = Vector3f::Zero();
	avg_state.orientation = Quaterniond::Identity();
	avg_state.velocity = Vector3d::Zero();
	avg_state.accel_bias = Vector3f::Zero();
	avg_state.clock_bias = 0;

	cov = Matrix<float, 12, 12>::Zero();
	pt_cov = Matrix<float, 4, 4>::Zero();

	// Default error bounds
	float gyro_bias_rms = 3.0 * M_PI / 180;
	gyro_bias_rms *= gyro_bias_rms;

	// default accelerometer bias: 15 mg
	float accel_bias_rms = 15.0 * 9.81 / 1e3;
	accel_bias_rms *= accel_bias_rms;

	// default clock rms noise  1 us/sqrt(s), at 3e8 m/s
	float clock_bias_rms = 300;
	clock_bias_rms *= clock_bias_rms;

	// default position error: 100 km
	float initial_pos_err = 100e3;
	initial_pos_err *= initial_pos_err;

	cov.block<3, 3>(0, 0) = Matrix3f::Identity()*gyro_bias_rms; // gyro bias
	cov.block<3, 3>(3, 3) = Matrix3f::Identity()*M_PI*M_PI*0.5; // attitude
	cov.block<3, 3>(6, 6) = Matrix3f::Identity() * 100; // velocity
	cov.block<3, 3>(9, 9) = Matrix3f::Identity() * accel_bias_rms;

	pt_cov.block<3, 3>(0, 0) = Matrix3f::Identity() * initial_pos_err;
	pt_cov(3, 3) = clock_bias_rms;

	assert(invariants_met());
}

void
pseudorange_ins_qkf::init_attitude(const Quaterniond& attitude,
		const Matrix3f& attitude_error)
{
	avg_state.orientation = attitude;
	clear_covariance_block(3, attitude_error);
}

void
pseudorange_ins_qkf::init_velocity(const Vector3d& vel, const Vector3f& vel_error)
{
	avg_state.velocity = vel;
	clear_covariance_block(6, vel_error.asDiagonal());
}

void
pseudorange_ins_qkf::init_position(const Vector3d& pos, const Vector3f& pos_error)
{
	avg_state.position = pos;
	clear_covariance_block(12, pos_error.asDiagonal());
}

void
pseudorange_ins_qkf::clear_covariance_block(size_t rowcol,
	   const Matrix3f& repl)
{
	if (rowcol <= 9) {
		// Zero out the self-covariance of this 3x3 block, as well as all
		// cross-covariance terms associated with it.
		cov.block<3, 12>(rowcol, 0) = Matrix<float, 3, 12>::Zero();
		cov.block<12, 3>(0, rowcol) = Matrix<float, 12, 3>::Zero();
		// Initialize the self variance to the provided value.
		cov.block<3, 3>(rowcol, rowcol) = repl;
	}
	else {
		// default the position and GPS clock bias
		float clock_bias_rms = 300;

		clock_bias_rms *= clock_bias_rms;
		pt_cov = Matrix<float, 4, 4>::Zero();
		pt_cov.block<3, 3>(0, 0) = repl;
		pt_cov(3, 3) = clock_bias_rms;
	}
}

namespace {
/**
 * dst[dst_row, dst_col] += mult*src[src_row, src_col] + src[src_col, src_row]*mult.transpose()
 */
void ssyr2k(Matrix<float, 12, 12>& dst, int dst_row, int dst_col,
		const Matrix3f& mult, const Matrix<float, 12, 12>& src, int src_row, int src_col)
{
	dst.block<3, 3>(dst_row, dst_col) += mult * src.block<3, 3>(src_row, src_col)
			+ src.block<3, 3>(src_col, src_row)*mult.transpose();
}

/**
 * dst[dst, dst] += mult * src[src, src] * mult'
 * dst and src both symmetric
 */
void sgemmm(Matrix<float, 12, 12>& dst, int dst_diag,
		const Matrix3f& mult,
		const Matrix<float, 12, 12>& src, int src_diag)
{
	dst.block<3, 3>(dst_diag, dst_diag).part<Eigen::SelfAdjoint>() +=
			mult * src.block<3, 3>(src_diag, src_diag).part<Eigen::SelfAdjoint>() * mult.transpose();
}
} // !namespace anon

void
pseudorange_ins_qkf::predict_ecef(const Vector3f& gyro_meas,
		const Vector3f& accel_meas,
		float dt)
{
	Quaternionf attitude_conj = avg_state.orientation.cast<float>().conjugate();

	// Rotate the sensible acceleration into the inertial frame
	Vector3f accel_sensible_ecef = attitude_conj*(accel_meas - avg_state.accel_bias);
	// The local gravity vector, in the ECEF frame
	Vector3f accel_gravity = avg_state.position.cast<float>().normalized() * accel_gravity_norm;

	avg_state.inertial_accel = accel_sensible_ecef - accel_gravity;

	Matrix<float, 3, 3> accel_cov = cross<float>(-accel_sensible_ecef);

	// Take a full copy of the covariance matrix
	const Matrix<float, 12, 12> cov = this->cov;
	const Matrix<float, 4, 4> pt_cov = this->pt_cov;

	// Some convenience blocks
	const Matrix3f dtR = dt * attitude_conj.toRotationMatrix();
	const Matrix3f dtQ = dt * accel_cov;

	// For space savings, the matrix updates fall into the following basic forms:
	// C += A*B*A'     // where B == B'
	// C += A*B + B'A' // SSYR2K
	// C += A * B'     // SGEMM
	// C += A * B      // SGEMM (no transpose option)
	// They can be refactored to use the above instead to save some flash space

	// Compute the next covariance matrix using sparse blockwise operations.
	// TODO: Consider factoring this out into primitive BLAS-like ops to save space
	// TODO: This part certainly isn't right, thanks to changes in addressing
	// nop
	// this->cov.block<3, 3>(0, 0) = this->cov.block<3, 3>(0, 0);
	this->cov.block<3, 3>(0, 3) -= cov.block<3,3>(0, 0)*dtR.transpose();
	this->cov.block<3, 3>(0, 6) -= cov.block<3, 3>(0, 3)*dtQ.transpose() + cov.block<3, 3>(0, 9)*dtR.transpose();
	// nop
	// this->cov.block<3, 3>(0, 9) = this->cov.block<3, 3>(0, 9);

	sgemmm(this->cov, 3, dtR, cov, 0);
	ssyr2k(this->cov, 3, 3, -dtR, cov, 0, 3);

	this->cov.block<3, 3>(3, 6) += -dtR * cov.block<3, 3>(0, 6) - cov.block<3, 3>(3, 3)*dtQ.transpose()
			+ dtR*cov.block<3, 3>(0, 3)*dtQ.transpose() - cov.block<3, 3>(3, 9)*dtR.transpose()
			+ dtR*cov.block<3, 3>(0, 9)*dtR.transpose();
	this->cov.block<3, 3>(3, 9) -= dtR*cov.block<3, 3>(0, 9);

	ssyr2k(this->cov, 6, 6, -dtQ, cov, 3, 6);
	ssyr2k(this->cov, 6, 6, -dtR, cov, 9, 6);
	{
		Matrix3f tmp = dtR * (cov.block<3, 3>(3, 9) * dtQ).transpose();
		this->cov.block<3, 3>(6, 6).part<Eigen::SelfAdjoint>() += tmp + tmp.transpose();
	}
	sgemmm(this->cov, 6, dtQ, cov, 3);
	sgemmm(this->cov, 6, dtR, cov, 9);
	this->cov.block<3, 3>(6, 9) -= dtQ*cov.block<3, 3>(3, 9) + dtR*cov.block<3, 3>(9, 9);

	// nop
	// this->cov.block<3, 3>(9, 9) = this->cov.block<3, 3>(9, 9);

	this->pt_cov.block<3, 3>(0, 0).part<Eigen::SelfAdjoint>() += dt*dt*cov.block<3, 3>(6, 6);

	// Maintain symmetric form
	struct blockaddr_t {
		int row, col;
	} block_addr[] = {
			{ 3, 0 },
			{ 6, 0 },
			{ 6, 3 },
			{ 9, 0 },
			{ 9, 3 },
			{ 9, 6 },
	};
	for (int i = 0; i < 6; ++i ) {
		int row = block_addr[i].row;
		int col = block_addr[i].col;
		this->cov.block<3, 3>(row, col) = this->cov.block<3, 3>(col, row).transpose();
	}

	// Add Q-matrix state estimate noise blocks
	this->cov.block<3, 3>(0, 0) += gyro_stability_noise.asDiagonal() * dt;
	this->cov.block<3, 3>(3, 3) += gyro_white_noise.asDiagonal() * dt;
	this->cov.block<3, 3>(6, 6) += accel_white_noise.asDiagonal() * dt;
	this->cov.block<3, 3>(9, 9) += accel_stability_noise.asDiagonal() * dt;
	this->pt_cov.block<3, 3>(0, 0) += accel_white_noise.asDiagonal() * 0.5*dt*dt;
	this->pt_cov(3, 3) += clock_stability_noise * dt;

	// Project the mean forward
	Vector3d accel = avg_state.inertial_accel.cast<double>();
	avg_state.body_rate = gyro_meas - avg_state.gyro_bias;
	Quaterniond orientation = exp<float>(avg_state.body_rate * dt).cast<double>()
			* avg_state.orientation;
	Vector3d position = avg_state.position + avg_state.velocity * dt + 0.5*accel*dt*dt;
	Vector3d velocity = avg_state.velocity + accel*dt;

	avg_state.position = position;
	avg_state.velocity = velocity;
	// Note: Renormalization occurs during all measurement updates.
	avg_state.orientation = orientation;


	assert(invariants_met());
}

bool
pseudorange_ins_qkf::state::has_nan(void)const
{
	return hasNaN(gyro_bias)
		|| hasNaN(orientation.coeffs())
		|| hasNaN(position)
		|| hasNaN(velocity)
		|| hasNaN(accel_bias)
		|| hasNaN(inertial_accel)
		|| hasNaN(body_rate)
		|| std::isnan(clock_bias);
}

bool
pseudorange_ins_qkf::state::has_inf(void)const
{
	return hasInf(gyro_bias)
		|| hasInf(orientation.coeffs())
		|| hasInf(position)
		|| hasInf(velocity)
		|| hasInf(accel_bias)
		|| hasInf(inertial_accel)
		|| hasInf(body_rate)
		|| std::isinf(clock_bias);
}

bool
pseudorange_ins_qkf::state::is_real(void) const
{
	return !has_nan() && !has_inf();
}

Quaterniond
pseudorange_ins_qkf::state::apply_kalman_vec_update(const Matrix<float, 12, 1>& update)
{
	gyro_bias += update.segment<3>(0);
	Quaterniond posterior_update = exp<float>(update.segment<3>(3)).cast<double>();
	orientation = incremental_normalized(orientation * posterior_update);
	velocity += update.segment<3>(6).cast<double>();
	accel_bias += update.segment<3>(9);
	return posterior_update;
}

void
pseudorange_ins_qkf::state::apply_kalman_vec_update(const Matrix<float, 4, 1>& update)
{
	position += update.segment<3>(0).cast<double>();
	clock_bias += update(3);
}

float
pseudorange_ins_qkf::angular_error(const Quaterniond& q) const
{
	return q.angularDistance(avg_state.orientation);
}

float
pseudorange_ins_qkf::gyro_bias_error(const Vector3f& gyro_bias) const
{
	return (avg_state.gyro_bias - gyro_bias).norm();
}

float
pseudorange_ins_qkf::accel_bias_error(const Vector3f& accel_bias) const
{
	return (avg_state.accel_bias - accel_bias).norm();
}

float
pseudorange_ins_qkf::mahalanobis_distance(const state& q) const
{
	state_error_t delta = sigma_point_difference(avg_state, q);

	Matrix<float, 12, 1> main_err = delta.start<12>(), inv_delta;
	Matrix<float, 4, 1> pos_err = delta.end<4>(), inv_delta_end;

	cov.lu().solve(delta.start<12>(), &inv_delta);
	pt_cov.lu().solve(delta.end<4>(), &inv_delta_end);
	return std::sqrt(main_err.dot(inv_delta) + pos_err.dot(inv_delta_end));
}

pseudorange_ins_qkf::state_error_t
pseudorange_ins_qkf::sigma_point_difference(
		const pseudorange_ins_qkf::state& mean,
		const pseudorange_ins_qkf::state& point) const
{
	state_error_t ret;
	ret.segment<3>(0) = point.gyro_bias - mean.gyro_bias;
	if (mean.orientation.coeffs().dot(point.orientation.coeffs()) < 0) {
		// q == -q, but the covariance relation doesn't work without this step.
		// Force the point to lie on the same hemisphere as the mean.
		Quaterniond neg_orientation(point.orientation);
		neg_orientation.coeffs() *= -1;
		ret.segment<3>(3) = log<double>(mean.orientation.conjugate() * neg_orientation).cast<float>();
	}
	else {
		ret.segment<3>(3) = log<double>(mean.orientation.conjugate() * point.orientation).cast<float>();
	}
	ret.segment<3>(6) = (point.velocity - mean.velocity).cast<float>();
	ret.segment<3>(9) = point.accel_bias - mean.accel_bias;
	ret.segment<3>(12) = (point.position - mean.position).cast<float>();
	ret(15) = point.clock_bias - mean.clock_bias;

	assert(!hasNaN(ret));
	return ret;
}


bool
pseudorange_ins_qkf::invariants_met(void) const
{
	// The whole thing breaks down if NaN or Inf starts popping up
	return is_real() &&
		// Incremental normalization is working
		std::abs(1 - 1.0/avg_state.orientation.norm()) < std::sqrt(Eigen::machine_epsilon<double>());
}

bool
pseudorange_ins_qkf::is_real(void) const
{
	return !(hasNaN(cov)       ||
				hasInf(cov)    ||
				hasNaN(pt_cov) ||
				hasInf(pt_cov)) &&
			avg_state.is_real();
}
