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

basic_ins_qkf::basic_ins_qkf(
		const Vector3d& estimate,
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
	assert(is_real());
}

void
basic_ins_qkf::counter_rotate_cov(const Quaterniond&)
{
	// Rotate the principle axes of the angular error covariance by the
	// mean update.

	// TODO: This is only required in the case that the system covariance is
	// right-multiplied by the mean. The current design left-multiplies the
	// covariance by the mean.
	return;

	// TODO: There should be an expression that makes both this matrix's
	// construction and multiplication much more efficient.
	// Matrix<double, 12, 12> counter_rot = Matrix<double, 12, 12>::Identity();
	// counter_rot.block<3, 3>(3, 3) = update.cast<double>().conjugate().toRotationMatrix();
	// cov = (counter_rot * cov.cast<double>() * counter_rot.transpose()).cast<float>();
}
#if 0
basic_ins_qkf::state
basic_ins_qkf::average_sigma_points(const std::vector<state, aligned_allocator<state> >& points)
{
	state ret;
	Vector3d sum;
#define avg_vector_field(field) \
	sum = Vector3d::Zero(); \
	for (auto i = points.begin(), end = points.end(); i != end; ++i) { \
		sum += i->field; \
	} \
	ret.field = sum / points.size()

	avg_vector_field(gyro_bias);
	avg_vector_field(velocity);

	Vector3d p_sum = Vector3d::Zero();
	for (auto i = points.begin(), end = points.end(); i != end; ++i) {
		p_sum += i->position;
	}
	ret.position = p_sum / (double)points.size();

	std::vector<Quaterniond> quat_points;
	quat_points.reserve(points.size());
	for (auto i = points.begin(), end = points.end(); i != end; ++i) {
		quat_points.push_back(i->orientation);
	}
	ret.orientation = quaternion_avg_johnson(quat_points).normalized();
	return ret;
}
#endif

bool
basic_ins_qkf::is_real(void) const
{
	return !(hasNaN(cov) || hasInf(cov)) && avg_state.is_real();
}

Quaterniond
basic_ins_qkf::state::apply_kalman_vec_update(const Matrix<double, 12, 1> update)
{
	// std::cout << "***update available***\n"
	// 		<< "\tstate: "; print(std::cout);
	// std::cout << "\n\tupdate: " << update.transpose() << "\n";
	gyro_bias += update.segment<3>(0);
	Quaterniond posterior_update = exp<double>(update.segment<3>(3));
	orientation = (orientation * posterior_update).normalized();
	position += update.segment<3>(6);
	velocity += update.segment<3>(9);
	assert(is_real());
	return posterior_update;
}

#if 0
Quaterniond
basic_ins_qkf::state::apply_left_kalman_vec_update(const Matrix<double, 12, 1> update)
{
	// std::cout << "***update available***\n"
	// 		<< "\tstate: "; print(std::cout);
	// std::cout << "\n\tupdate: " << update.transpose() << "\n";
	gyro_bias += update.segment<3>(0);
	Quaterniond posterior_update = exp<double>(update.segment<3>(3));
	orientation = (posterior_update * orientation).normalized();
	position += update.segment<3>(6);
	velocity += update.segment<3>(9);
	assert(is_real());
	return posterior_update;
}
#endif

bool
basic_ins_qkf::state::has_nan(void)const
{
	return hasNaN(gyro_bias) || hasNaN(orientation.coeffs())
			|| hasNaN(position) || hasNaN(velocity);
}

bool
basic_ins_qkf::state::is_real(void) const
{
	return !(hasNaN(gyro_bias) || hasNaN(orientation.coeffs())
			|| hasNaN(position) || hasNaN(velocity))
			&& !(hasInf(gyro_bias) || hasInf(orientation.coeffs())
				|| hasInf(position) || hasInf(velocity));
}

void
basic_ins_qkf::state::print(std::ostream& str)
{
	str << "gyro_bias: " << gyro_bias.transpose()
		<< " orientation: " << orientation.coeffs().transpose()
		<< " position: " << position.transpose()
		<< " velocity: " << velocity.transpose();
}


double
basic_ins_qkf::angular_error(const Quaterniond& q) const
{
	return q.angularDistance(avg_state.orientation);
}

double
basic_ins_qkf::mahalanobis_distance(const state& q) const
{
	state_error_t delta = sigma_point_difference(avg_state, q);
	state_error_t inv_delta;
	cov.lu().solve(delta, &inv_delta);
	return std::sqrt(delta.dot(inv_delta));
}


basic_ins_qkf::state_error_t
basic_ins_qkf::sigma_point_difference(
	const basic_ins_qkf::state& median,
	const basic_ins_qkf::state& point) const
{
	// TODO: Optimization opportunity.  > 80% of the cost of this method is consumed
	// by taking the quaternion log.  However, many of these logs are taken of the identity.
	// Come up with some logic to make that path faster.
	state_error_t ret;
	ret.segment<3>(0) = point.gyro_bias - median.gyro_bias;
	if (median.orientation.coeffs().dot(point.orientation.coeffs()) < 0) {
		// q == -q, but the covariance relation doesn't work without this step.
		// Force the point to lie on the same hemisphere as the mean.
		Quaterniond neg_orientation(point.orientation);
		neg_orientation.coeffs() *= -1;
		ret.segment<3>(3) = log<double>(median.orientation.conjugate() * neg_orientation);
	}
	else {
		ret.segment<3>(3) = log<double>(median.orientation.conjugate() * point.orientation);
	}
	ret.segment<3>(6) = point.position - median.position;
	ret.segment<3>(9) = point.velocity - median.velocity;
	assert(!hasNaN(ret));
	return ret;
}
