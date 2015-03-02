/*
 * diagnostics.cpp
 *
 *  Created on: Jul 3, 2010
 *      Author: jonathan
 *
 *  Member functions that are not expected to be needed in a release build
 *  are in this file.  This way, a compiler that does not place each func
 *  defn in a separate section will be able to elide these functions at
 *  link time.
 */

#include <limits>

#include "ins_qkf.hpp"
#include "assertions.hpp"

using namespace Eigen;

basic_ins_qkf::state_error_t
basic_ins_qkf::sigma_point_difference(const basic_ins_qkf::state& median,
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

void
basic_ins_qkf::state::print(std::ostream& str)
{
	str << "gyro_bias: " << gyro_bias.transpose()
		<< " orientation: " << orientation.coeffs().transpose()
		<< " position: " << position.transpose()
		<< " velocity: " << velocity.transpose();
}

double
basic_ins_qkf::mahalanobis_distance(const state& q) const
{
	state_error_t delta = sigma_point_difference(avg_state, q);
	state_error_t inv_delta = cov.lu().solve(delta);
	return std::sqrt(delta.dot(inv_delta));
}

double
basic_ins_qkf::angular_error(const Quaterniond& q) const
{
	return q.angularDistance(avg_state.orientation);
}

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

bool
basic_ins_qkf::is_real(void) const
{
	return !(hasNaN(cov) || hasInf(cov)) && avg_state.is_real();
}

bool
basic_ins_qkf::invariants_met(void) const
{
	// The whole thing breaks down if NaN or Inf starts popping up
	return is_real() &&
		// Incremental normalization is working
		std::abs(1 - 1.0/avg_state.orientation.norm()) < 
			std::sqrt(std::numeric_limits<double>::epsilon());

}

