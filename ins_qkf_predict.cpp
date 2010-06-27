/*
 * ins_qkf_predict.cpp
 *
 *  Created on: Sep 2, 2009
 *      Author: Jonathan Brandmeyer
 */

#include "ins_qkf.hpp"
#include "assertions.hpp"
#include "timer.hpp"

using namespace Eigen;

namespace {

void
dump_state_errors(Matrix<double, 12, Dynamic> state_errors)
{
	// Dump only a few of the state error vectors
	Matrix<double, 6, 6> dump;
	dump << state_errors.col(0).start<6>()
		, state_errors.col(2).start<6>()
		, state_errors.col(4).start<6>()
		, state_errors.col(6).start<6>()
		, state_errors.col(8).start<6>()
		, state_errors.col(10).start<6>();
	std::cout << "state errors:\n" << dump << "\n";
}

Matrix<double, 3, 3>
axis_scale(const Vector3d& axis, double scale)
{
	return (scale - 1) * axis * axis.transpose() + Matrix3d::Identity();
}

void
linear_predict(basic_ins_qkf& _this, const Vector3d& gyro_meas, const Vector3d& accel_meas, double dt)
{
	// The two components of rotation that do not spin about the gravity vector
	// have an influence on the position and velocity of the vehicle.
	// Let r be an error axis of rotation, and z be the gravity vector.
	// Increasing r creates increasing error in the direction _|_ to r and z.
	// By the small angle theorem, the amount of error is ~ abs(r)*abs(z).
	// Increasing r also creates increasing error in the direction || to -z.
	// By the small angle theorem, the amount of error is ~ zero.
	// Therefore, rotate the error block about the z axis by -90 degrees, and
	// zero out the error vector in the z direction.
	// accel_cov is the relationship between error vectors in the tangent space
	// of the vehicle orientation and the translational reference frame.
	Vector3d accel_body = _this.avg_state.orientation.conjugate()*accel_meas;
	Vector3d accel_gravity = _this.avg_state.position.normalized()*9.81;
	Vector3d accel_resid = accel_body - accel_gravity;
#if 0
	// This form works well with zero static acceleration.
	Matrix<double, 3, 3> accel_cov =
		Eigen::AngleAxisd(-M_PI*0.5, _this.avg_state.position.normalized())
		* axis_scale(_this.avg_state.position.normalized(), 0) * 9.81;
#elif 1
	Matrix<double, 3, 3> accel_cov =
		Eigen::AngleAxisd(-M_PI*0.5, accel_body.normalized())
		* axis_scale(accel_body.normalized(), 0) * accel_meas.norm();
#else
	// The following form ends up being identical to the simpler one
	// above
	Matrix<double, 3, 3> accel_cov = 
		Eigen::AngleAxisd(-M_PI*0.5, _this.avg_state.position.normalized())
		* axis_scale(_this.avg_state.position.normalized(), 0) * 9.81
		+ Eigen::AngleAxisd(-M_PI*0.5, accel_resid.normalized())
		* axis_scale(accel_resid.normalized(), 0)*accel_resid.norm();
#endif
	// TODO: Optimization opportunity: the accel_cov doesn't change much over
	// the life of a mission. Precompute it once and then retain the original.
	// Then, only one 3x3 block ever gets updated in the A matrix below.

	// The linearized Kalman state projection matrix.
#if 0
	Matrix<double, 12, 12> A;
	     // gyro bias row
	A << Matrix<double, 3, 3>::Identity(), Matrix<double, 3, 9>::Zero(),
		 // Orientation row
		 _this.avg_state.orientation.conjugate().toRotationMatrix()*-dt,
			 Matrix<double, 3, 3>::Identity(), Matrix<double, 3, 6>::Zero(),
		 // Position row
		 Matrix<double, 3, 3>::Zero(), -accel_cov*0.5*dt*dt,
			 Matrix<double, 3, 3>::Identity(), Matrix<double, 3, 3>::Identity()*dt,
		 // Velocity row
		 Matrix<double, 3, 3>::Zero(), -accel_cov * dt,
			 Matrix<double, 3, 3>::Zero(), Matrix<double, 3, 3>::Identity();

	// 800x realtime, with vectorization
	_this.cov.part<Eigen::SelfAdjoint>() = A * _this.cov * A.transpose();
#else
	// 1500x realtime, without vectorization, on 2.2 GHz Athlon X2
	const Matrix<double, 12, 12> cov = _this.cov;
	const Matrix3d dtR = dt * _this.avg_state.orientation.conjugate().toRotationMatrix();
	const Matrix3d dtQ = accel_cov * dt;

	_this.cov.block<3, 3>(0, 3) -= cov.block<3,3>(0, 0)*dtR.transpose();
	_this.cov.block<3, 3>(0, 6) += dt * cov.block<3, 3>(0, 9);
	_this.cov.block<3, 3>(0, 9) -= cov.block<3, 3>(0, 3) * dtQ.transpose();
	_this.cov.block<3, 3>(3, 3).part<Eigen::SelfAdjoint>() += dtR*cov.block<3, 3>(0, 0)*dtR.transpose()
			- dtR*cov.block<3, 3>(0, 3) - cov.block<3, 3>(3, 0)*dtR.transpose();
	_this.cov.block<3, 3>(3, 6) += -dtR * (cov.block<3, 3>(0, 6) + dt*cov.block<3, 3>(0, 9))
			+ dt*cov.block<3, 3>(3, 9);
	_this.cov.block<3, 3>(3, 9) += -dtR*( -cov.block<3, 3>(0, 3)*dtQ.transpose() + cov.block<3, 3>(0, 9))
			- cov.block<3, 3>(3, 3)*dtQ.transpose();
	_this.cov.block<3, 3>(6, 6).part<Eigen::SelfAdjoint>() += dt*cov.block<3, 3>(6, 9) + dt*dt*cov.block<3, 3>(9, 9)
			+ dt*cov.block<3, 3>(9, 6);
	_this.cov.block<3, 3>(6, 9) += -cov.block<3, 3>(6, 3)*dtQ.transpose() + dt*cov.block<3, 3>(9, 9)
			- dt*cov.block<3, 3>(9, 3)*dtQ.transpose();
	_this.cov.block<3, 3>(9, 9).part<Eigen::SelfAdjoint>() += dtQ*cov.block<3, 3>(3, 3)*dtQ.transpose()
			- dtQ*cov.block<3, 3>(3, 9) - cov.block<3, 3>(9, 3)*dtQ.transpose();

	_this.cov.block<3, 3>(3, 0) = _this.cov.block<3, 3>(0, 3).transpose();
	_this.cov.block<3, 3>(6, 0) = _this.cov.block<3, 3>(0, 6).transpose();
	_this.cov.block<3, 3>(6, 3) = _this.cov.block<3, 3>(3, 6).transpose();
	_this.cov.block<3, 3>(9, 0) = _this.cov.block<3, 3>(0, 9).transpose();
	_this.cov.block<3, 3>(9, 3) = _this.cov.block<3, 3>(3, 9).transpose();
	_this.cov.block<3, 3>(9, 6) = _this.cov.block<3, 3>(6, 9).transpose();
#endif

	_this.cov.block<3, 3>(0, 0) += _this.gyro_stability_noise.asDiagonal() * dt;
	_this.cov.block<3, 3>(3, 3) += _this.gyro_white_noise.asDiagonal() * dt;
	_this.cov.block<3, 3>(6, 6) += _this.accel_white_noise.asDiagonal() * 0.5*dt*dt;
	_this.cov.block<3, 3>(9, 9) += _this.accel_white_noise.asDiagonal() * dt;

	Quaterniond orientation = exp<double>((gyro_meas - _this.avg_state.gyro_bias) * dt)
			* _this.avg_state.orientation;
	Vector3d accel = accel_body - _this.avg_state.position.normalized() * 9.81;
	Vector3d position = _this.avg_state.position + _this.avg_state.velocity * dt + 0.5*accel*dt*dt;
	Vector3d velocity = _this.avg_state.velocity + accel*dt;

	_this.avg_state.position = position;
	_this.avg_state.velocity = velocity;
	_this.avg_state.orientation = orientation;
}

}


void
basic_ins_qkf::predict(const Vector3d& gyro_meas, const Vector3d& accel_meas, double dt)
{
#ifdef TIME_OPS
	timer clock;
	clock.start();
#endif
#if 1
	// Always use linearized prediction
	linear_predict(*this, gyro_meas, accel_meas, dt);
#elif 0
	// Use linearized prediction only after UKF methods
	if (cov.block<3, 3>(3, 3).diagonal().maxCoeff() < 0.25*M_PI*M_PI) {
		// The following provides a ~5x speedup in the filter as a whole.
		linear_predict(*this, gyro_meas, accel_meas, dt);
#ifdef TIME_OPS
		double time = clock.stop()*1e6;
		std::cout << "linear predict time: " << time << "\n";
#endif
		return;
	}
#else

	process_state::container_t points;
	decompose_sigma_points(points, gyro_meas, accel_meas, dt);
	state::container_t projected_points;
	projected_points.reserve(points.size());
	for (auto i = points.begin(), end = points.end(); i != end; ++i) {
		projected_points.push_back(project_sigma_point(*i, dt));
	}
	state mean = average_sigma_points(projected_points);
#if 0
	// Test to ensure that what is computed by the arithmetic mean matches
	// what is computed by projecting the true mean forward through the state
	// equations.
	std::cout << "t-1 mean: "; avg_state.print(std::cout);
	std::cout << "\na priori mean: "; mean.print(std::cout);
	std::cout << "\nexpected angular change: " << gyro_meas.norm()*dt;
	std::cout << "\nabs angular change: " << avg_state.orientation.angularDistance(mean.orientation);
	std::cout << "\nexpected vec angular change: " << gyro_meas.transpose() * dt;
	std::cout << "\nalt angular change: " << (avg_state.orientation.conjugate()*-(gyro_meas*dt)).transpose();
	std::cout << "\nvec angular change: " << log<double>((mean.orientation.conjugate()*avg_state.orientation).normalized()).transpose() << "\n";

//	state alt_mean = project_sigma_point(process_state(avg_state, gyro_meas, accel_meas), dt);
//	std::cout << "\nprojected sigma point mean: "; alt_mean.print(std::cout); std::cout << "\n";
#endif
	Matrix<double, 12, 42> state_errors;
	for (int i = 0, end = points.size(); i != end; ++i) {
		state_errors.col(i) = sigma_point_difference(mean, projected_points[i]);
	}
	assert(!hasNaN(state_errors));

	// dump_state_errors(state_errors);
	cov = (state_errors * state_errors.transpose()) * 0.5;
	// std::cout << "cov: " << cov << std::endl;
	avg_state = mean;
	assert(is_real());
#endif
#ifdef TIME_OPS
	double time = clock.stop()*1e6;
	std::cout << "unscented predict time: " << time << "\n";
#endif
}

#if 0
void
basic_ins_qkf::decompose_sigma_points(
		basic_ins_qkf::process_state::container_t& points,
		const Vector3d& gyro_meas,
		const Vector3d& accel_meas,
		double)
{
	// Perform the decomposition in double-precision due to the deliberately poor
	// conditioning on the matrix.
	// TODO: Factor the covariance into D*cov, where D is a diagonal
	// scaling matrix that improves the conditioning of cov.  Then
	// multiply the result by sqrt(D).  This way the decomposition can be
	// stably performed using single-precision math.
	// cov <- 1/D*cov
	// LLT <- sqrt(D)*LLT(cov)

	// WARNING: none of the sigma points are multiplied by their weights. The
	// assumption is that the nonlinear propagation function follows
	// the relation:
	// lambda*x_1 = f(lambda*x_0)
	// for arbitrary scalars lambda.
	points.clear();
	points.reserve(42);
#if 0
	Matrix<double, 21, 21> aug_cov;
	aug_cov << cov,
		Matrix<double, 12, 9>::Zero(),
		Matrix<double, 3, 12>::Zero(), gyro_stability_noise.asDiagonal(), Matrix<double, 3, 6>::Zero(),
		Matrix<double, 3, 15>::Zero(), gyro_white_noise.asDiagonal(), Matrix<double, 3, 3>::Zero(),
		Matrix<double, 3, 18>::Zero(), accel_white_noise.asDiagonal();
	LLT<Matrix<double, 21, 21> > cov_sqrt(aug_cov);

	// std::cout << "Augmented covariance: " << aug_cov << std::endl;
	// std::cout << "LLT: " << cov_sqrt.matrixL() << std::endl;
	assert(!hasNaN(cov_sqrt.matrixL()));
	assert(!hasInf(cov_sqrt.matrixL()));
	for (int i = 0; i < 21; ++i) {
		points.push_back(process_state(avg_state, Matrix<double, 21, 1>(
				cov_sqrt.matrixL().col(i)), gyro_meas, accel_meas));
		// points.back().print(std::cout); std::cout << "\n";
		points.push_back(process_state(avg_state, Matrix<double, 21, 1>(
				-cov_sqrt.matrixL().col(i)), gyro_meas, accel_meas));
		// points.back().print(std::cout); std::cout << std::endl;
	}
#else
	// An optimization.  Since the lower 9 rows of the augmented covariance
	// matrix are zero (except for the diagonal elemnts towards the right)
	// only perform the LLT over the upper-left block.  This saves many
	// multiplications by zero.
	LLT<Matrix<double, 12, 12> > cov_sqrt(cov);
	assert(!hasNaN(cov_sqrt.matrixL()));
	assert(!hasInf(cov_sqrt.matrixL()));

	for (int i = 0; i < 12; ++i) {
		points.push_back(process_state(avg_state, Matrix<double, 12, 1>(
				cov_sqrt.matrixL().col(i)), gyro_meas, accel_meas, false));
		// points.back().print(std::cout); std::cout << "\n";
		points.push_back(process_state(avg_state, Matrix<double, 12, 1>(
				-cov_sqrt.matrixL().col(i)), gyro_meas, accel_meas, false));
		// points.back().print(std::cout); std::cout << std::endl;
	}

	for (int i = 0; i < 3; ++i) {
		double v = std::sqrt(gyro_stability_noise[i]);
		points.push_back(process_state(avg_state, gyro_meas, accel_meas));
		points.back().gyro_instability[i] = v;
		points.push_back(process_state(avg_state, gyro_meas, accel_meas));
		points.back().gyro_instability[i] = -v;
	}

	for (int i = 0; i < 3; ++i) {
		double v = std::sqrt(gyro_white_noise[i]);
		points.push_back(process_state(avg_state, gyro_meas, accel_meas));
		points.back().gyro_meas[i] += v;
		points.push_back(process_state(avg_state, gyro_meas, accel_meas));
		points.back().gyro_meas[i] += -v;
	}

	for (int i = 0; i < 3; ++i) {
		double v = std::sqrt(accel_white_noise[i]);
		points.push_back(process_state(avg_state, gyro_meas, accel_meas));
		points.back().accel_meas[i] += v;
		points.push_back(process_state(avg_state, gyro_meas, accel_meas));
		points.back().accel_meas[i] += -v;
	}
#endif
}


basic_ins_qkf::state
basic_ins_qkf::project_sigma_point(
	const basic_ins_qkf::process_state& p,
	double dt)
{
	// Optimization opportunity:
	// TODO: Much of the cost (~45%) in this algorithm is associated with taking the
	// exponential of the model rate.  However, only 6 of the 42 sigma points
	// have a model rate that is affected by the gyro bias.
	state ret;
	ret.gyro_bias = p.gyro_bias + p.gyro_instability * dt;
	
	Vector3d model_rate = (p.gyro_meas - p.gyro_bias)*dt;
#if 0
	Vector3d rev_model_rate = log<double>(exp<double>(model_rate));
	if (!model_rate.isApprox(rev_model_rate)) {
		std::cout << "Warning: " << model_rate.transpose()
				<< " round-trips to " << rev_model_rate.transpose() << "\n";
	}
#endif

	Quaterniond step(exp<double>(model_rate));
	ret.orientation = (step * p.orientation);

	// TODO: Refine the accelleration due to gravity by the ECEF distance, or some other model
	Vector3d accel = p.orientation.conjugate() * p.accel_meas -
			p.position.normalized() * 9.81;
	ret.velocity = (p.velocity + accel * dt);
	ret.position = p.position + p.velocity*dt + (0.5 * accel*dt*dt);
	assert(ret.is_real());
	return ret;
}

#endif
