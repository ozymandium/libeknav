#ifndef LIBEKNAV_PR_INS_QKF_HPP
#define LIBEKNAV_PR_INS_QKF_HPP
/*
 * pr_ins_qkf.hpp
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

#include <sigma_points.hpp>
#include <quaternions.hpp>
#include <Eigen/StdVector>

using Eigen::Vector3f;
using Eigen::Vector3d;
using Eigen::Vector2d;
using Eigen::Quaterniond;
using Eigen::Matrix;
using Eigen::aligned_allocator;


/**
 * A partially split GPS/INS integration class.  The covariance of the position estiamte
 * is forced to be separate from the covariance of the rest of the state.  In effect,
 * the position estimate itself is a separate filter.
 *
 * This class does not use position estimates to correct the terms in the velocity
 * solution, or to estimate the bias in the accelerometer.  Four satellite observations
 * per period are absolutely required, in both deltarange by doppler and pseudorange
 * by timing.
 *
 * Floats are used for error estimates, while doubles are used for the state vector
 * itself.
 */
struct pseudorange_ins_qkf
{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	/**
	 * The covariance of the zero-mean gaussian white noise that is added to
	 * the gyro bias at each time step.  This value is treated as a diagonal
	 *  matrix, in units of radians^2/second^2.
	 */
	Vector3f gyro_stability_noise;

	/**
	 * The covariance of the zero-mean gaussian white noise that is added to
	 * the gyro measurement at each time step, in rad^2/second^2
	 */
	Vector3f gyro_white_noise;

	/**
	 * The covariance of the zero-mean gaussian white noise that is added to
	 * the accelerometer measurement at each time step, in (m/s/s)^2
	 */
	Vector3f accel_white_noise;

	/**
	 * The covariance of the zero-=mean gaussian white noise that is added to
	 * the gyro bias error at each time step
	 */
	Vector3f accel_stability_noise;

	/**
	 * The covariance of the GPS rcv clock bias noise, in m^2/s.   Default to 1 us/sqrt(s)
	 */
	float clock_stability_noise;

	/**
	 * The magnitude of the local acceleration due to gravity.
	 */
	float accel_gravity_norm;

	/**
	 * A term for the basic state of the system
	 */
	struct state
	{
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW

		/**
		 * Default-construct an undefined state object.
		 * @return
		 */
		state(){}

		/**
		 * Provided a kalman update vector, apply the vector as an offset to this
		 * state.
		 * @param update A 12-vector to be applied
		 * @return The rotation applied to the mean orientation
		 */
		Quaterniond apply_kalman_vec_update(const Matrix<float, 12, 1>& update);

		/**
		 * Provide a kalman update vector, applied only to the position-time solution.
		 */
		void apply_kalman_vec_update(const Matrix<float, 4, 1>& update);

		/**
		 * An estimate of the bias error in the rate gyros, in radians/second
		 */
		Vector3f gyro_bias;

		/**
		 * An estimate of the orientation of the vehicle.  This quaternion represents
		 * a transformation from ECEF coordinates to the vehicle body frame.
		 */
		Quaterniond orientation;

		/// Velocity in Earth-centered Earth-fixed reference frame, in m/s
		Vector3d velocity;

		/// Estimated bias of the accelerometer, in m/s/s
		Vector3f accel_bias;

		/// Position in Earth-centered Earth-fixed reference frame, in meters
		Vector3d position;

		/// Clock bias term, in meters of light distance traveled.
		float clock_bias;

		/// The acceleration estimate of the vehicle, in the ECEF frame
		Vector3f inertial_accel;

		// The angular rate estimate of the vehicle, in the body frame
		Vector3f body_rate;

		/**
		 * @return True if the state vector contains any NaNs
		 */
		bool has_nan(void) const;

		/**
		 * @return True if the state vector contains any Infs
		 */
		bool has_inf(void) const;

		/**
		 * @return True if the state vector does not contain any NaNs or Infs
		 */
		bool is_real(void) const;

		/**
		 * Print a representation of this object to the stream str.
		 * @param str An output stream.
		 */
		void print(std::ostream& str);
	};

	/// The average state of the filter at any time t.
	state avg_state;

	/// Covariance term.  Elements are ordered exactly as in struct state
	Matrix<float, 12, 12> cov;
	/// Covariance of position and time
	Matrix<float, 4, 4> pt_cov;

	/**
	 * Initialize a new GPS/INS integrator
	 */
	pseudorange_ins_qkf();

	/**
	 * Initialize the attitude to @paramref attitude, and error matrix to
	 * attitude_error.  All cross-covariance terms between this value and others
	 * in the system covariance matrix are set to zero.
	 */
	void init_attitude(const Quaterniond& attitude, const Eigen::Matrix3f& attitude_error);
	/// Initialize the vehicle velocity to vel, with diagonal error vel_error
	void init_velocity(const Vector3d& vel, const Vector3f& vel_error);
	/// Initialize the vehicle position to vel, with diagonal error pos_error
	void init_position(const Vector3d& pos, const Vector3f& pos_error);


	/**
	 * Report an INS observation, to propagate the filter forward by one time
	 * step. The coordinate system is maintained in ECEF coordinates.
	 * TODO: Provide a workspace parameter for storage that should be
	 * carried forward to an observation function carried out in this
	 * time step.
	 * @param gyro_meas The measured angular velocity, in rad/sec
	 * @param accel_meas The measured inertial reference frame acceleration, in m/s
	 * @param dt The elapsed time since the last measurement, in seconds.
	 */
	void predict_ecef(const Vector3f& gyro_meas, const Vector3f& accel_meas, float dt);

	/**
	 * Make a single vector observation, with some angular uncertainty.
	 * Warning: only one of obs_vector or obs_gps_pv_report can be called
	 * after a single call to predict().
	 *	@param ref Reference vector, when the orientation is the identity.
	 *		Must be a unit vector.
	 *	@param obs Vector observation, should not be a unit vector
	 *	@param error one-sigma squared magnitude error in the observation
	 */
	void obs_vector(const Vector3f& ref, const Vector3f& obs, float error);

	/**
	 * Make a single satellite observation
	 * @param sat_pos The position of the satellite, based on ephemeris data.
	 * @param pseudorange The time difference between the satellites transmission
	 * and the reception of the signal, in meters
	 * @param error The estimated one-sigma-squared error in the measurement, exclusive
	 * of the clock error, in m^2.
	 */
	void obs_gps_pseudorange(const Vector3d& sat_pos, double pseudorange, float error);

	/**
	 * Make a single satellite speed observation
	 * @param sat_velocity The velocity vector of the satellite, in ECEF coords,
	 * 	based on ephemeris data
	 * @param deltarange The apparent relative velocity of the satellite, based on
	 * 	GPS doppler signals
	 * @param error The nominal error in the deltarange measurement.
	 */
	void obs_gps_deltarange(const Vector3d& sat_velocity, double deltarange, float error);

	/**
	 * Measure the total angular error between the filter's attitude estimate
	 * and some other orientation.
	 * @param orientation The attitude to compare against
	 * @return The angular difference between them, in radians
	 */
	float angular_error(const Quaterniond& orientation) const;

	/**
	 * Measure the total gyro bias error between the filter's estimate and
	 * some other bias
	 * @param gyro_bias The gyro bias vector to compare against
	 * @return The vector difference between them, in radians/second
	 */
	float gyro_bias_error(const Vector3f& gyro_bias) const;

	/**
	 * Measure the total accel bias error between the filter's estimate and
	 * some other bias
	 * @param accel_bias The accelerometer bias vector to compare against
	 * @return The vector difference between them, in m/s/s
	 */
	float accel_bias_error(const Vector3f& accel_bias) const;

	/**
	 * Determine the statistical distance between an example sample point
	 * and the distribution computed by the estimator.
	 */
	float mahalanobis_distance(const state& sample) const;

private:

	/**
	 * The type of an error term between two state vectors.
	 */
	typedef Eigen::Matrix<float, 16, 1> state_error_t;
	/// Compute the error difference between a sigma point and the mean as: point - mean
	state_error_t sigma_point_difference(const state& mean, const state& point) const;

	/**
	 * A predicate that returns true when all of the class invariants are
	 * met
	 */
	bool invariants_met() const;

	void clear_covariance_block(size_t row, const Eigen::Matrix3f& repl);

public:

	/**
	 * Verify that the covariance and average state are niether NaN nor Inf
	 * @return True, iff no element in the covariance or mean are NaN or Inf
	 */
	bool is_real(void) const;
};

#endif
