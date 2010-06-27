#ifndef AHRS_INS_QKF_HPP
#define AHRS_INS_QKF_HPP

#include "sigma_points.hpp"
#include "quaternions.hpp"
#include <Eigen/StdVector>

using Eigen::Vector3f;
using Eigen::Vector3d;
using Eigen::Vector2d;
using Eigen::Quaterniond;
using Eigen::Matrix;
using Eigen::aligned_allocator;

/// A satellite vehicle raw observation
struct sv_obs
{
	typedef std::vector<sv_obs, aligned_allocator<sv_obs> > container_t;
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	/// The position of the satellite, in ECEF coords, in meters
	Vector3d position;
	/// The velocity of the satellite, in ECEF coords, in meters/second
	Vector3d velocity;
	/** The doppler-determined relative magnitude of velocity between the antenna
	 * and the satellite, in meters/second
	 */
	double dv;
	/// The RMS error of the satellite dopper, in (m/s)^2
	double dv_error;
	/// The distance to the satellite, in meters
	double range;
	/// The RMS error of the satellite range, in meters^2
	double range_error;
};

struct sv_error
{
	double dv_error;
	double range_error;
};

struct basic_ins_qkf
{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	/// The maximum number of satellites that may be tracked by the filter
	static const size_t max_sv = 12;
	/**
	 * The covariance of the zero-mean gaussian white noise that is added to
	 * the gyro bias at each time step.  This value is treated as a diagonal
	 *  matrix, in units of radians^2/second^2.
	 */
	const Vector3d gyro_stability_noise;

	/**
	 * The covariance of the zero-mean gaussian white noise that is added to
	 * the gyro measurement at each time step, in rad^2/second^2
	 */
	const Vector3d gyro_white_noise;

	/**
	 * The covariance of the zero-mean gaussian white noise that is added to
	 * the accelerometer measurement at each time step, in (m/s/s)^2
	 */
	const Vector3d accel_white_noise;

	/**
	 * A term for the basic state of the system
	 */
	struct state
	{
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		typedef std::vector<state, aligned_allocator<state> > container_t;

		/**
		 * Construct a state variable from mean and error terms
		 * @param mean The average state
		 * @param error An error vector of length 12 or greater
		 */
		template <typename Vector_T>
		state(const state& mean, const Vector_T& error)
			: gyro_bias(mean.gyro_bias + error.template segment<3>(0))
			, orientation(mean.orientation * exp<double>(error.template segment<3>(3)))
			, position(mean.position + error.template segment<3>(6).template cast<double>())
			, velocity(mean.velocity + error.template segment<3>(9))
		{
			assert(!has_nan());
		}

#if 0
		template <typename Vector_T>
		state(const state& mean, const Vector_T& error, bool)
			: gyro_bias(mean.gyro_bias + error.segment(0, 3))
			, orientation(mean.orientation * exp<double>(error.segment(3, 3)))
			, position(mean.position + error.segment(6, 3).template cast<double>())
			, velocity(mean.velocity + error.segment(9, 3))
		{
			assert(!has_nan());
		}
#endif
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
		Quaterniond apply_kalman_vec_update(const Matrix<double, 12, 1> update);
		Quaterniond apply_left_kalman_vec_update(const Matrix<double, 12, 1> update);

		/**
		 * An estimate of the bias error in the rate gyros, in radians/second
		 */
		Vector3d gyro_bias;

		/**
		 * An estimate of the orientation of the vehicle.  This quaternion represents
		 * a transformation from ECEF coordinates to the vehicle body frame.
		 */
		Quaterniond orientation;

		/// Position in Earth-centered Earth-fixed reference frame, in meters
		Vector3d position;

		/// Velocity in Earth-centered Earth-fixed reference frame, in m/s
		Vector3d velocity;

		/**
		 * @return True if the state vector contains any NaNs
		 */
		bool has_nan(void) const;

		/**
		 * @return True if the state vector does not contain any NaNs or Infs
		 */
		bool is_real(void) const;

		/**
		 * Print a representation of this object to the stream str.
		 * @param str An output stream.
		 */
		virtual void print(std::ostream& str);
	};

	/// The average state of the filter at any time t.
	state avg_state;

	/// Covariance term.  Elements are ordered exactly as in struct state
	Matrix<double, 12, 12> cov;

	/**
	 * Initialize a new basic INS QFK
	 * @param pos_estimate Initial estimate of the position
	 * @param pos_error one-sigma initial bounds for position error
	 * @param bias_error one-sigma initial bounds for bias error in the gyros
	 * @param v_error one-sigma bounds for velocity error (initial v == 0)
	 * @param gyro_white_noise The diagonal matrix of gyro white noise
	 * @param gyro_stability_noise The diagonal matrix of gyro instability noise
	 * @param accel_white_noise The diagonal matrix of accelerometer white noise
	 */
	basic_ins_qkf(const Vector3d& pos_estimate,
			double pos_error, double bias_error, double v_error,
			const Vector3d& gyro_white_noise,
			const Vector3d& gyro_stability_noise,
			const Vector3d& accel_white_noise,
			const Vector3d& vel_estimate = Vector3d::Zero());

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
	void predict(const Vector3d& gyro_meas, const Vector3d& accel_meas, double dt);

	/**
	 * Report an INS observation, to propagate the filter forward by one time
	 * step. This function differs from predict() in that it uses the NED frame
	 * instead of the ECEF frame.
	 */
	void predict_ned(const Vector3d& gyro_meas, const Vector3d& accel_meas, double dt);

	/**
	 * Make a single vector observation, with some angular uncertainty.
	 * Warning: only one of obs_vector or obs_gps_pv_report can be called
	 * after a single call to predict().
	 *	@param ref Reference vector, when the orientation is the identity.
	 *		Must be a unit vector.
	 *	@param obs Vector observation, should not be a unit vector
	 *	@param error one-sigma squared magnitude error in the observation
	 */
	void obs_vector(const Vector3d& ref, const Vector3d& obs, double error);

	/**
	 * Incorporate a GPS PVT report.
	 * @param pos The sensor position, in earth-centered, earth-fixed coords, meters
	 * @param vel The sensor velocity, in meters/second
	 * @param p_error The RMS position error, (m)^2
	 * @param v_error The RMS velocity error, (m/s)^2
	 */
	void obs_gps_pv_report(const Vector3d& pos, const Vector3d& vel, const Vector3d& p_error, const Vector3d v_error);

	/**
	 * Incorporate a GPS position report, in either ECEF or NED coordinates.
	 * @param pos The position, in meters
	 * @param p_error The position error, in meters.
	 */
	void obs_gps_p_report(const Vector3d& pos, const Vector3d& p_error);

	/**
	 * Incorporate a GPS velocity report, in ECEF 3d coordinates.
	 * @param vel The 3d velocity, relative to the fixed earth frame, in (m/s).
	 * @param v_error The one-sigma RMS velocity error (m/s)^2
	 */
	void obs_gps_v_report(const Vector3d& vel, const Vector3d& v_error);

	/**
	 * Observe a GPS vector track over ground report, in north-east-down coordinates
	 * @param vel The 2d velocity value, parallel to the ground, in m/s
	 * @param v_error The one-sigma RMS velocity error (m/s)^2
	 */
	void obs_gps_vtg_report(const Vector2d vel, const double v_error);

	/**
	 * Directly observe the gyro sensor bias. In practice, we cannot do this. However,
	 * the true bias is not a random walk. It tends to return towards zero when it is
	 * farther away from zero (not temperature dependant). Therefore, we can
	 * incorporate this extra knowledge through a periodic "observation" of zero 
	 * bias with a large error.
	 *
	 * Use with extreme caution. It is almost certainly better to just clamp the
	 * bias term after making an observation.
	 *
	 * @param bias The observed bias, in radians/sec
	 * @param bias_error The one-sigma estimate of the gyro bias error, in radians/sec
	 */
	void obs_gyro_bias(const Vector3d& bias, const Vector3d& bias_error);

	/**
	 * Incorporate a "raw" satellite observation.  The observation should be
	 * fully filtered for ionospheric, atmospheric, and other error models.
	 *
	 * @param sv_pos list of satellite positions, in ECEF
	 * @param sv_vel list of satellite velocities, in ECEF
	 * @param sv_range list of satellite distances, in meters
	 * @param sv_range_error list of satellite distance error 1-sigma^2 bounds,
	 * 	in meters^2
	 * @param sv_dv list of relative velocities, in meters/second
	 * @param sv_dv_error list of relative velocity error 1-sigma^2 bounds,
	 *  in (m/s)^2
	 */
	void obs_gps_sv(const sv_obs::container_t& sv_pos);

	double angular_error(const Quaterniond& orientation) const;
	double gyro_bias_error(const Vector3d& gyro_bias) const;

	/**
	 * Determine the statistical distance between an example sample point
	 * and the distribution computed by the estimator.
	 */
	double mahalanobis_distance(const state& sample) const;

private:
	struct process_state : state
	{
		typedef std::vector<process_state, aligned_allocator<process_state> > container_t;
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		/**
		 * Construct a process augmented state from mean and error terms
		 * @param mean The t-1 mean.
		 * @param error An error vector of length 21 or greater
		 * @param gyro_mean The observed gyro measurement
		 * @param accel_mean The observed accellerometer measurement
		 */
		template <typename Matrix_T>
		process_state(const state& mean, const Matrix_T& error, const Vector3d& gyro_mean, const Vector3d& accel_mean)
			: state(mean, error)
			, gyro_instability(error.template segment<3>(12))
			, gyro_meas(gyro_mean + error.template segment<3>(15))
			, accel_meas(accel_mean + error.template segment<3>(18))
		{
		}

		/**
		 * Construct a process augmented state vector from mean and error terms,
		 * but do not include any error terms for the augmented state.
		 * @param mean The t-1 mean
		 * @param error An error vector of length 12 or greater
		 * @param gyro_mean The observed gyro measurement
		 * @param accel_mean The observed accellerometer measurement
		 * @param (anon) A dummy parameter used to disambiguate the signature of this function from the
		 * primary constructor
		 */
		template <typename Matrix_T>
		process_state(const state& mean, const Matrix_T& error, const Vector3d& gyro_mean, const Vector3d& accel_mean, bool)
			: state(mean, error)
			, gyro_instability(Vector3d::Zero())
			, gyro_meas(gyro_mean)
			, accel_meas(accel_mean)
		{
		}

		/**
		 * Construct a process augmented state vector that is nearly a copy of
		 * the average state, without any error terms.
		 * @param mean The t-1 mean
		 * @param gyro_mean The observed gyro measurement
		 * @param accel_mean The observed accelerometer measurement
		 * @return
		 */
		process_state(const state& mean, const Vector3d& gyro_mean, const Vector3d& accel_mean)
			: state(mean)
			, gyro_instability(Vector3d::Zero())
			, gyro_meas(gyro_mean)
			, accel_meas(accel_mean)
		{
		}

		/// The amount by which the gyro null point has drifted in this time step.
		Vector3d gyro_instability;
		/// The measurement of the angular rate gyro's, in rad/sec
		Vector3d gyro_meas;
		/// The measurement of the accelerometer, in m/sec/sec
		Vector3d accel_meas;
		void print(std::ostream& str);
	};

	/**
	 * The type of an error term between two state vectors.
	 */
	typedef Eigen::Matrix<double, 12, 1> state_error_t;

	/** Produce a set of 24 equally-weighted sigma points from avg_state and cov
	 * and load them into @paramref points.
	 */
	void decompose_sigma_points(process_state::container_t& points, const Vector3d& gyro_meas, const Vector3d& accel_meas, double dt);
	/** Project one sigma point through the non-linear process equation.
	 * TODO: Should this be a member function of process_state?
	 */
	state project_sigma_point(const process_state& p, double dt);
	/// Compute the median value of a set of sigma points.
	state average_sigma_points(const state::container_t& points);
	/// Compute the error difference between a sigma point and the mean as: point - mean
	state_error_t sigma_point_difference(const state& mean, const state& point) const;

	/**
	 * Augmented state vector that incorporates linear sensor noise from the
	 * vector observer.
	 */
	struct vector_obs_state : state
	{
		typedef std::vector<vector_obs_state, aligned_allocator<vector_obs_state> > container_t;
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		/**
		 * Construct an augmented state for making a vector observation
		 * @param mean The a priori mean
		 * @param error An error vector of length 15 or greater
		 */
		template <typename Matrix_T>
		vector_obs_state(const state& mean, const Matrix_T& error)
			: state(mean, error)
			, meas_error(error.template segment<3>(12))
		{
		}
		/// The measurement error, in the units of the observer.
		Vector3d meas_error;
	};
	/**
	 * Produce a set of sigma points and their error vectors when performing a
	 * 	vector observation.
	 * @param points A container to receive the set of sigma points
	 * @param errors A matrix to receive the set of error vectors (deviation
	 * 	from the sigma points)
	 * @param error The sigma-squared uncertainty in the vector observation,
	 * 	in the units of the sensor.
	 */
	void decompose_sigma_points(vector_obs_state::container_t& points, Matrix<double, 15, 30>& errors, double error);
	/// Compute the minimal quaternion that rotates ref into its expected position
	Quaterniond observe_vector(const vector_obs_state& point, const Vector3d& ref_vector);

	struct gps_raw_state : state
	{
		typedef std::vector<gps_raw_state, aligned_allocator<gps_raw_state> > container_t;
		template<typename VectorT>
		gps_raw_state(const state& mean, const VectorT& error);

		sv_error obs_errors[max_sv];
	};

	struct gps_estimated_obs
	{
		double range[max_sv];
		double dv[max_sv];
	};

	void decompose_sigma_points(gps_raw_state::container_t& points, const sv_obs::container_t& sv_pos);
	gps_estimated_obs observe_sv(const gps_raw_state& point, const sv_obs::container_t& sv_pos);
	gps_estimated_obs average_sigma_points(const std::vector<gps_estimated_obs>& points, int n_sv);
	Matrix<double, Eigen::Dynamic, 1> sigma_point_difference(const gps_estimated_obs& median, const gps_estimated_obs& point, int n_sv);

public:
	/** Perform the posterior counter-rotation of the covariance matrix by
	  the update that gets applied to the estimated state.
	  */
	void counter_rotate_cov(const Quaterniond& update);

	/**
	 * Verify that the covariance and average state are niether NaN nor Inf
	 * @return True, iff no element in the covariance or mean are NaN or Inf
	 */
	bool is_real(void) const;
};

template<typename VectorT>
basic_ins_qkf::gps_raw_state::gps_raw_state(const state& mean, const VectorT& error)
	: state(mean, error)
{
	const unsigned base = 12;
	const unsigned stride = 2;
	for (unsigned sv = 0; sv < max_sv && sv*stride + base < (unsigned)error.rows(); ++sv) {
		obs_errors[sv].dv_error = error(base + sv*stride);
		obs_errors[sv].range_error = error(base + sv*stride + 1);
	}
	// Clear remaining error terms with something loud
	for (unsigned sv = (error.rows() - 12)/2; sv != max_sv; ++sv) {
		obs_errors[sv].dv_error = NAN;
		obs_errors[sv].range_error = NAN;
	}
}

#endif
