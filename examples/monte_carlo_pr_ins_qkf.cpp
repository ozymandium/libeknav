/*
 * monte_carlo_ins_qkf.cpp
 *
 *  Created on: Feb 18, 2010
 *      Author: Jonathan Brandmeyer
 */

#include <eknav/pr_ins_qkf.hpp>
#include <eknav/random_vector.hpp>
#include <eknav/posix/random_seed.hpp>

#include <boost/random.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/bind.hpp>
#include <boost/thread.hpp>

#include <iostream>
#include <fstream>
#include <cstdlib>

#include <getopt.h>

struct program_arguments
{
	unsigned n_iterations;
	unsigned n_tests;
	const char* output_dir;
	bool random_seed;
	uint32_t forced_seed;
	const char* seed_list;

	std::vector<uint32_t> seeds;

	program_arguments(int argc, char **argv)
		: n_iterations(10000), n_tests(10), output_dir("./"), random_seed(false)
		, forced_seed(0), seed_list(NULL)
	{
		int opt = 0;
		opterr = 0;
		while ((opt = getopt(argc, argv, "-rf:o:i:n:")) != -1) {
			switch (opt) {
			case 'r':
				random_seed = true;
				break;

			case 'f':
				random_seed = false;
				n_tests = 1;
				forced_seed = boost::lexical_cast<uint32_t>(optarg);
				break;

			case 'o':
				output_dir = optarg;
				break;

			case 'i':
				n_iterations = boost::lexical_cast<unsigned>(optarg);
				break;

			case 'n':
				n_tests = boost::lexical_cast<unsigned>(optarg);
				break;

			case 'l':
				seed_list = optarg;
				random_seed = false;
				forced_seed = 0;
				break;

			case '?':
				std::cerr << "Usage: monte_carlo_test [options]" << "\n"
					<< "\t-r Use truely random seeds\n"
					<< "\t-o output_dir Send output to directory 'output_dir'\n"
					<< "\t-i n_iterations Run each test for n_iterations\n"
					<< "\t-nfor  n_iterations Run n_iterations of the test\n"
					<< "\t-f initial_seed Initialize the PRNG to $initial_seed\n";
				std::exit(EXIT_FAILURE);
			}
		}
	}

	/**
	 * Initialize the random seeds in accordance with the selected program
	 * options.
	 */
	void initialize_seeds(void)
	{
		if (random_seed) {
			// Initialize the seeds from /dev/urandom
			uint32_t test_count = n_tests;
			while (test_count --> 0) {
				seeds.push_back(::random_seed());
			}
		}
		else if (seed_list) {
			// Take the random seeds from a file.  Also take the number of tests
			// to run by the length of this file, too.
			std::ifstream file(seed_list);
			std::string line;
			while (std::getline(file, line)) {
				std::istringstream filter(line);
				uint32_t seed;
				filter >> seed;
				if (filter) {
					seeds.push_back(seed);
				}
			}
		}
		else {
			// Generate them from a consistent PRNG
			boost::mt19937 top_sequence;
			uint32_t test_count = n_tests;
			while (test_count --> 0) {
				seeds.push_back(top_sequence());
			}
		}
	}
};

void print_orientation_decomp(const pseudorange_ins_qkf& model)
{
	Eigen::SelfAdjointEigenSolver<Matrix<float, 12, 12> > soln(
		model.cov);
	std::cout << "\tcovariance matrix:\n" << model.cov;
	std::cout << "\n\teigenvalues:\n" << soln.eigenvalues().transpose();
	std::cout << "\n\teigenvectors:\n" << soln.eigenvectors();
}

struct test_run
{
	double dt;
	bool warned;
	boost::mt19937 entropy;

	// Draw from uniform random number between -pi and pi rad/sec, in all
	// three dimensions
	Vector3d dc_angular_velocity;
	// Draw from uniform random number between -pi and pi rad/sec, in
	// all three dimentions
	Vector3d ac_angular_velocity_amplitude;
	// Draw from log-scaled frequency distribution between 1/20th and 10 Hz
	double ac_angular_velocity_frequency;
	// Initial phase angle, between -pi and pi
	double ac_angular_velocity_phase;

	// Draw from the uniform unit sphere, in four dimensions.
	Quaterniond initial_orientation;

	// Draw from a normal random distribution, with a one-sigma distance of
	// 0.2 radians
	Vector3d gyro_bias;

	Vector3d accel_bias;

	// Draw from the uniform unit sphere, in three dimentions, of diameter
	// 6.378e6 meters
	Vector3d initial_position;
	// Draw from uniform random number, between 0 and 30 m/s
	Vector3d velocity;
	// Uniform reandom static acceleration (aside from gravity), from 0-3G
	Vector3d static_acceleration;

	uint32_t random_seed;

	// The first reference vector's true direction
	Vector3d reference_vector_0;
	// The second reference vector's true direction
	Vector3d reference_vector_1;

	RandomVector<double, 3> gyro_noise;
	RandomVector<double, 3> vector_noise;
	RandomVector<double, 3> gps_position_noise;
	RandomVector<double, 3> gps_dv_noise;
	RandomVector<double, 3> accelerometer_noise;

	test_run(uint32_t random_seed);
	/**
	 * Execute the test run
	 */
	void perform_test(int n_iterations, int imu_interval, int angle_interval, int gps_interval, std::ostream& output);

private:
	struct _: pseudorange_ins_qkf::state
	{
		int time;
	} true_state;

	void print_initialization_params(std::ostream& output);
	void initialize_state(void);
	void increment_state(void);

	void sample_imu(pseudorange_ins_qkf& tester, int interval);
	void sample_angle_0(pseudorange_ins_qkf& tester);
	double sample_gps(pseudorange_ins_qkf& tester);
	void print_state(int time, pseudorange_ins_qkf& model, std::ostream& output);
};

// ADIS 16354, calibrated
// double gyro_bias = 0.25 * M_PI / 180.0;
// ADXRS 610, uncalibrated
double gyro_bias = 2.0 * M_PI / 180.0;
// ADXRS 610, partially calibrated.
// double gyro_bias = 1.0 * M_PI / 180.0;
// ADIS 16365, default, 50 Hz bandwidth, 256 SPS
double gyro_noise = 0.4 * M_PI / 180.0;
// double gyro_bias = 0.025 * M_PI / 180.0;
double gyro_stability = 0.008 * M_PI / 180.0;

double accel_bias = 0.1;
double accel_stability = 1e-4 * accel_bias;

// WAAS-enabled GPS error
double gps_err = 2;
double mag_err = 1.0 / 180.0 * M_PI;
double gps_dv_err = 0.1;

test_run::test_run(uint32_t random_seed)
	: dt(0.001), warned(false), entropy(random_seed), random_seed(random_seed)
	, gyro_noise(Vector3d::Zero(), (Vector3d::Ones()*::gyro_noise*::gyro_noise).asDiagonal(), &entropy)
	, vector_noise(Vector3d::Zero(), (Vector3d::Ones()*std::pow(mag_err, 2)).asDiagonal(), &entropy)
	, gps_position_noise(Vector3d::Zero(), (Vector3d::Ones()*gps_err*gps_err).asDiagonal(), &entropy)
	, gps_dv_noise(Vector3d::Zero(), (Vector3d::Ones()*gps_dv_err*gps_dv_err).asDiagonal(), &entropy)
	, accelerometer_noise(Vector3d::Zero(), (Vector3d::Ones()*0.04*0.04).asDiagonal(), &entropy)

{
	std::cout << random_seed << std::endl;
	boost::uniform_01<boost::mt19937, double> uniform_entropy(entropy);
	using namespace boost;
	using Eigen::VectorXd;

	uniform_real<double> pi_to_pi(-M_PI, M_PI);
	dc_angular_velocity << pi_to_pi(entropy), pi_to_pi(entropy), pi_to_pi(entropy);
	ac_angular_velocity_amplitude << pi_to_pi(entropy), pi_to_pi(entropy), pi_to_pi(entropy);
	ac_angular_velocity_amplitude *= 0.1;
	dc_angular_velocity *= 0.25;

	uniform_real<double> freq(log(1.0/20), log(10));
	ac_angular_velocity_frequency = exp(freq(entropy));
	ac_angular_velocity_phase = pi_to_pi(entropy);

	std::vector<double> orientation = uniform_on_sphere<double>(3)(uniform_entropy);
	initial_orientation = Quaterniond(Eigen::AngleAxis<double>( pi_to_pi(entropy), VectorXd::Map(&orientation[0], orientation.size())));

	uniform_on_sphere<double> d3(3);
	std::vector<double> d3_tmp = d3(uniform_entropy);
	initial_position = VectorXd::Map(&d3_tmp[0], d3_tmp.size()) * 6.378e6;
	d3_tmp = d3(uniform_entropy);
	reference_vector_0 = VectorXd::Map(&d3_tmp[0], d3_tmp.size()).normalized();
	// TODO: Do we need to specify that the reference vectors are farther apart
	// than some minimum threshold?
	d3_tmp = d3(uniform_entropy);
	reference_vector_1 = VectorXd::Map(&d3_tmp[0], d3_tmp.size()).normalized();

	// TODO: Does sigma need to be 0.2, or 0.2*0.2?
	// TODO: Switch to bounded bias calc.
	uniform_real<double> bias(-::gyro_bias, ::gyro_bias);
	gyro_bias << bias(entropy), bias(entropy), bias(entropy);

	uniform_real<double> bias2(-::accel_bias, ::accel_bias);
	accel_bias << bias2(entropy), bias2(entropy), bias2(entropy);

	normal_distribution<double> vel(0, 30);
	velocity << vel(uniform_entropy), vel(uniform_entropy), vel(uniform_entropy);

	double static_accel_mag = uniform_entropy()*0.5*9.81;
	Vector3d accel_direction;
	do {
		d3_tmp = d3(uniform_entropy);
		accel_direction = VectorXd::Map(&d3_tmp[0], 3);
		// Pick an acceleration vector that isn't directly inline with the reference vector.
	} while (fabs(accel_direction.dot(reference_vector_0)) > 0.9848);

	static_acceleration = static_accel_mag * accel_direction;
}

void
test_run::initialize_state(void)
{
	true_state.orientation = initial_orientation.normalized();
	true_state.accel_bias = accel_bias.cast<float>();
	true_state.gyro_bias = gyro_bias.cast<float>();
	true_state.position = initial_position;
	true_state.velocity = velocity;
	true_state.body_rate = (dc_angular_velocity
			+ ac_angular_velocity_amplitude * std::sin(ac_angular_velocity_phase)).cast<float>();
	true_state.time = 0;
}

Quaterniond true_attitude;
Vector3d	true_accel;
Vector3d	true_pos;
Vector3d	true_vel;

void
test_run::increment_state(void)
{
	true_state.time++;
	true_state.position += true_state.velocity * dt + static_acceleration * dt * dt * 0.5;
	true_state.velocity += static_acceleration * dt;

	// gyro_bias remains constant. TODO: Implement a guass-markov
	// chain for the gyro_bias progression.
	true_state.body_rate = (dc_angular_velocity + ac_angular_velocity_amplitude
			* std::sin(true_state.time * dt * ac_angular_velocity_frequency*(2*M_PI) + ac_angular_velocity_phase)).cast<float>();
	true_state.orientation = exp<double>(true_state.body_rate.cast<double>() * dt) * true_state.orientation;
	true_state.orientation = true_state.orientation.normalized();
	true_attitude = true_state.orientation;
	true_accel = static_acceleration;
	true_pos = true_state.position.normalized()*9.81;
	true_vel = true_state.velocity;
}

void
test_run::sample_imu(pseudorange_ins_qkf& tester, int interval)
{
	Vector3d gyro = gyro_noise() + (gyro_bias + true_state.body_rate.cast<double>());
	Vector3d accel = accelerometer_noise() + true_state.accel_bias.cast<double>() + true_state.orientation *
			(static_acceleration + true_state.position.normalized()*9.81);
	// std::cout << "sample_imu, gyro: " << gyro.transpose() << " accel: " << accel.transpose() << "\n";
	tester.predict_ecef(gyro.cast<float>(), accel.cast<float>(), interval * dt);
}

void
test_run::sample_angle_0(pseudorange_ins_qkf& tester)
{
	Vector3d vec0 = vector_noise() + true_state.orientation * reference_vector_0;
	// std::cout << "sample_angle: " << vec.transpose() << "\n";
	// Vector3d vec1 = vector_noise() + true_state.orientation * reference_vector_1;
	tester.obs_vector(reference_vector_0.cast<float>(),
			vec0.normalized().cast<float>(),
			std::pow(mag_err, 2));
	// tester.obs_vector(reference_vector_1,
	//      vec1.normalized(),
	//		std::pow(5 / 180.0 * M_PI, 2));
}

double
test_run::sample_gps(pseudorange_ins_qkf& tester)
{
	Vector3d pos = gps_position_noise() + true_state.position;
	Vector3d vel = gps_dv_noise() + true_state.velocity;
	// std::cout << "sample_gps, pos: " << pos.transpose() << ", vel: " << vel.transpose() << "\n";
	tester.obs_gps_pv_report(pos, vel,
			Vector3f::Ones()*gps_err*gps_err,
			Vector3f::Ones()*gps_dv_err*gps_dv_err);
	return 1.0;
}

void
test_run::perform_test(int n_iterations,
        int imu_interval,
        int angle_interval,
        int gps_interval,
        std::ostream& output)
{
	pseudorange_ins_qkf tester;
	tester.init_position(initial_position + gps_position_noise(), ::gps_err*3 * Vector3f::Ones());
	tester.init_velocity(velocity + gps_dv_noise(), ::gps_dv_err*3 * Vector3f::Ones());
	tester.gyro_white_noise = Vector3f::Ones()*::gyro_noise*::gyro_noise;
	tester.accel_white_noise = Vector3f::Ones()*0.04*0.04;
	tester.gyro_stability_noise = Vector3f::Ones()*pow(::gyro_stability, 2);
	tester.accel_stability_noise = Vector3f::Ones()*accel_stability*accel_stability;
	tester.accel_gravity_norm = 9.81f;
	tester.clock_stability_noise = 300*300;

	print_initialization_params(output);
	initialize_state();
	tester.init_attitude(true_state.orientation, std::pow(20 * M_PI/180, 2)*Eigen::Matrix3f::Identity());
	print_state(-1, tester, output);
	pseudorange_ins_qkf alternate(tester);

	alternate.avg_state.orientation =
			Eigen::AngleAxis<double>(M_PI, reference_vector_0) * alternate.avg_state.orientation;
	// print_orientation_decomp(tester); std::cout << std::endl;
	// tester.avg_state.print(std::cout); std::cout << std::endl;

	n_iterations -= (n_iterations % imu_interval-1);
	pseudorange_ins_qkf *chosen = NULL;
	int i;
	for (i = 0; i < n_iterations; ++i) {
		if (chosen)
			print_state(i, *chosen, output);
		else
			print_state(i, tester, output);

		increment_state();

		if (i % imu_interval == 0) {
			sample_imu(tester, imu_interval);
			// sample_imu(alternate, imu_interval);
		}
		if (i % angle_interval == 0) {
			sample_angle_0(tester);
			// sample_angle_0(alternate);
			output << "alt state: "; print_state(i, alternate, output);

			output << "Sampled angular error\n";
		}
		if (i != 0 && (i % gps_interval) == 0) {
#if 0
			// Choose which INS to proceed with based on which one is farther
			// off.
			double main_err = fabs(sample_gps(tester));
			// tester.limit_gyro_bias(::gyro_bias);
			double alt_err = fabs(sample_gps(alternate));
			if (!chosen && i == gps_interval*2) {
				if (main_err > alt_err)
					chosen = &alternate;
				else
					chosen = &tester;
			}
#endif
			// alternate.limit_gyro_bias(::gyro_bias);
			output << "Sampled GPS\n";
		}
	}
	print_state(i, tester, output);

//	print_orientation_decomp(tester); std::cout << std::endl;
}

void
test_run::print_initialization_params(std::ostream& output)
{
	output << "random_seed: " << random_seed
		<< " gyro bias: " << gyro_bias.transpose()
		<< " reference vector: " << reference_vector_0.transpose()
		<< " reference angle: " << std::acos(reference_vector_0.dot(initial_position.normalized()))
		<< " dc av: " << dc_angular_velocity.transpose()
		<< " dc av angle: " << std::acos(dc_angular_velocity.normalized().dot(reference_vector_0))
		<< " ac av: " << ac_angular_velocity_amplitude.transpose()
		<< " ac av freq: " << ac_angular_velocity_frequency
		<< " net accel: " << (static_acceleration + initial_position.normalized()*9.81).transpose()
		<< " angle error: " << initial_orientation.angularDistance(Quaterniond::Identity())
		<< " p0: " << initial_position.normalized().transpose() << std::endl;
}

void
test_run::print_state( int time,
		pseudorange_ins_qkf& model,
		std::ostream& output)
{
	double mahalanobis_distance = model.mahalanobis_distance(true_state);

	double angular_error = model.avg_state.orientation.angularDistance(
			true_state.orientation) * 180.0/M_PI;
	double ref_angular_error = std::acos((model.avg_state.orientation * reference_vector_0)
			.dot(true_state.orientation * reference_vector_0)) * 180.0 / M_PI;
	double bias_error = (model.avg_state.gyro_bias - true_state.gyro_bias).norm() * 180 / M_PI;
	double position_error = (model.avg_state.position - true_state.position).norm();
	double velocity_error = (model.avg_state.velocity - true_state.velocity).norm();
	double bias2_error = (model.avg_state.accel_bias - true_state.accel_bias).norm();

	output << time
		<< ", " << mahalanobis_distance
		<< ", " << angular_error
		<< ", " << ref_angular_error
		<< ", " << bias_error
		<< ", " << position_error
		<< ", " << velocity_error
		<< ", " << bias2_error
		<< "\n";
#if 0
	if (!warned && bias_error > 0.8) {
		warned = true;
		std::cout << random_seed << "; WARNING: high bias error at time: " << time << "\n";
	}
#endif
	// true_state.print(output); output << "\n\test state: ";
	// model.avg_state.print(output); output << "\n";
	// output << "\ttrue bias: " << true_state.gyro_bias.transpose()
	// 	<< ", model bias: " << model.avg_state.gyro_bias.transpose() << "\n";
}

void
run_some(const program_arguments& args)
{
	// Procedure: For each individual test, the random number generator shall
	// be reseeded from the host os. The seed shall be saved, and each
	// of the initial conditions shall be drawn from their appropriate
	// distributions
	boost::mt19937 top_sequence;

	for (unsigned i = 0; i < args.n_tests; ++i) {
		uint32_t seed = 0;
		if (args.random_seed) {
			seed = random_seed();
		}
		else if (args.forced_seed) {
			seed = args.forced_seed;
		}
		else {
			seed = top_sequence();
		}


		std::ofstream output(args.output_dir + boost::lexical_cast<std::string>(seed) + ".csv");
		test_run run(seed);
		run.perform_test(args.n_iterations, 10, 50, 500, output);
		output << std::flush;
		output.close();
	}
}

int
main(int argc, char **argv)
{
	std::ios::sync_with_stdio(false);

	program_arguments args(argc, argv);
	if (args.random_seed) {
		program_arguments half(args);
		half.n_tests /= 2;
		args.n_tests -= half.n_tests;

		boost::thread rest(boost::bind(run_some, boost::ref(half)));
		run_some(args);
		rest.join();
	}
	else {
		run_some(args);
	}
	return 0;
}
