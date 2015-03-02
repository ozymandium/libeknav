/*
 * test_ins_qkf.cpp
 *
 *  Created on: Aug 17, 2009
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

#include "ins_qkf.hpp"
#include <iostream>
#include <iomanip>
#include "random_vector.hpp"
#include <boost/lexical_cast.hpp>
#include "timer.hpp"

Quaterniond
initial_orientation(void)
{
	Matrix<double, 3, 3> orientation;
	orientation.col(2) = Vector3d(0.707, 0.707, 0).normalized();
	orientation.col(1) = Vector3d(-0.707, 0.707, 0).normalized();
	orientation.col(0) = orientation.col(1).cross(orientation.col(2));
	return Quaterniond(orientation);
}

void
print_eigenvalues(const basic_ins_qkf& model)
{
	for (int i = 0; i < 12; i += 3) {
		Eigen::SelfAdjointEigenSolver<Matrix<double, 3, 3> > soln(model.cov.block<3, 3>(i, i));
		std::cout << soln.eigenvalues().transpose() << " ";
	}
}

void
print_orientation_decomp(const basic_ins_qkf& model)
{
	Eigen::SelfAdjointEigenSolver<Matrix<double, 12, 12> > soln(
		model.cov);
	std::cout << "eigenvalues: " << soln.eigenvalues().transpose();
	std::cout << "\n\teigenvectors:\n" << soln.eigenvectors();
}

int
main(int argc, char **argv)
{
	if (argc != 2) {
		std::cout << "usage: test_ins_qkf n_iterations\n";
		return 1;
	}
	const unsigned n_iterations = boost::lexical_cast<unsigned>(argv[1]);

	std::cout << std::setw(11) << std::scientific;
	auto entropy = boost::mt19937();
	const double gyro_interval = 0.001;
	Vector3d angular_velocity(0, M_PI * 0.85, 0);
	const double mag_error = std::pow(5 / 180.0 * M_PI, 2);
	const double gps_pos_error = 10*10;
	const double gps_dv_error = 0.1*0.1;
	const double accel_error = 0.04*0.04;
	const double gyro_error = 0.01;
	Vector3d position(1017.67e3, -5079.282e3, 3709.041e3);
	Vector3d gravity = 9.81 * position.normalized();
	RandomVector<double, 3, boost::mt19937> gyros(/*Vector3d::Zero(),*/Vector3d(0.15, 0.2, -0.2), (Vector3d::Ones()*gyro_error).asDiagonal(), &entropy);
	RandomVector<double, 3, boost::mt19937> obs_0(Vector3d::UnitZ(), (Vector3d::Ones()*mag_error).asDiagonal(), &entropy);
	RandomVector<double, 3, boost::mt19937> obs_1(Vector3d::UnitY(), (Vector3d::Ones()*mag_error).asDiagonal(), &entropy);
	RandomVector<double, 3, boost::mt19937> gps_pos(position, (Vector3d::Ones()*gps_pos_error).asDiagonal(), &entropy);
	RandomVector<double, 3, boost::mt19937> gps_dv(Vector3d::Zero(), (Vector3d::Ones()*gps_dv_error).asDiagonal(), &entropy);
	RandomVector<double, 3, boost::mt19937> accel(gravity, (Vector3d::Ones()*accel_error).asDiagonal(), &entropy);

	Quaterniond q_angular_velocity(exp<double>(angular_velocity * gyro_interval));
	Quaterniond d_orientation = Quaterniond::Identity();
	Quaterniond orientation = initial_orientation();
	// basic_ins_qkf tester(Vector3d(0, 0, 0), 4e6, 0.2, 10,
	basic_ins_qkf tester(position, 1e4, 0.447, 10,
			// TODO: Beware that none of these can be < sqrt(eps(double))
			Vector3d::Ones()*gyro_error,
			Vector3d::Ones()*0.00001,
			Vector3d::Ones()*accel_error);
	timer clock;

	std::cout << "Initial covariance:\n" << tester.cov << "\n";
	clock.start();

	basic_ins_qkf::state s;
	s.orientation = d_orientation * orientation;
	s.gyro_bias = gyros.get_mean();
	s.position = position;
	s.velocity = gps_dv.get_mean();
	std::cout << tester.mahalanobis_distance(s) << "\n";
	for (unsigned i = 0; i < n_iterations; ++i) {
		d_orientation = (q_angular_velocity * d_orientation).normalized();
		// TODO: The gravity vector isn't right here... Need to formalize the orientation
		// quaternion
		if (i % 4 == 0) {
			// INS observations at 250 Hz
			tester.predict(gyros() + angular_velocity, d_orientation * orientation * accel(), 0.005);
		}

		if (i % 50 == 0) {
			// Angular observation at 20 Hz
			// std::cout << "observed angle\n";
			// print_orientation_decomp(tester); std::cout << "\n";
			tester.obs_vector(Vector3d::UnitZ(), (d_orientation * orientation * obs_0()).normalized(), mag_error);
			// print_orientation_decomp(tester); std::cout << "\n";
			// tester.obs_vector(Vector3d::UnitY(), (d_orientation * orientation * obs_1()).normalized(), mag_error);
		}
		if (i % 100 == 0) {
			// GPS observations at 4 Hz
			// std::cout << "observed GPS\n";
			tester.obs_gps_pv_report(gps_pos(), gps_dv(), Vector3d::Ones()*gps_pos_error, Vector3d::Ones()*gps_dv_error);
		}
	}
	double obs_us = clock.stop() * 1e6;
	std::cout << "\nFilter time: " << obs_us << "us" << std::endl;
	std::cout << "*** Final State ***\n";
	// std::cout << "covariance eigenvalues: "; print_eigenvalues(tester);
	std::cout << "\norientation cov: "; print_orientation_decomp(tester);
	std::cout << "\norientation cov matrix:\n" << tester.cov;
	std::cout << "\n\nactual orientation: " << (d_orientation * orientation).coeffs().transpose();
	std::cout << "\nest orientation: " << tester.avg_state.orientation.coeffs().transpose();
	std::cout << "\norientation error: " << tester.avg_state.orientation.angularDistance(d_orientation * orientation) * 180.0/M_PI << " deg";
	if (tester.avg_state.orientation.dot(d_orientation * orientation) < 0) {
		tester.avg_state.orientation.coeffs() *= -1;
	}
	std::cout << "\nvector orientation error: " << log<double>((d_orientation * orientation).conjugate()*tester.avg_state.orientation).transpose() * 180.0/M_PI << " deg";
	std::cout << "\noriginal orientation angle: " << log<double>(orientation).transpose() * 180.0/M_PI;
	std::cout << "\norientation difference angle: " << log<double>(d_orientation).transpose()*180.0 / M_PI;
	std::cout << "\nest bias: " << tester.avg_state.gyro_bias.transpose();
	std::cout << "\nest position: " << tester.avg_state.position.transpose();
	std::cout << "\nest velocity: " << tester.avg_state.velocity.transpose() << std::endl;
	return 0;
}
