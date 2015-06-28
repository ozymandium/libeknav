#ifndef LIBEKNAV_QUEST_HPP
#define LIBEKNAV_QUEST_HPP

/*
 * quest.hpp
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

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Cholesky>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <Eigen/LU>
#include <Eigen/StdVector>
#include <cstdlib>

using Eigen::Quaternion;
using Eigen::Matrix;


/**
 * Optimal two-point Wahaba's problem solver using optimally weighted TRIAD.
 * @param covariance[out] The 3x3 covariance matrix of the result's confidence
 * @param obs0 The first observation vector.
 * @param ref0 The reference vector corresponding to the first observation vector.
 * @param weight0 The scalar weight corresponding to the first observation vector.
 * @param obs1 The second observation vector
 * @param ref1 The reference vector corresponding to the second observation vector.
 * @param weight1 The scalar weight corresponding to the second observation vector.
 *
 * @return the quaternion that optimally rotates the pair of reference vectors to the
 * pair of observation vectors.
 */
template <typename FloatT>
Quaternion<FloatT>
quaternion_soln_optimal2(Matrix<FloatT, 3, 3>& covariance,
	const Matrix<FloatT, 3, 1>& obs0,
	const Matrix<FloatT, 3, 1>& ref0,
	FloatT weight0,
	const Matrix<FloatT, 3, 1>& obs1,
	const Matrix<FloatT, 3, 1>& ref1,
	FloatT weight1)
{
	// Assume each weight is equal to 1/sigma^2 for that parameters error value.
	// Therefore, to scale the weights on 1-index basis...
	FloatT net_weight = weight0 + weight1;
	FloatT var0 = 1.0/weight0, var1 = 1.0/weight1, var_net = 1.0/net_weight,
		   var_net_squared = 1.0/(net_weight*net_weight);

	// solve Wahaba's problem by finding the attitude matrix that optimally
	// rotates from ref to obs
	FloatT lambda = std::sqrt(var_net_squared*(weight0*weight0 + weight1*weight1 +
		2*weight0*weight1*(obs0.dot(obs1)*ref0.dot(ref1)
			+ obs0.cross(obs1).norm() * ref0.cross(ref1).norm())));

	// Construct third pseudo observation as the cross products of the other
	// two, similar to TRIAD.
	Matrix<FloatT, 3, 1> obs2 = obs0.cross(obs1).normalized();
	Matrix<FloatT, 3, 1> ref2 = ref0.cross(ref1).normalized();

	Matrix<FloatT, 3, 3> attitude_matrix =
		var_net*weight0 / lambda * (
			obs0*ref0.transpose() +
			obs0.cross(obs2) * ref0.cross(ref2).transpose())
		+ var_net*weight1 / lambda * (
			obs1 * ref1.transpose() +
			obs1.cross(obs2) * ref1.cross(ref2).transpose())
		+ obs2 * ref2.transpose();

	// Assume for the moment that each weight has been chosen to be 1/sigma^2
	covariance = (ref0 * ref1.transpose() 
			+ ref1 * ref0.transpose()) * ref0.dot(ref1)*var_net;
	covariance += (var0 - var_net)*ref1*ref1.transpose();
	covariance += (var1 - var_net)*ref0*ref0.transpose();
	covariance *= 1.0/ref0.cross(ref1).squaredNorm();
	covariance += var_net * Matrix<FloatT, 3, 3>::Identity();
	return Quaternion<FloatT>(attitude_matrix);
}

/**
 * Compute the optimal quaternion that rotates from the reference set refs to
 * the observation set obs.  Each observation may be weighted by premultiplying
 * the observation by the scalar weighting constant.  The sum of the weights
 * need not equal unity.  This method uses Markley's SVD.
 */
template <typename FloatT>
Quaternion<FloatT>
quaternion_soln_markley_svd(
	const std::vector<Matrix<FloatT, 3, 1> >& refs,
	const std::vector<Matrix<FloatT, 3, 1> >& obs)
{
	using namespace Eigen;

	// Compute the attitude profile matrix B
	Matrix<FloatT, 3, 3> sum;
	sum.setZero();
	for (size_t i = 0; i != refs.size(); ++i) {
		sum += obs[i]*refs[i].transpose();
	}

	SVD<Matrix<FloatT, 3, 3> > svd(sum);

	Matrix<FloatT, 3, 3> result = svd.matrixU()
		* (Matrix<FloatT, 3, 1>() << 1, 1, svd.matrixU().determinant()*svd.matrixV().determinant()).finished().asDiagonal()
		* svd.matrixV().transpose();
	return Quaternion<FloatT>(result);
}

/**
 * Compute the optimal quaternion that rotates from the reference set refs to
 * the observation set obs.  Each observation may be weighted by premultiplying
 * the observation by the scalar weighting constant.  The sum of the weights
 * need not equal unity.  This method uses Shuster's QUEST.
 */
template <typename FloatT>
Quaternion<FloatT>
quaternion_soln_quest(
	const std::vector<Matrix<FloatT, 3, 1> >& refs,
	const std::vector<Matrix<FloatT, 3, 1> >& obs)
{
	using namespace Eigen;

	// Compute the attitude profile matrix B
	Matrix<FloatT, 3, 3> sum;
	sum.setZero();
	for (size_t i = 0; i != refs.size(); ++i) {
		sum += obs[i]*refs[i].transpose();
	}

	FloatT sigma = sum.trace();
	Matrix<FloatT, 3, 3> S = sum + sum.transpose();
	// The antisymmetic matrix composed of the column vector Z
	Matrix<FloatT, 3, 3> bar_Z = sum - sum.transpose();
	Matrix<FloatT, 3, 1> Z; Z << bar_Z(2, 1), bar_Z(0, 2), bar_Z(1, 0);

	FloatT delta = S.determinant();
	FloatT kappa = S.adjoint().trace();
	FloatT a = sigma*sigma - kappa;
	FloatT b = sigma*sigma + Z.dot(Z);
	FloatT c = delta + Z.dot(S * Z);
	FloatT d = Z.dot(S*(S*Z));
	
	// Perform Newton iteration on the characteristic polynomial
	// lambda^4 - (a+b)lambda^2 - c*lambda + (ab + c*gamma - d) = 0
	FloatT lambda_max = 1.0;

	// Find the gibbs vector associated with this eigenvalue
	Matrix<FloatT, 3, 1> Y;
   	((lambda_max + sigma)*Matrix<FloatT, 3, 3>::Identity() - S).lu().solve(Z, &Y);

	// Convert to a quaternion
	Quaternion<FloatT> ret;
	ret.vec() = Y;
	ret.w() = 1.0;
	ret.coeffs() *= 1.0 / std::sqrt(1 + Y.squaredNorm());
	return ret;
}

/**
 * Optimal two-point Wahaba's problem solver, using QUEST.
 * @param covariance[out] The 3x3 covariance matrix of the result's confidence
 * @param obs0 The first observation vector.
 * @param ref0 The reference vector corresponding to the first observation vector.
 * @param weight0 The scalar weight corresponding to the first observation vector.
 * @param obs1 The second observation vector
 * @param ref1 The reference vector corresponding to the second observation vector.
 * @param weight1 The scalar weight corresponding to the second observation vector.
 *
 * @return the quaternion that optimally rotates the pair of reference vectors to the
 * pair of observation vectors.
 */
template <typename FloatT>
Quaternion<FloatT>
quaternion_soln_quest2(Matrix<FloatT, 3, 3>& covariance,
	const Matrix<FloatT, 3, 1>& obs0,
	const Matrix<FloatT, 3, 1>& ref0,
	FloatT weight0,
	const Matrix<FloatT, 3, 1>& obs1,
	const Matrix<FloatT, 3, 1>& ref1,
	FloatT weight1)
{
	using namespace Eigen;

	// Compute the attitude profile matrix B from the two observations
	Matrix<FloatT, 3, 3> sum = obs0 * ref0.transpose() * weight0;
	sum += obs1 * ref1.transpose() * weight1;

	FloatT sigma = sum.trace();
	Matrix<FloatT, 3, 3> S = sum + sum.transpose();
	// The antisymmetic matrix composed of the column vector Z
	Matrix<FloatT, 3, 3> bar_Z = sum - sum.transpose();
	Matrix<FloatT, 3, 1> Z; Z << bar_Z(2, 1), bar_Z(0, 2), bar_Z(1, 0);

	// Perform Newton iteration on the characteristic polynomial
	// lambda^4 - (a+b)lambda^2 - c*lambda + (ab + c*gamma - d) = 0
	FloatT lambda_max = std::sqrt(weight0*weight0 + weight1*weight1
			+ 2.0 * weight0*weight1 * (ref0.dot(ref1)*obs0.dot(obs1) 
				+ ref0.cross(ref1).norm()*obs0.cross(obs1).norm()));

	// Find the gibbs vector associated with this eigenvalue
	Matrix<FloatT, 3, 1> Y;
   	((lambda_max + sigma)*Matrix<FloatT, 3, 3>::Identity() - S).lu().solve(Z, &Y);

	// Convert to a quaternion
	Quaternion<FloatT> ret;
	ret.vec() = Y;
	ret.w() = 1.0;
	ret.coeffs() *= 1.0 / std::sqrt(1 + Y.squaredNorm());

	// Assume for the moment that each weight has been chosen to be 1/sigma^2
	FloatT net_weight = weight0 + weight1;
	FloatT var0 = 1.0/weight0, var1 = 1.0/weight1, var_net = 1.0/net_weight;
	covariance = (obs0 * obs1.transpose() 
			+ obs0 * obs1.transpose()) * obs0.dot(obs1)*var_net;
	covariance += (var0 - var_net)*obs1*obs1.transpose();
	covariance += (var1 - var_net)*obs0*obs0.transpose();
	covariance *= 1.0/obs0.cross(obs1).squaredNorm();
	covariance += var_net * Matrix<FloatT, 3, 3>::Identity();
	covariance *= 0.25;

	return ret;
}
#endif // !defined (AHRS_QUEST_HPP)

