#ifndef LIBEKNAV_SIGMA_POINTS_HPP
#define LIBEKNAV_SIGMA_POINTS_HPP

/*
 * sigma_points.hpp
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
#include <vector>

using Eigen::Quaternion;
using Eigen::Matrix;

template <typename FloatT>
void decompose_sigma_points(std::vector<Quaternion<FloatT> >& ret,
		const Quaternion<FloatT>& avg, const Matrix<FloatT, 3, 3>& dev)
{
	using namespace Eigen;

	ret.clear();
	LLT<Matrix<FloatT, 3, 3> > matrix_sqrt(dev);
	ret.push_back(avg);
	for (int i = 0; i < 3; ++i) {
		// TODO: Excess math required due to sqrt in the norm?
		FloatT mag = matrix_sqrt.col(i).norm();
		Matrix<FloatT, 3, 1> axis = matrix_sqrt.col(i) / mag;
		ret.push_back(Quaternion<FloatT>(AngleAxis<FloatT>(mag, axis)) * avg);
		ret.push_back(Quaternion<FloatT>(AngleAxis<FloatT>(-mag, axis)) * avg);
	}
}

// Compute the quaternion average using Johnson's eigenvalue decomposition method
template <typename FloatT>
Quaternion<FloatT>
quaternion_avg_johnson(const std::vector<Quaternion<FloatT> >& points)
{
	using namespace Eigen;
	// Compute an average quaternion using the methods of Johnson's Ph. D.
	// thesis
	Matrix<FloatT, 4, Dynamic> cols(4, points.size());
	for (std::size_t i = 0, end = points.size(); i != end; ++i) {
		cols.col(i) = points[i].coeffs();
	}

	Quaternion<FloatT> ret;
	SelfAdjointEigenSolver<Matrix<FloatT, 4, 4> > soln(cols * cols.transpose());
	ret.coeffs() = soln.eigenvectors().col(3);
	return ret;
}


// Compute the quaternion average using the Markley SVD method
template <typename FloatT>
Quaternion<FloatT>
quaternion_avg_markley(const std::vector<Quaternion<FloatT> >& points)
{
	using namespace Eigen;
	// double scale = 1.0 / points.size();
	Matrix<FloatT, 3, 3> sum;
	sum.setZero();
	for (int i = 0, end = points.size(); i != end; ++i) {
		sum += points[i].toRotationMatrix();
	}

	SVD<Matrix<FloatT, 3, 3> > svd(sum);

	Matrix<FloatT, 3, 3> result = svd.matrixU()
		* (Matrix<FloatT, 3, 1>() << 1, 1, svd.matrixU().determinant()*svd.matrixV().determinant()).finished().asDiagonal()
		* svd.matrixV().transpose();
	return Quaternion<FloatT>(result);
}

// TODO: quaternion_avg_quest, quaternion_avg_foam

#endif
