#ifndef AHRS_SIGMA_POINTS_HPP
#define AHRS_SIGMA_POINTS_HPP

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
	// Note: SVD may be applicable here
	Quaternion<FloatT> ret;
#if 1
	SelfAdjointEigenSolver<Matrix<FloatT, 4, 4> > soln(cols * cols.transpose());
	// std::cout << "Eigenvectors: " << soln.eigenvectors() << "\n";
	// std::cout << "Eigenvalues: " << soln.eigenvalues().transpose() << "\n";
	ret.coeffs() = soln.eigenvectors().col(3);
#else
	SVD<Matrix<FloatT, 4, 4> > svd(cols * cols.transpose());
	// std::cout << "SVD.U: " << svd.matrixU() << std::endl;
	// std::cout << "SVD.sigma: " << svd.singularValues().transpose() << std::endl;
	ret.coeffs() = svd.matrixU().col(0);
#endif
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
	for (auto i = points.begin(), end = points.end(); i != end; ++i) {
		sum += i->toRotationMatrix();
	}

	SVD<Matrix<FloatT, 3, 3> > svd(sum);

	Matrix<FloatT, 3, 3> result = svd.matrixU()
		* (Matrix<FloatT, 3, 1>() << 1, 1, svd.matrixU().determinant()*svd.matrixV().determinant()).finished().asDiagonal()
		* svd.matrixV().transpose();
	return Quaternion<FloatT>(result);
}

// TODO: quaternion_avg_quest, quaternion_avg_foam

#endif
