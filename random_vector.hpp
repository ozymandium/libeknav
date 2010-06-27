#ifndef AHRS_RANDOM_VECTOR_HPP
#define AHRS_RANDOM_VECTOR_HPP

#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/StdVector>

// Note: Using Boost.Random by default, but compatible with TR1::random.
#include <boost/random.hpp>

/**
 * A normal random vector generator
 * @param FloatT Generate Eigen::Vector's of this type of float
 * @param Length The number of dimensions of the vector.  Eigen::Dynamic
 * is allowed
 * @param Prng The underlying continuous random vector type.  Defaults to
 * a linear congruential generator.
 */
template<typename FloatT, int Length, typename Prng = boost::mt19937>
class RandomVector
{
public:
	/// Required to ensure alignment requirements of this type when allocated
	/// on the heap.
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	/// The type of a sample.
	typedef Eigen::Matrix<FloatT, Length, 1> result_t;
	/// The type of the sample covariance matrix
	typedef Eigen::Matrix<FloatT, Length, Length> cov_t;

private:
	/// The distribution mean
	result_t mean;
	/// The matrix square root of the covariance matrix.  cov_sqrt*cov_sqrt.transpose()
	/// will return the original covariance matrix.
	cov_t cov_sqrt;
	/// A standard normal random variable.  Successive draws from this
	/// generator are all independent.
	boost::variate_generator<Prng*, boost::normal_distribution<FloatT> > std_norm;

public:
	/**
	 * Create a normal random vector generator
	 * @param mean The mean vector mu of the distribution
	 * @param cov The covariance matrix P of the distribution
	 * @param engine A borrowed pointer to a common entropy pool.  The PRNG should be a model
	 * 		of UniformRandomNumberGenerator.  See
	 * 		http://www.boost.org/doc/libs/1_40_0/libs/random/random-concepts.html
	 * The lifetime of @paramref *engine must exceed the lifetime of *this.
	 */
	RandomVector(const result_t& mean, const cov_t& cov, Prng* engine);

	/**
	 * Draw a new sample from the random vector.
	 */
	result_t operator()(void);

	const result_t& get_mean() const { return mean; }
};

template<typename FloatT, int Length, typename Prng>
RandomVector<FloatT, Length, Prng>::RandomVector(const result_t& mean, const cov_t& cov, Prng* engine)
	: std_norm(engine, boost::normal_distribution<FloatT>(0, 1))
{
	this->mean = mean;
	// Compute the matrix square root by the Cholesky decomposition
	cov_sqrt = cov.llt().matrixL();
}

template<typename FloatT, int Length, typename Prng>
typename RandomVector<FloatT, Length, Prng>::result_t
RandomVector<FloatT, Length, Prng>::operator()(void)
{
	// TODO: Technically, this is supposed to be a (template) function of
	// the prng, that takes a prng pointer as as an argument.
	result_t norms;
	for (size_t i = 0; i != Length; ++i) {
		norms[i] = std_norm();
	}
	// See http://en.wikipedia.org/wiki/Multivariate_normal_distribution#Drawing_values_from_the_distribution
	// for justification.
	return mean + cov_sqrt * norms;
}

#endif
