#include "quaternions.hpp"

#if SPECIALIZE_ON_DOUBLE
template<>
Eigen::Matrix<float, 3, 1>
log<float>(const Quaternion<float>& q)
{
	std::cout << "going through specialization\n";
	Eigen::AngleAxis<double> res(q.cast<double>());
	return (res.axis() * res.angle()).cast<float>();
}
#endif
