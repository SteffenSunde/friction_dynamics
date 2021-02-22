#ifndef DATA_STRUCTURES_HPP
#define DATA_STRUCTURES_HPP

#define _USE_MATH_DEFINES
#include <math.h>

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#include <Eigen/Dense>


using Vec3 = Eigen::Matrix<long double, 3, 1>;
using Vec = Eigen::VectorXd;

#endif