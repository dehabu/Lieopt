#ifndef _SE3_HPP_
#define _SE3_HPP_

	#include <iostream>
	#include <Eigen/Dense>

	using namespace Eigen;
	using namespace std;

	typedef Matrix<double, 6, 1> Vector6d;


	Matrix4d se3Exp(Vector6d V);

	Matrix4d* JExpSE3(Vector6d V, int n);

	void getDiffTime(timeval t1,  timeval t2, char* msg);

#endif

