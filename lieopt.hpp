#ifndef _LIEOPT_H_
	#define _LIEOPT_H_

	#include <iostream>
	#include <Eigen/Dense>
	#include "optimizer.hpp"
	#include "mytime.hpp"

	using namespace std;
	using namespace Eigen;

	typedef Matrix<double, 4, 8> Matrix48d;
	typedef Matrix<double, 8, 8> Matrix8d;
	typedef Matrix<double, 6, 1> Vector6d;

	MatrixXd lieOpt(MatrixXd Xsrc, MatrixXd Xdst, double threshold, int maxits=1000, double epsx=0, double epsg=0, double epsf=0);

	Matrix4d recalculateTrInliersSE3(Matrix4d Tr, MatrixXd Xsrc, MatrixXd Xdst, double threshold);

	Matrix4d horn_fit(MatrixXd Xsrc, MatrixXd Xdst);

	Vector4d disambiguateQuaternion(Vector4d q);

	Matrix4d denormalizeTransformation(const Matrix4d Tr, const Matrix4d Tsrc, const Matrix4d Tdst);

	pair<MatrixXd, Matrix4d> normalizeHomoData(MatrixXd D);

#endif
