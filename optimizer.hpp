#ifndef _OPTIMIZER_H_
	#define _OPTIMIZER_H

	#include <iostream>
	#include <Eigen/Dense>
	#include "alglib/stdafx.h"
	#include "alglib/optimization.h"
	#include "SE3.hpp"
	#include "mytime.hpp"

	using namespace std;
	using namespace Eigen;
	using namespace alglib;

	typedef Matrix<double, 4, 8> Matrix48d;
	typedef Matrix<double, 8, 8> Matrix8d;
	typedef Matrix<double, 6, 1> Vector6d;

	VectorXd run_opt(MatrixXd Data, VectorXd seed, int maxits, double epsx, double epsg, double epsf);

	void cost_function_grad(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr);

	/**
	* General structure for functors
	*/
	struct FunctorOpti
	{
		virtual VectorXd cost(VectorXd& x) = 0;
		virtual MatrixXd grad(VectorXd& x) = 0;
		virtual int setData(MatrixXd Data) = 0;
		virtual int inputs() const = 0;
		virtual int values() const = 0;	
	};

#endif
