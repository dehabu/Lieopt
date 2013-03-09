#include "SE3.hpp"

Matrix4d se3Exp(Vector6d V)
{
	Matrix4d S, SSq, SCub;
	S << 0, -V(2), V(1), V(3), V(2), 0, -V(0), V(4), -V(1), V(0), 0, V(5), 0, 0, 0, 0;

	double theta = V.head<3>().norm();
	double thetaSq = theta * theta;
	double thetaCub = thetaSq * theta;

	if(theta == 0)
		return Matrix4d::Identity() + S;
	else
	{
		SSq = S * S;
		SCub = SSq * S;

		return Matrix4d::Identity() + S + (((1.0 - cos(theta)) / thetaSq) * SSq) + (((theta - sin(theta)) / thetaCub) * SCub);
	}
}

/**
* It returns an array of 6 4x4 matrices (a tensor in matrix form)
*/
Matrix4d* JExpSE3(Vector6d V, int n)
{
	Matrix4d* derivative = new Matrix4d[6];

	double da1[] = {	0, 0, 0, 0, 
				0, 0, 1, 0,  
				0, -1, 0, 0, 
				0, 0, 0, 0
			};

	double da2[] = {	0, 0, -1, 0,
				0, 0, 0, 0,
				1, 0, 0, 0,
				0, 0, 0, 0
			};

	double da3[] = {	0, 1, 0, 0,
				-1, 0, 0, 0,
				0, 0, 0, 0,
				0, 0, 0, 0
			};

	double da4[] = {	0, 0, 0, 0,
				0, 0, 0, 0,
				0, 0, 0, 0,
				1, 0, 0, 0
			};

	double da5[] = {	0, 0, 0, 0,
				0, 0, 0, 0,
				0, 0, 0, 0,
				0, 1, 0, 0
			};

	double da6[] = {	0, 0, 0, 0,
				0, 0, 0, 0,
				0, 0, 0, 0,
				0, 0, 1, 0
			};

	//basis for the derivative
	Matrix4d DA1(da1);
	Matrix4d DA2(da2);
	Matrix4d DA3(da3);
	Matrix4d DA4(da4);
	Matrix4d DA5(da5);
	Matrix4d DA6(da6);

	Matrix4d A, Apower;
	Matrix4d Je1, Je2, Je3, Je4, Je5, Je6;
	Matrix4d DAcc1, DAcc2, DAcc3, DAcc4, DAcc5, DAcc6;
	double fact = 1;

	// A = hat(V)
	A << 0, -V(2), V(1), V(3), V(2), 0, -V(0), V(4), -V(1), V(0), 0, V(5), 0, 0, 0, 0;

	DAcc1 = DA1;
	DAcc2 = DA2;
	DAcc3 = DA3;
	DAcc4 = DA4;
	DAcc5 = DA5;
	DAcc6 = DA6;
	Apower = A;

	Je1 = DA1;
	Je2 = DA2;
	Je3 = DA3;
	Je4 = DA4;
	Je5 = DA5;
	Je6 = DA6;	

	for(int i=2; i<=n; i++)
	{
		//i-th power (wrt a,b,c)
		DAcc1 = DA1*Apower + A*DAcc1;
		DAcc2 = DA2*Apower + A*DAcc2;
		DAcc3 = DA3*Apower + A*DAcc3;
		DAcc4 = A*DAcc4;
		DAcc5 = A*DAcc5;
		DAcc6 = A*DAcc6;

		fact *= i;

		Je1 = Je1 + (DAcc1 / fact);
		Je2 = Je2 + (DAcc2 / fact);
		Je3 = Je3 + (DAcc3 / fact);
		Je4 = Je4 + (DAcc4 / fact);
		Je5 = Je5 + (DAcc5 / fact);
		Je6 = Je6 + (DAcc6 / fact);

		if(i < n)
			Apower = Apower * A;
	}

	derivative[0] = Je1;
	derivative[1] = Je2;
	derivative[2] = Je3;
	derivative[3] = Je4;
	derivative[4] = Je5;
	derivative[5] = Je6;

	return derivative;
}
