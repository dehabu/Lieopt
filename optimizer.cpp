#include "optimizer.hpp"

//Timerx timer;

/**
* The specific functor I'm interesting in
*/
struct l2costFunctor : public FunctorOpti
{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW 

	//The Reduced Measurement Matrix
	// It is calculated from data assocciations between images --> constraints
	Matrix8d RM;
	//The point on the Lie group
	Matrix4d expX;
	//M2 is [se3Exp(x), -eye(4)], a 4x8 matrix
	Matrix48d M2;
	
	double Error;

	l2costFunctor() : FunctorOpti() {}

	VectorXd cost(VectorXd& x)
	{
		expX = se3Exp(x);
	
		//build M2 as [se3Exp(x), -eye(4)] 4x8 matrix
		M2.block<4, 4>(0, 0) = expX;
		M2.block<4, 4>(0, 4) = -Matrix4d::Identity();

		Error = (M2*RM*M2.transpose()).trace();

		VectorXd F(1);
		F(0) = Error;

		return F;
	}
 
	MatrixXd grad(VectorXd& x)
	{
		//derivative of the exponential map at x
		// x can be something different from 0! :D
		// at 0 the derivative is rather simple
		Matrix4d* dexp = JExpSE3(x, 20);

		//we assume that Error = trace(M2*RM*M2') has been computed before
		//when calling the operator (). (there is no need for computing it twice)
		MatrixXd Lf = (RM.transpose() + RM)*M2.transpose();
		Matrix4d L = Lf.block(0, 0, 4, 4);
	
		//The Jacobian, a 1x6 vector
		MatrixXd J(1, 6);
		//It is 1 because we are using the Reduced Measurement Matrix RM
		for(int i=0; i<6; i++)
			J(0, i) = (dexp[i]*L).trace();

		delete dexp;

		return J;
	}
 
	/**
	* This method assumes the data is already normalized an in the format
	  [ SRC_POINTS (Nx3) | ones (Nx1) | DST_POINTS (Nx3) | ones (Nx1) ]
	*/
	int setData(MatrixXd D)
	{
		if(D.cols() != 8)
		{
			cerr << "The format of the matrix Data is incorrect. Wrong number of columns!" << endl;
			return -1;		
		}

		RM = D.transpose()*D;

		return 0;
	}

	// size of the state
	int inputs() const { return 6; }

	// number of constraints
	int values() const { return 1; }
};


void cost_function_grad(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr)
{
	//timer.start();
	FunctorOpti* object = (FunctorOpti*)ptr;

	
	//we need to conver types
	VectorXd X(x.length());
	for(int i=0; i<x.length(); i++)
		X(i) = x[i];

	VectorXd c = object->cost(X);
	MatrixXd J = object->grad(X);

	func = c(0);
	for(int i=0; i<x.length(); i++)
		grad[i] = J(0, i);

	//timer.stop("Cost function + gradient");
}

VectorXd run_opt(MatrixXd Data, VectorXd seed, int maxits, double epsx, double epsg, double epsf)
{
	int D = seed.size();

	//initializa the functor with the input Data
	FunctorOpti* functor = new l2costFunctor();
	functor->setData(Data);

	real_1d_array x;
	x.setcontent(D, seed.data());

	mincgstate state;
	mincgreport rep;

	mincgcreate(x, state);
	mincgsetcond(state, epsg, epsf, epsx, maxits);
	mincgoptimize(state, cost_function_grad, NULL, (void*)functor);
	mincgresults(state, x, rep);

	//cout << " Iterations = " << rep.iterationscount << endl;
	//cout << " Cost function calls = " << rep.nfev << endl;

	VectorXd sol(D);
	for(int i=0; i<D; i++)
		sol(i) = x[i];

   	return sol;
}
