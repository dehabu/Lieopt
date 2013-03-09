#include "lieopt.hpp"

pair<MatrixXd, Matrix4d> normalizeHomoData(MatrixXd D)
{
	int N = D.rows();
	MatrixXd Y(N, 4);

	RowVectorXd T = D.colwise().mean();
	
	double avgDist = 0.0;
	for(int i=0; i<N; i++)
		avgDist += (D.row(i) - T).norm();
	avgDist /= double(N);

	double sc = sqrt(3) / avgDist;

	//transformation for normalizing points
	Matrix4d Tr;
	Tr << 	sc,  0,   0,   -sc*T(0),
		0,   sc,  0,   -sc*T(1),
		0,   0,   sc,  -sc*T(2),
		0,   0,   0,    1;	

	Vector4d aux, aux2;
	for(int i=0; i<N; i++)
	{
		aux.head<3>() = D.row(i);
		aux(3) = 1.0;
		aux2 = Tr*aux;

		aux2(0) /= aux2(3);
		aux2(1) /= aux2(3);
		aux2(2) /= aux2(3);
		aux2(3) /= aux2(3);

		Y.row(i) = aux2;
	}

	return make_pair<MatrixXd, Matrix4d>(Y, Tr);	
}

Matrix4d denormalizeTransformation(const Matrix4d Tr, const Matrix4d Tsrc, const Matrix4d Tdst)
{
	return Tdst.inverse()*Tr*Tsrc;
}

Vector4d disambiguateQuaternion(Vector4d q)
{
	double lastMax = -1;
	double sign = 0.0;
	for(int i=0; i<4; i++)
	{
		if(abs(q(i)) > lastMax)
		{
			lastMax = abs(q(i));
			sign = (q(i) > 0)? 1.0 : -1.0;
		}
	}

	return q*sign;
}

Matrix4d horn_fit(MatrixXd Xsrc, MatrixXd Xdst)
{
	int N = Xsrc.rows();
	MatrixXd XlN(N, 3);
	MatrixXd XrN(N, 3);

	//demeanization of points
	RowVectorXd mean_src = Xsrc.colwise().mean();
	RowVectorXd mean_dst = Xdst.colwise().mean();

	for(int i=0; i<N; i++)
	{
		Vector3d v = Xsrc.row(i) - mean_src;
		XlN.row(i) = v.head<3>();

		v = Xdst.row(i) - mean_dst;
		XrN.row(i) = v.head<3>();
	}

	double S_XX = 0, S_XY = 0, S_XZ = 0;
	double S_YX = 0, S_YY = 0, S_YZ = 0;
	double S_ZX = 0, S_ZY = 0, S_ZZ = 0;

	for(int i=0; i<N; i++)
	{
		S_XX = S_XX + (XlN(i, 0)*XrN(i, 0));
		S_XY = S_XY + (XlN(i, 0)*XrN(i, 1));
		S_XZ = S_XZ + (XlN(i, 0)*XrN(i, 2));

		S_YX = S_YX + (XlN(i, 1)*XrN(i, 0));
		S_YY = S_YY + (XlN(i, 1)*XrN(i, 1));
		S_YZ = S_YZ + (XlN(i, 1)*XrN(i, 2));

		S_ZX = S_ZX + (XlN(i, 2)*XrN(i, 0));
		S_ZY = S_ZY + (XlN(i, 2)*XrN(i, 1));
		S_ZZ = S_ZZ + (XlN(i, 2)*XrN(i, 2));	
	}

	Matrix4d W;
	W << 	S_XX + S_YY + S_ZZ,	S_YZ - S_ZY, 		S_ZX - S_XZ,		S_XY - S_YX,		
		S_YZ - S_ZY,		S_XX - S_YY - S_ZZ,	S_XY + S_YX,		S_ZX + S_XZ,
		S_ZX - S_XZ,		S_XY + S_YX,		-S_XX + S_YY - S_ZZ,	S_YZ + S_ZY,
		S_XY - S_YX,		S_ZX + S_XZ,		S_YZ + S_ZY,		-S_XX - S_YY + S_ZZ;

	//let's calculate an eigen decomposition
	SelfAdjointEigenSolver<Matrix4d> eigensolver(W);
	if (eigensolver.info() != Success)
		abort();

	Vector4d D = eigensolver.eigenvalues();
	Matrix4d V = eigensolver.eigenvectors();

	//Gets eigenvector corresponding to maximum eigenvalue
	Vector4d q = V.rightCols<1>();
	Vector4d q2 = disambiguateQuaternion(q);

	//Now the quaternion has to be mapped to an orthogonal matrix
	Vector4d q2norm = q2 / q2.norm();

	Matrix3d Z;
	Vector3d v = q2norm.tail<3>();

	Z << 	q2norm(0), -q2norm(3), q2norm(2),
		q2norm(3), q2norm(0), -q2norm(1),
		-q2norm(2), q2norm(1), q2norm(0);

	Vector3d mu_l = mean_src.head<3>();
	Vector3d mu_r = mean_dst.head<3>();

	Matrix3d R = v*v.transpose() + Z*Z;
	Vector3d T = mu_r - R*mu_l;

	Matrix4d result;
	result.block(0, 0, 3, 3) = R;
	result.block(0, 3, 3, 1) = T;
	result(3, 0) = 0;
	result(3, 1) = 0;
	result(3, 2) = 0;
	result(3, 3) = 1;

	return result;
}


Matrix4d recalculateTrInliersSE3(Matrix4d Tr, MatrixXd Xsrc, MatrixXd Xdst, double threshold)
{
	int N = Xsrc.rows();
	int inliers = 0;

	//I'm gonna do this in a nasty way...
	MatrixXd Xsrc_in(N, 4);
	MatrixXd Xdst_in(N, 4);

	for(int i=0; i<N; i++)
	{
		Vector4d v = Tr*Xsrc.row(i).transpose();
		v(0) /= v(3);
		v(1) /= v(3);
		v(2) /= v(3);
		v(3) /= v(3);

		double nn = (v - Xdst.row(i).transpose()).norm();
		if(nn < threshold)
		{
			Xsrc_in.row(inliers) = Xsrc.row(i);
			Xdst_in.row(inliers) = Xdst.row(i);
			inliers++;
		}
	}

	if(inliers < 3)
		return Tr;
	else
	{

		MatrixXd Xsrc_in2 = Xsrc_in.block(0, 0, inliers, 3);
		MatrixXd Xdst_in2 = Xdst_in.block(0, 0, inliers, 3);

		return horn_fit(Xsrc_in2, Xdst_in2);
	}
}

MatrixXd lieOpt(MatrixXd Xsrc, MatrixXd Xdst, double threshold, int maxits, double epsx, double epsg, double epsf)
{
	int N = Xsrc.rows();
	MatrixXd DataNorm(N, 8);
	
	//Timerx timer;
	
	//timer.start();
	pair<MatrixXd, Matrix4d> normXsrc = normalizeHomoData(Xsrc);
	//timer.stop("Normalization 1");

	//timer.start();
	pair<MatrixXd, Matrix4d> normXdst = normalizeHomoData(Xdst);
	//timer.stop("Normalization 2");

	//timer.start();
	DataNorm.block(0, 0, N, 4) = normXsrc.first;
	DataNorm.block(0, 4, N, 4) = normXdst.first;
	//timer.stop("Block copy data");
	
	VectorXd seed = VectorXd::Zero(6);

	//timer.start();
	VectorXd sol_t = run_opt(DataNorm, seed, maxits, epsx, epsg, epsf);
	//timer.stop("Run_opt");

	//timer.start();
	Matrix4d Tr = se3Exp(sol_t);
	//timer.stop("SE3Exp");

	//timer.start();
	Matrix4d Tr_in = recalculateTrInliersSE3(Tr, normXsrc.first, normXdst.first, threshold);
	//timer.stop("Recalculate Tr inliers");

	//timer.start();
	Matrix4d Tr2 = denormalizeTransformation(Tr_in, normXsrc.second, normXdst.second);
	//timer.stop("Denormalize Tr");

	return Tr2;

}


