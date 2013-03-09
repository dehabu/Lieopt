#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include "lieopt.hpp"
#include "mytime.hpp"

using namespace std;
using namespace Eigen;

pair<MatrixXd, MatrixXd> readData(const char* file)
{
	ifstream fd(file);
	if(!fd.is_open())
	{
		cerr << "Error reading the input file" << endl;
		throw -2;
	}

	int N;	
	//number of matches
	fd >> N;

	MatrixXd Xsrc(N, 3);
	MatrixXd Xdst(N, 3);
	for(int i=0; i<N; i++)
	{
		fd >> Xsrc(i, 0);
		fd >> Xsrc(i, 1);
		fd >> Xsrc(i, 2);

		fd >> Xdst(i, 0);
		fd >> Xdst(i, 1);
		fd >> Xdst(i, 2);
	}

	fd.close();

	return make_pair<MatrixXd, MatrixXd>(Xsrc, Xdst);
}

int main(int argc, char **argv)
{
	Timerx timer;
	double threshold = 0.3;

	if(argc != 2)
	{
		cerr << "Calling error. You must provide a data file, e.g.: " << argv[0] << " datafile.txt" << endl;
		return -1;
	}

	pair<MatrixXd, MatrixXd> matches = readData(argv[1]);
	MatrixXd Xsrc = matches.first;
	MatrixXd Xdst = matches.second;

	//add other parameters later
	timer.start();
	MatrixXd sol = lieOpt(Xsrc, Xdst, threshold, 100, 1e-6, 1e-6, 1e-6);
	timer.stop("Total time");

	cout << "The solution found is: " << endl << sol << endl << endl;
	
	return 0;	   
}
