#include <iostream>
#include <vector>
#include <chrono>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

using namespace Eigen;

MatrixXd quadraticSpline(const VectorXd &x, const VectorXd &y) {
	int N = x.size()-1;
	MatrixXd out(3,N);

	// TODO: Solve exercise 1.d

	return out;
}

MatrixXd quadraticSplineFast(const VectorXd &x, const VectorXd &y) {
	int N = x.size()-1;
	MatrixXd out(3,N);

	// TODO: Solve exercise 1.e

	return out;
}




/***** TESTING ******/


int main() {
	const double TOL = 10e-5;
	int testnum;
	std::cin >> testnum;
	for(int i=1; i<testnum; i++){
		int N;
		std::cin >> N;
		VectorXd x(N+1), y(N+1);
		MatrixXd ans(N,3);

		for(int j=0; j<N+1; j++)
			std::cin >> x(j);
		for(int j=0; j<N+1; j++)
			std::cin >> y(j); 

		for(int j=0; j<N; j++){
			std::cin >> ans(j,0);
			std::cin >> ans(j,1);
			std::cin >> ans(j,2);
		}
		
		int flag;
		std::cin >> flag;

		auto start = std::chrono::steady_clock::now();	
		MatrixXd test;
		if(flag)
			test = quadraticSpline(x,y);
		else 
			test = quadraticSplineFast(x,y);
		auto end = std::chrono::steady_clock::now();
		double difftime = std::chrono::duration<double,std::milli>(end-start).count();

		if((test-ans).norm() > TOL * ans.norm()){
			std::cout << "Test " << i << " INCORRECT.\n";
			std::cerr << "Your answer:\n" << test << "\n\n" << "Correct answer:\n" << ans << "\n\n";
			break;
		}
		else
			std::cout << "Test " << i << " with N=" << N << " correct in " << difftime << "ms.\n";
	}
}
