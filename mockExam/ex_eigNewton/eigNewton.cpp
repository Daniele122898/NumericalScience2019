#include <iostream>

#include <Eigen/Dense>

void eigNewton(const Eigen::MatrixXd &A, double tol, int maxItr, Eigen::VectorXd &z) {
	// TODO: Task (c)
}

int main() {
	int n = 2;
	int m = n + 1;

	Eigen::VectorXd x = Eigen::VectorXd::Ones(n);
	x << 1.0, 0.0;
	double lambda = 0.0;

	Eigen::VectorXd z(m);
	z.head(n) = x;
	z(n) = lambda;

	Eigen::MatrixXd A = Eigen::MatrixXd::Ones(n, n);
	A << 2, 1, 1, 2;

	eigNewton(A, .0001, 10, z);

	std::cout << "eigenvector:" << std::endl;
	std::cout << z.head(n) << std::endl;
	std::cout << "eigenvalue: " << z(n) << std::endl;

	return 0;
}
