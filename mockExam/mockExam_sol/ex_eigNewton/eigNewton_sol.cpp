#include <iostream>

#include <Eigen/Dense>

void eigNewton(const Eigen::MatrixXd &A, double tol, int maxItr, Eigen::VectorXd &z) {
	int m = z.size();
	int n = m - 1;

	Eigen::MatrixXd DF(m, m);
	Eigen::VectorXd F(m);
	Eigen::VectorXd F_old(m);

	for (int i = 0; i < maxItr; ++i) {
		Eigen::VectorXd x = z.head(n);
		F.head(n) = A * x - z(n) * x;
		F(n) = 1.0 - 0.5 * x.squaredNorm();

		if (F.squaredNorm() < tol) {
			std::cout << "tol reached with i = " << i << std::endl;
			return;
		}

		DF.topLeftCorner(n, n) = A - z(n) * Eigen::MatrixXd::Identity(n,  n);
		DF.col(n) = -z;
		DF.row(n) = -z.transpose();
		DF(n, n) = 0;

		z += -DF.fullPivLu().solve(F);
	}

	std::cout << "maxItr reached" << std::endl;
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
