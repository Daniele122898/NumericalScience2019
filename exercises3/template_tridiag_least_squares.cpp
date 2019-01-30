#include <iostream>

#include <eigen3/Eigen/Dense>

/* @brief 
 * @param[in] z An $n$-dimensional vector containing one side of input data
 * @param[in] c An $n$-dimensional vector containing the other side of input data
 * @param[out] x The vector of parameters $(\alpha,\beta)$, intercept and slope of the line fitted
 */
Eigen::Vector2d lsqEst(const Eigen::VectorXd &z, const Eigen::VectorXd &c) {
	Eigen::Vector2d x;
	int n = z.rows();
	Eigen::MatrixXd A(n, 2);
	// First lets construct A
	for (unsigned int i=0; i<n; i++) {
		A(i,0) = z(i);
		if(i==0) {
			A(i,1) = z(i+1);
		} else if (i == n-1) {
			A(i,1) = z(i-1);
		} else {
			A(i,1) = z(i-1)+z(i+1);
		}
	}
	// now lets solve this shit
	x = (A.transpose()*A).llt().solve(A.transpose()*c);

	return x;
}

/* @brief 
 * @param[in] z An $n$-dimensional vector containing one side of input data
 * @param[in] c An $n$-dimensional vector containing the other side of input data
 * @param[out] x The vector of parameters $(\alpha,\beta)$, intercept and slope of the line fitted
 */
Eigen::Vector2d lsqEstCorrect(const Eigen::VectorXd &z, const Eigen::VectorXd &c) {
	Eigen::Vector2d x;

	int n = z.size();
	Eigen::MatrixXd A(n,2);
	A.col(0) = z;
	A(0,1) = z(1);
	for(int i = 1; i < n - 1; ++i) {
		A(i, 1) = z(i - 1) + z(i + 1);
	}
	A(n - 1, 1) = z(n - 2);
	
	// Normal equations
	Eigen::MatrixXd lhs = A.transpose() * A;
	Eigen::VectorXd rhs = A.transpose() * c;
	x = lhs.fullPivLu().solve(rhs);

	return x;
}

int main() {
    int n = 10;
    Eigen::VectorXd z(n), c(n);
    for(int i = 0; i < n; ++i) {
		z(i) = i + 1;
		c(i) = n - i;
	}

	Eigen::Vector2d x = lsqEst(z, c);
	Eigen::Vector2d x2 = lsqEstCorrect(z, c);

	std::cout << "alpha = " << x(0) << std::endl;
	std::cout << "beta = "  << x(1) << std::endl;
	std::cout << "alpha2 = " << x2(0) << std::endl;
	std::cout << "beta2 = "  << x2(1) << std::endl;

	return 0;
}
