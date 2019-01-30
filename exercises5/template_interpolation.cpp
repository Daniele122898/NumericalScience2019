#include <iostream>

#include <eigen3/Eigen/Dense>

struct Newton {
	Newton(const Eigen::VectorXd &x) : _x(x), _a(x.size()) { }
	void Interpolate(const Eigen::VectorXd &y);
	double operator()(double x) const;

public:
	Eigen::VectorXd geta() {
		return _a;
	}

private:
	Eigen::VectorXd _x;	// nodes
	Eigen::VectorXd _a;	// coefficients
};

// Compute the coefficients in the Newton basis.
void Newton::Interpolate(const Eigen::VectorXd &y) {
	// compute _a
	// First lets compute A for Aa=y
	int n = y.size();
	Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n,n);
	A.col(0) = Eigen::VectorXd::Ones(n);
	for (int i=1; i<n; i++) {
		Eigen::VectorXd temp = Eigen::VectorXd::Ones(n);
		temp.head(i) = Eigen::VectorXd::Zero(i);
		for (int j=i; j<n; j++) {
			temp(j) = _x(j)-_x(i-1);
		}
		A.col(i) = A.col(i-1).cwiseProduct(temp);
	}
	// Now we can solve for _a
	_a = A.partialPivLu().solve(y);
}

// Evaluate the interpolant at x.
double Newton::operator()(double x) const {
	int n = _a.size();
	double y = _a(n-1);
	for (int i=n-2; i >= 0; --i) {
		y = y*(x-_x(i))+_a(i);
	}
	return y; 
}

struct Lagrange {
	Lagrange(const Eigen::VectorXd &x);
	void Interpolate(const Eigen::VectorXd &y) { _y = y; }
	double operator()(double x) const;

public:
Eigen::VectorXd getl() {
	return _l;
}

private:
	Eigen::VectorXd _x;	// nodes
	Eigen::VectorXd _l;	// weights
	Eigen::VectorXd _y;	// coefficients
};

// Compute the weights l for given nodes x.
Lagrange::Lagrange(const Eigen::VectorXd &x) : _x(x), _l(x.size()), _y(x.size()) {
	int n = _l.size();
	for (int i=0; i < n; i++) {
		_l(i) = 1;
		for (int j=0; j < n; j++) {
			if (i==j) continue;
			_l(i) *= 1.0/(_x(i) - _x(j));
		}
	}
}

// Evaluate the interpolant at x.
double Lagrange::operator()(double x) const {
	// first lets calculate the w(x):=prod(0->n, x-xj)
	int n = _x.size();
	double w = 1;
	for (int i=0; i<n; i++) {
		w *= (x-_x(i));
	}
	// now lets calculate Li(x)= w(x)*(li/(x-xi))
	Eigen::VectorXd Li(n);
	for (int i=0; i<n; i++) {
		Li(i) = w*(_l(i)/(x-_x(i)));
	}
	// now the lagrange interpolant q(x) := sum(0->n, yi*Li(x))
	double result = 0;
	for(int i=0; i<n; i++) {
		result += _y(i)*Li(i);
	}
	return result;
}

// Runge function
Eigen::VectorXd r(const Eigen::VectorXd &x) {
	return (1.0 / (1.0 + 25.0 * x.array() * x.array())).matrix();
}

int main() {
	int n = 5;
	Eigen::VectorXd x;
	x.setLinSpaced(5, -1.0, 1.0);
	Eigen::VectorXd y = r(x);

	Newton p(x);
	p.Interpolate(y); // correct result: p._a = [0.0384615, 0.198939, 1.5252, -3.31565, 3.31565]
	std::cout << "p._a is = " << std::endl
		<< p.geta().transpose() << std::endl
		<< "0.0384615, 0.198939,   1.5252, -3.31565,  3.31565 correct result" << std::endl;

	Lagrange q(x);    // correct result: p._l = [0.666667, -2.66667, 4, -2.66667, 0.666667]
	q.Interpolate(y);

	std::cout << "q._l is = " << std::endl
		<< q.getl().transpose() << std::endl
		<< "0.666667 -2.66667        4 -2.66667 0.666667 correct result" << std::endl;

	// Compute difference of p and q.
	int m = 22;
	double offset = 0.08333333333;
	x.setLinSpaced(m, -1.0 + offset, 1.0 - offset);
	double norm2 = .0;
	for (int i = 0; i < m; ++i) {
		double d = p(x(i)) - q(x(i));
		norm2 += d * d;
	}

	// By uniquenss of the interpolation polynomial, we expect p = q.
	std::cout << "This number should be close to zero: " << norm2 << std::endl;

	return 0;
}

