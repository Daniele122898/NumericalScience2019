#include <iostream>

#include <Eigen/Dense>

struct Newton {
	Newton(const Eigen::VectorXd &x) : _x(x), _a(x.size()) { }
	void Interpolate(const Eigen::VectorXd &y);
	double operator()(double x) const;

private:
	Eigen::VectorXd _x;	// nodes
	Eigen::VectorXd _a;	// coefficients
};

// Compute the coefficients in the Newton basis.
void Newton::Interpolate(const Eigen::VectorXd &y) {
	_a = y;
	int n = _a.size();
	for (int j = 0; j < n - 1; ++j) {
		for (int i = n - 1; i > j; --i) {
			_a(i) = (_a(i) - _a(i - 1)) / (_x(i) - _x(i - 1 - j));
		}
	}
}

// Evaluate the interpolant at x.
double Newton::operator()(double x) const {
	int n = _a.size();
	double y = _a(n - 1);
	for (int i = n - 2; i >= 0; --i) {
		y = y * (x - _x(i)) + _a(i);
	}
	return y;
}

struct Lagrange {
	Lagrange(const Eigen::VectorXd &x);
	void Interpolate(const Eigen::VectorXd &y) { _y = y; }
	double operator()(double x) const;

private:
	Eigen::VectorXd _x;	// nodes
	Eigen::VectorXd _l;	// weights
	Eigen::VectorXd _y;	// coefficients
};

// Compute the weights l for given nodes x.
Lagrange::Lagrange(const Eigen::VectorXd &x) : _x(x), _l(x.size()), _y(x.size()) {
	int n = _x.size();
	for (int j = 0; j < n; ++j) {
		double dw = 1.0;
		for (int i = 0; i < n; ++i) {
			if (i != j) dw *= _x(j) - _x(i);
		}
		_l(j) = 1.0 / dw;
	}
}

// Evaluate the interpolant at x.
double Lagrange::operator()(double x) const {
	int n = _x.size();
	Eigen::VectorXd L(n);
	double wx = 1.0;
	for (int i = 0; i < n; ++i) {
		wx *= x - _x(i);
	}
	for (int i = 0; i < n; ++i) {
		L(i) = wx * _l(i) / (x - _x(i));
	}
	return _y.dot(L);
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

	Lagrange q(x);    // correct result: p._l = [0.666667, -2.66667, 4, -2.66667, 0.666667]
	q.Interpolate(y);

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

