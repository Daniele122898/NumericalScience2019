#include <cmath>
#include <iostream>
#include <functional>
#include <vector>

double GaussLegendre5(const std::function<double(double)> &f, double a, double b) {
	const std::vector<double> c = { -0.90617984593, -0.53846931010, 0.0, 0.53846931010, 0.90617984593 };
	const std::vector<double> w = { 0.23692688505, 0.47862867049, 0.56888888888, 0.47862867049, 0.23692688505 };

	int m = c.size();
	std::vector<double> x(m);
	double d = b - a;
	for (int i = 0; i < m; ++i) {
		x[i] = d * (c[i] + 1.0) / 2.0 + a;
	}

	double q = .0;
	for (int i = 0; i < x.size(); ++i) {
		q += w[i] * f(x[i]);
	}

	return q * d / 2;
}

double CompositeSimpson(const std::function<double(double)> &f, const std::vector<double> &x) {
	int n = x.size();
	int m = n - 1;

	double q = 1.0 / 6.0 * (x[1] - x[0]) * f(x[0]);
	for (int j = 1; j < m; ++j) {
		q += 1.0 / 6.0 * (x[j + 1] - x[j - 1]) * f(x[j]);
	}
	for (int j = 1; j <= m; ++j) {
		q += 2.0 / 3.0 * (x[j] - x[j - 1]) * f((x[j] + x[j - 1]) / 2);
	}
	q += 1.0 / 6.0 * (x[m] - x[m - 1]) * f(x[m]);

	return q;
}

std::vector<double> LinSpace(int n, double a, double b) {
	std::vector<double> x(n);
	double d = (b - a) / (n - 1);
	for (int i = 0; i < n; ++i) {
		x[i] = i * d + a;
	}
	return x;
}

double f(double x) {
	return std::sqrt(x);
}

double F(double x) {
	double y = std::sqrt(x);
	return 2.0 / 3.0 * y * y * y;
}

int main() {
	int n = 5;
	double a = 0.0;
	double b = 1.0;

	std::cout << "Gauss-Legendre: " << GaussLegendre5(f, a, b) << std::endl;
	std::cout << "Simpson: " << CompositeSimpson(f, LinSpace(n, a, b)) << std::endl;
	std::cout << "Exact value: " << F(b) - F(a) << std::endl;

	return 0;
}

