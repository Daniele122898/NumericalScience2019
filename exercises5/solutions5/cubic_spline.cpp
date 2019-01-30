// to compile run: g++ -std=gnu++11 CubicSplines_sol.cpp -lmgl
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <mgl2/mgl.h>

using namespace Eigen;

MatrixXd cubicSpline(const VectorXd &T, const VectorXd &Y) {
	// returns the matrix representing the spline interpolating the data
	// with abscissae T and ordinatae Y. Each column represents the coefficients
	// of the cubic polynomial on a subinterval.
	// Assumes T is sorted, has no repeated elements and T.size() == Y.size().

	int n = T.size() - 1; // T and Y have length n+1

	VectorXd h = T.tail(n) - T.head(n); // vector of lengths of subintervals (tj-tj-1)

	// build the matrix of the linear system associated to the second derivatives
	MatrixXd A = MatrixXd::Zero(n-1, n-1);
	A.diagonal()   = (T.segment(2,n-1) - T.segment(0,n-1))/3;
	A.diagonal(1)  = h.segment(1,n-2)/6;
	A.diagonal(-1) = h.segment(1,n-2)/6; 

	// build the vector of the finite differences of the data Y
	VectorXd slope = (Y.tail(n) - Y.head(n)).cwiseQuotient(h); 
				
	// right hand side vector for the system with matrix A
	VectorXd r = slope.tail(n-1) - slope.head(n-1);

	// solve the system and fill vector of second derivatives
	VectorXd sigma(n+1);
	sigma.segment(1,n-1) = A.partialPivLu().solve(r);
	sigma(0) = 0; // "simple" boundary conditions
	sigma(n) = 0; // "simple" boundary conditions

	// build the spline matrix with polynomials' coefficients
	MatrixXd spline(4, n);
	spline.row(0) = Y.head(n);
	spline.row(1) = slope - h.cwiseProduct(2*sigma.head(n) + sigma.tail(n))/6;
	spline.row(2) = sigma.head(n)/2;
	spline.row(3) = (sigma.tail(n) - sigma.head(n)).cwiseQuotient(6*h);

	return spline;
}

VectorXd evalCubicSpline(const MatrixXd &S, const VectorXd &T, const VectorXd &evalT) {
	// Returns the values of the spline S calculated in the points evalT.
	// Assumes T is sorted, with no repetitions.

	int n = evalT.size();
	VectorXd out(n);

	for (int i=0; i < n; i++) {
		for (int j=0; j < T.size()-1; j++) {
			if (evalT(i) < T(j+1) || j==T.size()-2) {
				double x = evalT(i) - T(j);
				out(i) = S(0,j) + x*(S(1,j) + x*(S(2,j) + x*S(3,j)));
				break;
			}
		}
	}

	return out;
}

int main() {
	// tests
	VectorXd T(9);
	VectorXd Y(9);
	T << 0, 0.4802, 0.7634, 1, 1.232, 1.407, 1.585, 1.879, 2;
	Y << 0., 0.338, 0.7456, 0, -1.234, 0 , 1.62, -2.123, 0;

	int len = 1 << 9;
	VectorXd evalT = VectorXd::LinSpaced(len, T(0), T(T.size()-1));

	VectorXd evalSpline = evalCubicSpline(cubicSpline(T, Y), T, evalT);

 	mglData datx, daty;
	datx.Link(evalT.data(), len);
	daty.Link(evalSpline.data(), len);
	mglGraph gr;
	gr.SetRanges(0, 2, -3, 3);
	gr.Plot(datx, daty, "0");
	gr.WriteFrame("spline.eps");
}
	
