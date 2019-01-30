// to compile run: g++ -std=gnu++11 cubic_spline.cpp -lmgl
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
	// we define h that is (tj-tj-1)
	VectorXd h = T.tail(n) - T.head(n); // basically the n elements counted from the right - the n elements counted from the left.

	// build the Matrix of the Linear system associated to the second derivatives. 
	// this is the tridiagonal matrix as seen on page 157
	MatrixXd A = MatrixXd::Zero(n-1,n-1);
	// let's now fill in the diagonals
	A.diagonal() = (T.segment(2, n-1) - T.segment(0, n-1))/3;
	A.diagonal(1) = h.segment(1, n-2)/6;
	A.diagonal(-1) = h.segment(1, n-2)/6;

	// build the vector of the finite differences of the data Y.
	// this must be yj -yj-1 / hj
	VectorXd slope = (Y.tail(n) - T.head(n)).cwiseQuotient(h);

	// With this we can calculate rj
	VectorXd r = slope.tail(n-1) - slope.head(n-1);

	// now lets solve this System for sigma (second derivatives)
	VectorXd sigma(n+1);
	sigma.segment(1,n-1) = A.partialPivLu().solve(r);
	sigma(0) = 0; // these are the extra boundary conditions
	sigma(n) = 0; // imposed on the cubic spline

	// now we can build the spline matrix with polynomials' coefficients
	MatrixXd spline(4, n);
	spline.row(0) = Y.head(n); // this is aj.
	// this is going to be bj
	spline.row(1) = slope - h.cwiseProduct(2*sigma.head(n)+sigma.tail(n))/6;
	// this is going to be cj
	spline.row(2) = sigma.head(n)/2;
	// this is going to be dj
	spline.row(3) = (sigma.tail(n)-sigma.head(n)).cwiseQuotient(6*h);  

	return spline;
}

VectorXd evalCubicSpline(const MatrixXd &S, const VectorXd &T, const VectorXd &evalT) {
	// Returns the values of the spline S calculated in the points X.
	// Assumes T is sorted, with no repetetions.

	// these are points between the T points of the spline that need to be interpolated
	int n = evalT.size();
	VectorXd out(n);

	// looping through EACH point evalT to see between which tj and tj+1 of the spline it is
	for (int i=0; i < n; i++) {
		// now looping through all the tj's to see between which 2 points fo the spline it is
		for (int j=0; j < T.size()-1; j++) {
			// if its smaller than tj+1 it means its between tj and tj+1 thus its in sj
			if (evalT(i) < T(j+1) || j==T.size()-2) {
				// This takes the point evalT and substracts the next point of tj that is LEFT to evalT
				// This way we get the actual distance x for that specific polynomial
				double x = evalT(i) - T(j);
				// we then calculate the polynomail aj+x*bj+x²*cj+x³*dj
				out(i) = S(0,j) + x*(S(1,j) + x*(S(2,j) + x*S(3,j)));
				// since we found the right spline polynomial we can break out of the inner loop and calculate
				// (interpolate) the next point of interest
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
	
