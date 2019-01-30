// compile with: g++ template_structured_matrix_vector.cpp -I/usr/include/eigen3 -lmgl

#include <chrono>
#include <iostream>
#include <iomanip>
#include <limits>
#include <ratio>
#include <vector>

#include <eigen3/Eigen/Dense>
#include <mgl2/mgl.h>

using namespace Eigen;

/* \brief compute $\mathbf{A}\mathbf{x}$
 * \mathbf{A} is defined by $(\mathbf{A})_{i,j} := \min {i,j}$
 * \param[in] x vector x for computation of A*x = y
 * \param[out] y = A*x
 */
void multAminSlow(const VectorXd & x, VectorXd & y) {
    unsigned int n = x.size();

    VectorXd one = VectorXd::Ones(n);
    VectorXd linsp = VectorXd::LinSpaced(n,1,n);
    y = ( ( one * linsp.transpose() )
          .cwiseMin( linsp * one.transpose()) ) * x;
}

/* \brief compute $\mathbf{A}\mathbf{x}$
 * \mathbf{A} is defined by $(\mathbf{A})_{i,j} := \min {i,j}$
 * Instead of a "Matlab style" construcion of the product,
 * we use simple loops.
 * \param[in] x vector x for computation of A*x = y
 * \param[out] y = A*x
 */
void multAminLoops(const VectorXd & x, VectorXd & y) {
    unsigned int n = x.size();

    MatrixXd A(n,n);

    for(unsigned int i = 0; i < n; ++i) {
        for(unsigned int j = 0; j < n; ++j) {
            A(i,j) = std::min(i+1,j+1);
        }
    }
    y = A * x;
}

/* \brief compute $\mathbf{A}\mathbf{x}$
 * This function has optimal complexity.
 * \mathbf{A} is defined by $(\mathbf{A})_{i,j} := \min {i,j}$
 * \param[in] x vector x for computation of A*x = y
 * \param[out] y = A*x
 */
void multAmin(const VectorXd & x, VectorXd & y) {
	unsigned int n = x.size();

    if (n == 0) return;
    
    y = VectorXd::Zero(n);
    VectorXd v = VectorXd::Zero(n);

    v(n-1) = x(n-1);
    for (unsigned int i = n-1; i>0; --i) {
        v(i-1) = v(i) + x(i-1);
    } 

    y(0) = v(0);
    for (unsigned int i=1; i<n; ++i) {
        y(i) = y(i-1) + v(i);
    }

}

int main(void) {
    // Testing correctness of the code
    unsigned int M = 10;
    VectorXd xa = VectorXd::Random(M);
    VectorXd ys, yf;

    /* STUFF TO DELETE
    unsigned int n1 = xa.size();
    VectorXd one = VectorXd::Ones(n1);
    VectorXd linsp = VectorXd::LinSpaced(n1,1,n1);
    auto y1 = one * linsp.transpose();
    std::cout << y1 << std::endl;
    std::cout << "\n\n" << std::endl;
    std::cout << (linsp * one.transpose()) << std::endl;
    std::cout << "\n\n" << std::endl;
    std::cout << y1.cwiseMin( linsp * one.transpose()) << std::endl;
    std::cout << "\n\n" << std::endl;
    std::cout << (y1.cwiseMin( linsp * one.transpose())) * xa << std::endl;
    DELETE */

    multAmin(xa, yf);
    multAminSlow(xa, ys);

    // Error should be small
    std::cout << "||ys-yf|| = " << (ys - yf).norm() << std::endl;

    unsigned int nLevels = 9;
	unsigned int *n = new unsigned int[nLevels];
	double *minTime = new double[nLevels];
	double *minTimeLoops = new double[nLevels];
	double *minTimeEff = new double[nLevels];

	n[0] = 4;
	for (unsigned int i=1; i<nLevels; i++)
		n[i] = 2*n[i-1];

    unsigned int nruns = 10;

    std::cout << "--> Timings:" << std::endl;
    // Header, see iomanip documentation
    std::cout << std::setw(15)
              << "N"
              << std::scientific << std::setprecision(3)
              << std::setw(15) << "multAminSlown"
              << std::setw(15) << "multAminLoops"
              << std::setw(15) << "multAmin"
              << std::endl;
    // From $2^4$ to $2^{13}$
    for(unsigned int i = 0; i<nLevels; i++) {
        // Compute runtime many times
    	double min_slow = std::numeric_limits<double>::infinity();
    	double min_slow_loops = std::numeric_limits<double>::infinity();
		double min_fast = std::numeric_limits<double>::infinity();

        for(unsigned int r = 0; r < nruns; ++r) {
            VectorXd x = VectorXd::Random(n[i]);
            VectorXd y;

        	std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
        	std::chrono::duration<double, std::ratio<1>> duration;

            // Runtime of slow method
            start = std::chrono::high_resolution_clock::now();
            multAminSlow(x, y);
            end = std::chrono::high_resolution_clock::now();
            duration = end - start;
            min_slow = std::min(min_slow, (double)duration.count());

            // Runtime of slow method with loops
            start = std::chrono::high_resolution_clock::now();
            multAminLoops(x, y);
            end = std::chrono::high_resolution_clock::now();
            duration = end - start;
            min_slow_loops = std::min(min_slow_loops, (double)duration.count());

            // Runtime of fast method
            start = std::chrono::high_resolution_clock::now();
            multAmin(x, y);
            end = std::chrono::high_resolution_clock::now();
            duration = end - start;
            min_fast = std::min(min_fast, (double)duration.count());
        }
        
        minTime[i] = min_slow;
        minTimeLoops[i] = min_slow_loops;
        minTimeEff[i] = min_fast;

        std::cout << std::setw(15)
                  << n[i]
                  << std::scientific << std::setprecision(3)
                  << std::setw(15) << minTime[i]
                  << std::setw(15) << minTimeLoops[i]
                  << std::setw(15) << minTimeEff[i]
                  << std::endl;
    }
    
    // Plotting with MathGL
    double nMgl[nLevels];
    double ref1[nLevels], ref2[nLevels];
    for (int i=0; i<nLevels; i++) {
    	nMgl[i] = n[i];
    	ref1[i] = 1e-8*pow(n[i],2);
    	ref2[i] = 1e-7*n[i];
    }
    
    mglData matSize;
    matSize.Link(nMgl, nLevels);
    
    mglData data1, data2;
    mglData dataRef1, dataRef2;
  	data1.Link(minTime, nLevels);
  	data2.Link(minTimeEff, nLevels);
  	dataRef1.Link(ref1,nLevels);
  	dataRef2.Link(ref2,nLevels);
  	
  	mglGraph *gr = new mglGraph;
    gr->Title("Runtime of multAmin");
  	gr->SetRanges(n[0],n[0]*pow(2,nLevels-1),1e-6,1e+1);  gr->SetFunc("lg(x)","lg(y)");
  	gr->Axis();
  	gr->Plot(matSize,data1,"k +"); gr->AddLegend("slow","k +");
  	gr->Plot(matSize,data2,"r +"); gr->AddLegend("efficient","r +");
  	gr->Plot(matSize,dataRef1,"k"); gr->AddLegend("O(n^2)","k");
  	gr->Plot(matSize,dataRef2,"r"); gr->AddLegend("O(n)","r");
  	gr->Label('x',"Matrix size [n]",0);
  	gr->Label('y', "Runtime [s]",0);
    gr->Legend(2);
	gr->WriteFrame("multAmin_comparison.eps");


    // The following code is just for demonstration purposes.
    // Build Matrix B with dimension 10x10
    unsigned int nn = 10;
    MatrixXd B = MatrixXd::Zero(nn,nn);
    for(unsigned int i = 0; i < nn; ++i) {
        B(i,i) = 2;
        if(i < nn-1) B(i+1,i) = -1;
        if(i > 0) B(i-1,i) = -1;
    }
    B(nn-1,nn-1) = 1;

    std::cout << std::fixed;
    
    std::cout << "B =" << std::endl
        << B << std::endl;

    VectorXd sol = VectorXd::Zero(nn);
    
    for (int i=0; i<nn; i++) {
        VectorXd unitvec = VectorXd::Unit(nn, i);

        multAmin(B*unitvec, sol);

        std::cout << "y = " << std::endl
            << sol << std::endl;
        sol = VectorXd::Zero(nn);
    }

    std::cout << " A*B*ej gives us the same unit vector meaning that B is the inverse of A and thus A*B is equal to I." << std::endl;

}
