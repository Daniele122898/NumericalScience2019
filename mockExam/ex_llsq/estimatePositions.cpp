#include <iostream>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

using namespace Eigen;

//!
//! \brief buildDistanceLSQMatrix Efficiently build the system matrix
//! mapping positions to distances.
//! \param n Number of points (including $x_1$).
//! \return The system matrix $A$.
//!

SparseMatrix<double> buildDistanceLSQMatrix(int n) {
    SparseMatrix<double> A(n*(n-1)/2, n-1);
    // list of non-zero coefficients.
    std::vector<Triplet<double>> triplets;
    // two non-zeros per row at most 2n(n-1)/2 but  first n-1 rows
    // only 1 entry. ->(n-1)(n-1) entries
    triplets.reserve((n-1)*(n-1)); 
    // loop over vertical blocks. we removed the first column
    int row =0;
    for (int i=0; i<n-1; i++) {
        for(int j = i; j < n-1; j++)
        {
            // this will create the Identity matrix by going down the rows
            // and inserting 1 always shifted right by one
            triplets.push_back(Triplet<double>(row, j, 1));
            if (i>0) { // since we removed the first column we need to check this
                // This adds the -1 always at the current i-1 index meaning 1 behind 
                // the identity matrix and then just going down the rows
                triplets.push_back(Triplet<double>(row, i-1, -1));
            } 
            row++; // go to the next row.
        }
        
    }
    // Build matrix
    A.setFromTriplets(triplets.begin(), triplets.end());
    A.makeCompressed();
    return A;
}


//!
//! \brief estimatePointsPositions Return positions (without $x_1$).
//! The coordinate $x_1$ is assumed $x_1 = 0$.
//! \param D An $n \times n$ anti-symmetric matrix of distances.
//! \return Vector of positions $x_2, \dots, x_n$.
//!

VectorXd estimatePointsPositions(const MatrixXd& D) {
    VectorXd x;
    // Vector of sum of columns of A
    ArrayXd b = D.rowwise().sum().tail(D.cols()-1);
    // Vector 1
    ArrayXd one = ArrayXd::Constant(D.cols()-1, 1);
    // Apply SMW formula
    x =  (b + one * b.sum()) / D.cols();
    return x;
}


int main() {
    
    int n = 5;
    
    // PART 1: build and print system matrix A
    auto A = buildDistanceLSQMatrix(n);
    std::cout << "**** PART 1 ****" << std::endl;
    std::cout << "The Matrix A is:"
              << std::endl
              << A
            //   << std::endl
            //   << A.transpose()
            //   << std::endl
            //   << A.transpose() * A
              << std::endl;

    // PART 2: solve the LSQ system and find positions
    std::cout << "**** PART 2 ****" << std::endl;
    // Vector of positions
    n = 5;

    // Build D
    MatrixXd D(n,n);
    D << 0,   -2.1, -3,   -4.2, -5,
         2.1,  0,   -0.9, -2.2, -3.3,
         3,    0.9,  0,   -1.3, -1.1,
         4.2,  2.2,  1.3,  0,   -1.1,
         5,    3.3,  1.1,  1.1,  0;
    std::cout << "The matrix D is:"
              << std::endl
              << D
              << std::endl;

    // Find out positions
    VectorXd x_recovered = estimatePointsPositions(D);
    std::cout << "The positions [x_2, ..., x_n] obtained from the LSQ system are:"
              << std::endl
              << x_recovered
              << std::endl;
}
