#include <iostream>

#include <vector>

#include <Eigen/Sparse>
#include <Eigen/SparseLU>

using Triplet = Eigen::Triplet<double>;
using Triplets = std::vector<Triplet>;

using Vector = Eigen::VectorXd;
using Matrix = Eigen::SparseMatrix<double>;

//! \brief Efficiently construct the sparse matrix A given c, i_0 and j_0
//! \param[in] c contains entries c_i for matrix A
//! \param[in] i0 row index i_0
//! \param[in] j0 column index j_0
//! \return Sparse matrix A
Matrix buildA(const Vector & c, unsigned int i0, unsigned int j0) {
    assert(i0 > j0);
    
    unsigned int n = c.size() + 1;
    Matrix A(n,n);
    Triplets triplets;
    
    unsigned int ntriplets = 2*n;
    
    // Reserve space
    triplets.reserve(ntriplets);
    
    // Build triplets vector
    for(unsigned int i = 0; i < n; ++i) {
        triplets.push_back(Triplet(i,i,1));
        if(i < n-1) triplets.push_back(Triplet(i,i+1,c[i]));
    }
    triplets.push_back(Triplet(i0,j0,1));
    
    // Construct sparse matrix from its triplets
    
    A.setFromTriplets(triplets.begin(), triplets.end());
    return A;
}

//! \brief Solve system Ax = b with optimal complexity O(n)
//! \param[in] c Entries for matrix A
//! \param[in] b r.h.s. vector
//! \param[in] i0 index
//! \param[in] j0 index
//! \return Solution x, s.t. Ax = b
Vector solveLSE(const Vector & c, const Vector & b, unsigned int i0, unsigned int j0) {
    assert(c.size() == b.size()-1 && "Size mismatch!");
    assert(i0 > j0);
    
    //// PROBLEM 2b
    
    // Allocate solution vector
    Vector ret(b.size());
    
    unsigned int n = b.size();

    Triplets triplets;

    unsigned int ntriplets = 2*n;

    // Reserve space
    triplets.reserve(ntriplets);

    // Build triplets vector
    for(unsigned int i = 0; i < n; ++i) {
        triplets.push_back(Triplet(i,i,1));
        if(i < n-1) triplets.push_back(Triplet(i,i+1,c[i]));
    }
    
    // Construct sparse matrix from its triplets
    Matrix A(n,n);
    A.setFromTriplets(triplets.begin(), triplets.end());

    // Vectors u and v
    Vector u = Vector::Zero(n);
    Eigen::MatrixXd v = Eigen::MatrixXd::Zero(1,n);
    u( i0 ) = 1;
    v( 0, j0 ) = 1;
    
    // Apply SMW formula
    Vector Ainv_b = A.triangularView<Eigen::Upper>().solve( b );
    ret = Ainv_b - A.triangularView<Eigen::Upper>().solve( u * (v * Ainv_b)(0) ) / (1. + (v * A.triangularView<Eigen::Upper>().solve( u ))(0) );
    
    return ret;
}

int main(int, char**) {
    // Setup data for problem
    unsigned int n = 150;
    
    unsigned int i0 = 5, j0 = 4;
    
    Vector b = Vector::Random(n); // Random vector for b
    Vector c = Vector::Random(n-1); // Random vector for c
    
    //// PROBLEM 2a
    std::cout << "*** PROBLEM 2a:" << std::endl;
    
    Matrix A = buildA(c, i0, j0);
    
    // Solve sparse system using sparse LU and our own routine
    A.makeCompressed();
    Eigen::SparseLU<Matrix> splu;
    splu.analyzePattern(A); 
    splu.factorize(A);
    
    std::cout << "Error: " << std::endl << ( solveLSE(c,b,i0,j0) - splu.solve(b) ).norm() << std::endl;
}
