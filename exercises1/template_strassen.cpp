#include<bits/stdc++.h>
#include<eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

MatrixXf mult(const MatrixXf &A, const MatrixXf &B) {
	int N = A.rows();
	MatrixXf C = MatrixXf::Zero(N, N);

	for (unsigned int i=0; i<N; i++) {
		for (unsigned int j=0; j<N; j++) {
			for (unsigned int k=0; k<N; k++) {
				C(i,j) += A(i,k) * B(k,j);
			}
		}
	}

	return C;
}

MatrixXf mult2(const MatrixXf &A, const MatrixXf &B) {
	int N = A.rows();
	MatrixXf C = MatrixXf::Zero(N, N);

	for (unsigned int i=0; i<N; i++) {
		for (unsigned int j=0; j<N; j++) {
			for (unsigned int k=0; k<N; k++) {
				C(j,i) += A(j,k) * B(k,i);
			}
		}
	}

	return C;
}

#define b11 .block(0, 0, N/2, N/2)
#define b12 .block(0, N/2, N/2, N/2)
#define b21 .block(N/2, 0, N/2, N/2)
#define b22 .block(N/2, N/2, N/2, N/2)

MatrixXf mult_rec(const MatrixXf &A, const MatrixXf &B) {
	int N = A.rows();

	if (N <= 1) return A*B;

	MatrixXf C(N, N);

	C << 
		mult_rec(A b11, B b11) + mult_rec(A b12, B b21),
		mult_rec(A b11, B b12) + mult_rec(A b12, B b22),
		mult_rec(A b21, B b11) + mult_rec(A b22, B b21),
		mult_rec(A b21, B b12) + mult_rec(A b22, B b22);

	return C;
}
	

MatrixXf strassen(const MatrixXf &A, const MatrixXf &B) {
	int N = A.rows();

	if (N <= 1) return A*B;

	MatrixXf C(N, N);

	MatrixXf M1 = strassen(A b11 + A b22, B b11 + B b22);
	MatrixXf M2 = strassen(A b21 + A b22, B b11);
	MatrixXf M3 = strassen(A b11, B b12 - B b22);
	MatrixXf M4 = strassen(A b22, B b21 - B b11);
	MatrixXf M5 = strassen(A b11 + A b12, B b22);
	MatrixXf M6 = strassen(A b21 - A b11, B b11 + B b12);
	MatrixXf M7 = strassen(A b12 - A b22, B b21 + B b22);

	C << 
		M1 + M4 - M5 + M7,
		M3 + M5,
		M2 + M4,
		M1 - M2 + M3 + M6;

	return C;
}

int main() {
	srand(time(0));
	cout << setprecision(6) << setfill(' ');

	for (int i = 1; i < 9; i++) {
		int N = 1 << i;		
		cout << "Matrix size = " << N << endl;
		MatrixXd AA = MatrixXd::Random(N, N);
		MatrixXd BB = MatrixXd::Random(N, N);
		MatrixXd ans = AA*BB;
		MatrixXf A = AA.cast<float>();
		MatrixXf B = BB.cast<float>();

		auto start = std::chrono::steady_clock::now();
		MatrixXf W = mult(A, B);
		auto finish = std::chrono::steady_clock::now();
		cout << setw(24) << " " <<  setw(15) 
			<< "Time (s)" << setw(20) << "Error (l2-norm)"  << endl;
		cout << setw(24) << "Naive iterative "<< setw(15) 
			<< std::chrono::duration_cast<std::chrono::duration<double> > 
			(finish - start).count() << setw(20) << (W.cast<double>() - ans).norm() << endl;

		start = std::chrono::steady_clock::now();
		MatrixXf P = mult2(A, B);
		finish = std::chrono::steady_clock::now();
		cout << setw(24) << "Naive iterative 2"<< setw(15) 
			<< std::chrono::duration_cast<std::chrono::duration<double> > 
			(finish - start).count() << setw(20) << (P.cast<double>() - ans).norm() << endl;
	
		start = std::chrono::steady_clock::now();
		MatrixXf X = mult_rec(A, B);
		finish = std::chrono::steady_clock::now();
		cout << setw(24) << "Naive recursive "<< setw(15) 
			<< std::chrono::duration_cast<std::chrono::duration<double> > 
			(finish - start).count() << setw(20) << (X.cast<double>() - ans).norm() << endl;
		start = std::chrono::steady_clock::now();
		MatrixXf Y = strassen(A, B);
		finish = std::chrono::steady_clock::now();
		cout << setw(24) << "Strassen recursive " << setw(15) 
			<< std::chrono::duration_cast<std::chrono::duration<double> >
			(finish - start).count() << setw(20) << (Y.cast<double>() - ans).norm() << endl;
		start = std::chrono::steady_clock::now();	
		MatrixXf Z = A*B;
		finish = std::chrono::steady_clock::now();
		cout << setw(24) << "Eigen built-in "<< setw(15) 
			<< std::chrono::duration_cast<std::chrono::duration<double> >
			(finish - start).count() << setw(20) << (Z.cast<double>() - ans).norm() << "\n\n\n";
	}
}

