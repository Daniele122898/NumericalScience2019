#include <iostream>
#include <algorithm>
#include <vector>
#include <sstream>
#include <string>
#include <complex>

using namespace std;

const complex<double> I(0,1); // imaginary unit
const double PI = 3.14159265359;

template<class T> struct duplet {
	int ind;
	T val;

	duplet(int p, T v) {
		ind = p;
		val = v;
	}
};

template<class T> struct sparse_vec {
	double tol = 1e-6;
	vector<duplet<T>> duplets;
	int len=0;

	sparse_vec(int l) {
		len = l;
	}

	void append(int ind, T val) {
		if (abs(val)<tol)
		 return;
		duplet<T> newduplet(ind, val);
		duplets.push_back(newduplet);
	}

	void cleanup2() {
		std::sort(duplets.begin(), duplets.end(), 
			[] (duplet<T> x, duplet<T> y) { return x.ind < y.ind; });

		vector<duplet<T> > newduplets;
		T tmp = 0;
		for (int i=0; i<duplets.size(); i++) {
			if (i == duplets.size()-1 || duplets[i+1].ind != duplets[i].ind) {
				tmp += duplets[i].val;
				if (abs(tmp) >= tol) {
					duplet<T> newduplet(duplets[i].ind, tmp);
					newduplets.push_back(newduplet);
				}
				tmp = 0;
			}
			else {
				tmp += duplets[i].val;
			}
		}
		duplets = newduplets;
	}
			
	void cleanup() {
		std::sort(duplets.begin(), duplets.end(), 
			[] (const duplet<T>& d1, const duplet<T>& d2) { return d1.ind < d2.ind; });
		for (int i=0; i<duplets.size(); i++) {
			// check if we're at the point where the index is above the length
			if(duplets[i].ind >= len) {
				// this means we need to remove this and all subsequent elemts
				duplets.erase(duplets.begin()+i, duplets.end());
				break;
			}
			// lets check if the next element has the same index
			int next = (i+1);
			if (next < duplets.size() &&  duplets[i].ind == duplets[next].ind) {
				// the indices are equal. so add them up and remove the later one
				duplets[i].val += duplets[next].val;
				duplets.erase(duplets.begin()+next);
				// if the next one is equal too we need to restart the loop otherwise errors will occur.
				if (next < duplets.size() &&  duplets[i].ind == duplets[next].ind) {
					i--;
					continue;
				}
			}
			// this should be the final duplet with that index.
			// so now lets check if its value is in the tolerated zone
			if(abs(duplets[i].val)<tol) {
				// its value is too small so lets remove it
				duplets.erase(duplets.begin()+i);
				i--;
				continue;
			}

		}
	}

	T get_val(int ind) const {
		// get value by index using binary search
		if (duplets.empty())
			return 0;
		return _get_val(ind, 0, duplets.size());
	}

	T _get_val(int ind, int n1, int n2) const {
		if (n2 < n1)
			return 0;
		int m = (n1+n2)/2;
		if (duplets[m].ind == ind)
			return duplets[m].val;
		if (duplets[m].ind < ind)
			return _get_val(ind, m+1, n2);
		if (duplets[m].ind > ind)
			return _get_val(ind, n1, m-1);
	}

	static sparse_vec cwise_mult(const sparse_vec &a, const sparse_vec &b) {
		sparse_vec out(max(a.len,b.len));
		for (int i=0; i<out.len; i++) {
			out.append(i, a.get_val(i)*b.get_val(i));
		}
		return out;
	}

	static sparse_vec conv(const sparse_vec &a, const sparse_vec &b) {
		sparse_vec out(a.len + b.len - 1);
		for (auto x : a.duplets) {
			for (auto y : b.duplets) {
				out.append(x.ind + y.ind, x.val * y.val);
			}
		}
		return out;
	}

	static sparse_vec fft(const sparse_vec &x) {
		int n = x.len;
		if (n<= 1) return x;

		sparse_vec even(n/2);
		sparse_vec odd(n/2);
		for (auto duplet : x.duplets) {
			if (duplet.ind % 2 == 0)
				even.append(duplet.ind/2, duplet.val);
			else 
				odd.append((duplet.ind-1)/2, duplet.val);
		}

		sparse_vec f0 = fft(even);
		sparse_vec f1 = fft(odd);
		
		T omega = exp(-2.*PI/n*I);
		T s(1., 0.);
		sparse_vec tot(n);

		for (int k=0; k<n; k++) {
			tot.append(k, f0.get_val(k % (n/2)) + f1.get_val(k % (n/2))*s);
			s*= omega;
		}

		return tot;
	}

	static sparse_vec ifft(const sparse_vec &x) {
		double n = x.len;
		sparse_vec out(n);
		sparse_vec x_conj(n);
		for (auto duplet : x.duplets) {
			x_conj.append(duplet.ind, std::conj(duplet.val));
		}
		for (auto duplet : fft(x_conj).duplets) {
			out.append(duplet.ind, std::conj(duplet.val)/n);
		}

		return out;
	}

	static sparse_vec conv_fft(sparse_vec a, sparse_vec b) {
		int N = a.len + b.len -1;
		a.len = N;
		b.len = N;
		return ifft(cwise_mult(fft(a), fft(b)));
	}

	std::string to_string() const {
		std::stringstream ss;
		for (auto p : this->duplets) {
			ss << "(" << p.ind << "," << p.val << "),";
		}
		ss << "\n";
		std::string out = ss.str();	
		return out;
	}
			

};

void test() {
	cout << "TESTING MY OWN SHIT LUL" << endl;
	sparse_vec<complex<double> > x(5);
	x.append(0,complex<double>(8.2,0));
	x.append(1,complex<double>(1,-2));
	x.append(3,complex<double>(-3,4.66));
	x.append(4,complex<double>(0,4));
	x.append(4,complex<double>(0,4));
	x.append(5,complex<double>(0,4));
	x.append(6,complex<double>(0,4));
	x.cleanup();

	sparse_vec<complex<double> > y(5);
	y.append(0,complex<double>(8.2,0));
	y.append(1,complex<double>(1,-2));
	y.append(3,complex<double>(-3,4.66));
	y.append(4,complex<double>(0,4));
	y.append(4,complex<double>(0,4));
	y.append(5,complex<double>(0,4));
	y.append(6,complex<double>(0,4));
	y.cleanup2();

	cout << "TESTING CLEAN UP" << endl;
	cout << "x = " << x.to_string() << endl;
	cout << "y = " << y.to_string() << endl;

	cout << "\nDONE TESTING MY OWN SHIT\n" << endl; 
}




/***** TESTING ******/

int main() {

	test();

	sparse_vec<complex<double> > x(5);
	x.append(0,complex<double>(8.2,0));
	x.append(1,complex<double>(1,-2));
	x.append(3,complex<double>(-3,4.66));
	x.append(4,complex<double>(0,4));
	x.cleanup();

	sparse_vec<complex<double> > y(4);
	y.append(1,complex<double>(5,0));
	y.append(2,complex<double>(1.21,-4));
	y.append(3,complex<double>(4,2.4));
	y.cleanup();

	auto m = sparse_vec<complex<double> >::cwise_mult(x,y);
	m.cleanup();
	cout << "TESTS. Correct componentwise multiplication between x and y: (1,(5,-10)),(3,(-23.184,11.44)),\n";
	cout << "cwise_mult(x,y) = " << m.to_string();

	auto c = sparse_vec<complex<double> >::conv(x,y);
	c.cleanup();
	cout << "Correct exact discrete convolution between x and y: (1,(41,0)),(2,(14.922,-42.8)),(3,(26.01,13.26)),(4,(-6.2,17.7)),(5,(15.01,37.6386)),(6,(-7.184,16.28)),(7,(-9.6,16)),\n";
	cout << "conv(x,y) = " << c.to_string();
	auto cf = sparse_vec<complex<double> >::conv_fft(x,y);
	cf.cleanup();
	cout << "conv_fft(x,y) = " << cf.to_string();
}


