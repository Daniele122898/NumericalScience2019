#include<iostream>
#include<cmath>

using namespace std;

double power(double a, int b) {
    double tmp = 1;
	if (b < 0) {
		a = 1/a;
		b = -b;
	}
	while (b > 0) {
		tmp *= a;
		b--;
	}
	return tmp;
}

double fast_power(double a, int b) {
	if (b < 0) {
		a = 1/a;
		b = -b;
	}
	if (b == 0)
		return 1;
	if (b == 1)
		return a;
	if (b&1) // this means its odd
		return a * fast_power(a, b-1);
	return fast_power(a*a, b/2);
}

double fast_power(double a, double b) {
	return exp(b * log(a));
}

int main() {
    cout << power(2,8) << endl;
    cout << fast_power(2.,8) << endl;
    cout << fast_power(2.,8.) << endl;
    return 0;
}