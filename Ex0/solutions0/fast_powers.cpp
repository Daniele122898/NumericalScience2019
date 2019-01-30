#include<bits/stdc++.h>

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

double power_fast(double a, int b) {
	if (b < 0) {
		a = 1/a;
		b = -b;
	}
	if (b == 0)
		return 1;
	if (b == 1)
		return a;
	if (b & 1)
		return a * power_fast(a, b-1);
	return power_fast(a*a, b/2);
}

double power_fast(double a, double b) {
	return exp(b * log(a));
}

	

int main() {
	double err1 = 0, err2 = 0, err3 = 0;
	int count = 0;
	for (double a = 1; a < 100; a++) {
		for (int b = 0; b < 300; b++) {
			if(isinf(power(a,b)) || isinf(power_fast(a,b)) ||
				isinf(power_fast(a, (double) b)) || isinf(pow(a,b)))
				continue;
			err1 += abs((power(a,b) - pow(a,b)) / pow(a,b));
			err2 += abs((power_fast(a,b) - pow(a,b)) / pow(a,b));
			err3 += abs((power_fast(a, (double) b) - pow(a,b)) / pow(a,b));
			count++;
		}
	}
	cout << err1/count << " " << err2/count  << " " << err3/count;
}
