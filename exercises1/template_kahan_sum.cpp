#include<bits/stdc++.h>

using namespace std;

float naive_sum(vector<float> &v) {
	float sum = 0;
	for(float x : v){
		sum += x;
	}
	return sum;
}

float sum(vector<float> v) {
	float sum = 0;
	sort(v.begin(), v.end(), [] (auto x, auto y) { return abs(x) < abs(y); } );
	for(int i = 0; i < v.size(); i++) {
		sum += v[i];
	}
	return sum;
}

float acc_sum(vector<float> &v) {
	float sum = 0;
	float c = 0;

	for(float x : v){
		float y = x -c;
		float t = sum + y;
		c = (t-sum)-y;
		sum = t;
	}

	return sum;
}

int main() {
	srand(time(0));
	cout << setprecision(15);
	int N = 1e7;
	double corr_sum = 0;
	vector<float> v(N);
	for (int i = 0; i < N; i++) {
		double x = 1e-8 + (rand() % 10) * 1e-9;
		corr_sum += x;
		v[i] = (float) x;
	}
	corr_sum += 1;
	v[1]=1;

	cout << "naive_sum(v) = " << naive_sum(v) << endl
		 << "      sum(v) = " << sum(v) << endl 
		 << "  acc_sum(v) = " << acc_sum(v) << endl
		 << " correct sum = " << corr_sum << endl;
}
