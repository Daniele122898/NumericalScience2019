#include<bits/stdc++.h>

using namespace std;

struct vec {
	int capacity = 0;
	int size = 0;
	double* data = nullptr;

	vec() {
		capacity = 10;
		data = new double[capacity];
	}

	~vec() {
		delete[] data;
	}

	void push_back(double x) {
		if (size >= capacity) {
			capacity *= 2;
			double* newdata = new double[capacity];
			for (int i=0; i < size; i++) 
				newdata[i] = data[i];
			delete[] data;
			data = newdata;
		}
		data[size++] = x;
	}
};


int main() {
	vec v;
	for (int i=0; i < 25; i++){
		v.push_back(i);
		cout << v.data[i] << " ";
	}

	cout << endl << v.capacity;
	
}
	
