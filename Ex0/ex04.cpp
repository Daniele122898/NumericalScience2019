#include<bits/stdc++.h>

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
            double* new_data = new double[capacity];
            for (int i=0; i<size; i++)
                new_data[i] = data[i];
            delete[] data;
            data = new_data;
        }
        data[size++] = x;        
    }
};

int main() {
    vec v;
	for (int i=0; i < 25; i++){
		v.push_back(i);
		std::cout << v.data[i] << " ";
	}

	std::cout << std::endl << v.capacity << std::endl;
    return 0;
}