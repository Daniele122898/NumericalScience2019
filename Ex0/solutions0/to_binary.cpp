#include<bits/stdc++.h>

using namespace std;

void int_to_bits(int x) {
	for (int i = sizeof(int) * 8 - 1; i >= 0; i--)
		cout << (bool) (x & (1 << i));
}

void float_to_bits(float x) {
	assert(sizeof(int) == sizeof(float)); // to be sure that sizes actually match
	int* d = (int*) &x;
	int_to_bits(*d);
}

int main() {
	float x = -7943156.899;
	float_to_bits(x);

}
