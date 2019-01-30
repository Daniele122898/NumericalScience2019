#include<iostream>
#include<assert.h>

using namespace std;

void int_to_bits(int x) {
    for (int i = sizeof(int)*8-1; i>=0; i--) {
        cout << (bool) (x&(1<<i));
    }
    cout << endl;
}

void float_to_bits(float x) {
    assert(sizeof(float) == sizeof(int));
    int* num = (int*)&x;
    int_to_bits(*num);
}

int main() {
    int_to_bits(4);
    float_to_bits(3.5f);
    return 0;
}