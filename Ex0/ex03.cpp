#include<bits/stdc++.h>


template<class T> T factorial(T x) {
    T temp = x;
    while(x > 1) {
        temp *= --x;
    }
    return temp;
}

template<class T> double dbl_factorial(T x) {
    return factorial((double) std::round(x));
}

template<class T> void test(T max) {
	// To find the maximum correct value for the factorial, consider the
	// maximum value of the type considered and divide it by
	// 2, 3, 4, ... until you get something smaller than 1. 
	// Why this works? Hint: if a,b,c are integral types, then (a/b)/c == a/(b*c).
	T i = 1, x = max;
	while (x >= 1) {
		x /= ++i;
	}	
	std::cout << "Factorial of " << i-1 << " is " << factorial(i-1)
		 << ", but factorial of " << i << " is not " << factorial(i) << std::endl;
}

int main() {
    std::cout << factorial<int>(3) << std::endl;
    std::cout << factorial<int>(4) << std::endl;
    std::cout << factorial<long>(10) << std::endl;
	test(INT_MAX); test(LONG_MAX); test(DBL_MAX);

}