#include <stdio.h>

int main() {
	float A[10]; // This is an array of 10 floats
	int N = 10;

	// Set values for A
	for (int i = 0; i < N; i++) {
		// Set value for A[i]
		A[i] = (float)i;
	}

	// Now print each value
	for (int i = 0; i < N; i++) {
		printf("A[%d] = %g\n", i, A[i]);
	}

	return 0;
}