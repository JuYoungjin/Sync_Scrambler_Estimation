#include "general.h"

// Input  >> A : matrix, b : vector, rowdim : row dimension of A, coldim : column dimension of A
// Output >> x such that Ax = b
// By using Gaussian Elimination, find x that satisfies the condition Ax=b
uint64_t GE2(uint64_t* A, uint64_t *b, int rowdim, int coldim, int *deg)
{
	int i, j, rank = coldim;
	uint64_t temp;
	uint64_t res = 0;

	for (i = 0; i < rowdim; i++)
	{
		A[i] = A[i] ^ ((b[i / 64] >> (63 - i % 64)) & 1); 
	}

	for (i = 0; i < rowdim; i++)
	{
		for (j = i; j < rowdim; j++) 
		{
			if (((A[j] >> (63 - i)) & 1) != 0)
			{
				if (j != i) {
					temp = A[j];
					A[j] = A[i];
					A[i] = temp;
				}
				break;
			}
		}

		if (j == rowdim) 
		{
			rank = i;
			break;
		}

		for (j = 0; j < rowdim; j++)
		{
			if ((j != i) && (((A[j] >> (63 - i)) & 1) != 0)) 
			{
				A[j] = A[j] ^ A[i];
			}
		}
	}
	for (i = 0; i < rank; i++) 
	{
		res = res ^ ((A[i] & 1) << i);
	}
	*deg = rank;

	return res;
}

// Input  >> str : 64bit array, length : length of array 
// output >> minimal polynomial & degree of minimal polynomial
// Berlekamp-Massey Algorithm
uint64_t BM(uint64_t *str, int length, int *deg)
{
	uint64_t C, B;                                    
	uint64_t T = 0;
	uint64_t temp;
	int L, m, N, d;                                  
	int i;

	C = 1; B = 1;
	L = 0; m = -1; N = 0;                            

	while (N < length)
	{
		d = 0;
		temp = 0;
		if ((N) / 64 - (N - L) / 64 > 0)          
		{
			temp ^= (str[N / 64] >> (64 - ((N + 1) % 64)));
			temp ^= (str[(N / 64) - 1] << (((N + 1) % 64)));
		}
		else
		{
			temp ^= (str[N / 64] >> (64 - ((N + 1) % 64)));
		}
		temp = temp & C;                      
		for (i = 0; i <= L; i++)
		{
			d ^= (temp >> i) & 1;        
		}

		if (d == 1)
		{
			T = C;
			C = C ^ (B << (N - m));      
			if (L <= N / 2)
			{
				L = N + 1 - L;
				m = N;
				B = T;
			}
		}
		N++;
	}
	*deg = L;
	return C;
}

// gcd of a & b
uint64_t gcd(uint64_t a, uint64_t b)
{
	uint64_t tmp, n;

	if (a < b) {
		tmp = a;
		a = b;
		b = tmp;
	}

	while (b != 0)
	{
		n = a % b;
		a = b;
		b = n;
	}

	return a;
}

// Print Polynomial(deg is degree of polynomial)
void printpoly(uint64_t poly, int deg)
{
	int i;
	for (i = 0; i < deg; i++)
	{
		if (((poly >> i) & 1) == 1) printf("x^%d + ", i);
	}
	if (((poly >> deg) & 1) == 1) printf("x^%d", i);
	printf("\n");
}

// print Bit array as 0,1 representation
void printBIT(uint64_t *S, int start, int space, int length)
{
	int i;
	if (space == 0) space = INT32_MAX;
	for (i = 0; i < length; i++)
	{
		if (i%space == 0 && i != 0) printf(" . ");
		printf("%lld", BIT(S, start + i));
	}
	printf("\n");
}