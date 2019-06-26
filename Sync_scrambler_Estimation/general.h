#ifndef _GENERAL_H_
#define _GENERAL_H_

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>

#define SUCCESS 0
#define FAIL -1
#define SELFSYNCHFAIL 1

#define WALSH 1
#define BARKER 2
#define PN 3

#define SYNCH 1
#define SELFSYNCH 2

#define MAXDEGREE 50
#define WRONGSPREADINGCODELENGTH 0xFFFFFFFFFFFFFFFF
#define NOTSELFSYNCHRONOUS 0xFFFFFFFFFFFFFFFE

#define BIT(S,n) ((S[(n)/64]>>(63-(n)%64))&1)
#define InsertBIT(S,x,n) do { uint64_t TEMP=(x); TEMP=~TEMP&1; S[(n)/64]=((S[(n)/64]|((uint64_t)1<<(63-((uint64_t)n)%64)))^(TEMP<<(63-((uint64_t)n)%64))); }while(0)

// Input  >> A : matrix, b : vector, rowdim : row dimension of A, coldim : column dimension of A
// Output >> x such that Ax = b
// By using Gaussian Elimination, find x that satisfies the condition Ax=b
uint64_t GE2(uint64_t* A, uint64_t *b, int rowdim, int coldim, int *deg);


// Input  >> str : 64bit array, length : length of array 
// Output >> minimal polynomial & degree of minimal polynomial
// Calculate minimal polynomial which generates bit string str and save the degree of the polynomial in deg
// Berlekamp-Massey Algorithm
uint64_t BM(uint64_t *str, int length, int *deg);

// gcd of a & b
uint64_t gcd(uint64_t a, uint64_t b);

// Print Polynomial(deg is degree of polynomial)
void printpoly(uint64_t poly, int deg);

// print Bit array as 0,1 representation
void printBIT(uint64_t *S, int start, int space, int length);

#endif