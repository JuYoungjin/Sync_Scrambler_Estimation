#include "Descrambling_Synchronous_Method1.h"

// Input  >> S : Stream after spreading and scrambling, streamlength : bit length of stream S, DSlength : length of Spreading code, WALSHorBARKERorPN : Config that indicate Spreading Code
// Output >> Primitive polynomial for scrambling(stored in findpoly), Length of Spreading code(stored in DSlength), and descrambled stream S
// Descramble S after estimating primitive polynomial and initial value.
int DescramblingSync_Method1(uint64_t *S, int streamlength, int DSlength, int WALSHorBARKERorPN)
{
	uint64_t poly, IV;
	int deg;

	if (WALSHorBARKERorPN == WALSH)
		poly = FindSyncScramblerPolynomial_Method1_WALSH(S, streamlength, &deg);
	else if (WALSHorBARKERorPN == BARKER)
		poly = FindSyncScramblerPolynomial_Method1_BARKER(S, streamlength, DSlength, &deg);
	else if (WALSHorBARKERorPN == PN)
		poly = FindSyncScramblerPolynomial_Method1_PN(S, streamlength, DSlength, &deg);

	if (poly != FAIL) printf("\nEstimated primitive polynomial : "); printpoly(poly, deg);

	if (poly == FAIL) return FAIL;

	IV = FindInitialvlaue(S, poly, DSlength, &deg);
	printf("\nEstimated IV : %llx\n", IV);

	descramblingSync_poly(S, streamlength, poly, deg, IV);

	return poly;
}

// Input  >> S : Stream after spreading and scrambling, streamlength : bit length of stream S
// Output >> Primitive polynomial for scrambling & degree of primitive polynomial
// Estimate primitive polynomial and degree of primitive polynomial by using Method1 from S that is scrambled and spreaded by WALSH code.
uint64_t FindSyncScramblerPolynomial_Method1_WALSH(uint64_t *S, int streamlength, int *deg)
{
	uint64_t *U;												//U_i = S_i + S_(i+2) + S_(i+1) + S_(i+3) 
	uint64_t temp[3 * MAXDEGREE / 64 + 1] = { 0 };				//decimation sequence
	uint64_t poly[2];											
	uint64_t polytemp[2];
	int i, j;
	int deg_temp[4];

	U = (uint64_t *)calloc(streamlength / 64 + 1, sizeof(uint64_t));
	for (i = 0; i < streamlength / 64 + 1; i++)
	{
		U[i] = S[i] ^ ((S[i] << 2) ^ (S[i + 1] >> (64 - 2)));
	}
	printf("\nU = S + (S<<2) : \n");
	printBIT(U, 0, 0, 512);

	for (i = 0; i < streamlength / 64 + 1; i++)
	{
		U[i] = U[i] ^ ((U[i] << 1) ^ (U[i + 1] >> (63)));
	}
	printf("\nDelete WALSH & Message information T = U + (U<<1) : \n");
	for (i = 0; i < 256; )
	{
		printBIT(U, i, 2, 32);
		i += 32;
	}

	// calculate primitive polynomial of (U_0,U_2,...)		
	for (i = 0; i < 3 * MAXDEGREE; i++)
	{
		InsertBIT(temp, BIT(U, 2 * i), i); 
	}
	poly[0] = BM(temp, 3 * MAXDEGREE, &deg_temp[0]);			

	// calculate primitive polynomial of (U_1,U_3,...) 
	for (i = 0; i < 3 * MAXDEGREE; i++)
	{
		InsertBIT(temp, BIT(U, 2 * i + 1), i);
	}
	poly[1] = BM(temp, 3 * MAXDEGREE, &deg_temp[1]);


	// Check whether polynomial is right polynomial
	for (j = 20;; j += 20)
	{
		for (i = 0; i < 3 * MAXDEGREE / 64 + 1; i++)
		{
			temp[i] = 0;
		}

		// calculate minimal polynomial of (U_even number,...)
		for (i = 0; i < 3 * MAXDEGREE; i++)
		{
			InsertBIT(temp, BIT(U, 2 * i + j), i);
		}
		polytemp[0] = BM(temp, 3 * MAXDEGREE, &deg_temp[2]);


		for (i = 0; i < 3 * MAXDEGREE / 64 + 1; i++)
		{
			temp[i] = 0;
		}

		// calculate primitive polynomial of (U_odd number,...)
		for (i = 0; i < 3 * MAXDEGREE; i++)
		{
			InsertBIT(temp, BIT(U, 2 * i + j + 1), i);
		}
		polytemp[1] = BM(temp, 3 * MAXDEGREE, &deg_temp[3]);

		if (polytemp[0] == poly[0] && polytemp[1] == poly[1])
			continue;
		else if (polytemp[0] == poly[0] && polytemp[1] != poly[1]) 
		{
			free(U); 
			*deg = deg_temp[0];
			return poly[0];
		}
		else if (polytemp[0] != poly[0] && polytemp[1] == poly[1])
		{
			free(U); 
			*deg = deg_temp[1];
			return poly[1];
		}
		else
		{
			free(U); return FAIL;
		}
	}
}

// Input  >> S : Stream after spreading and scrambling, streamlength : bit length of stream S, DSlength : length of spreading code
// Output >> Primitive polynomial for scrambling & degree of primitive polynomial
// Estimate primitive polynomial and degree of primitive polynomial by using Method1 from length of spreading code and S that is scrambled and spreaded by BARKER code.
uint64_t FindSyncScramblerPolynomial_Method1_BARKER(uint64_t *S, int streamlength, int DSlength, int *deg)
{
	int xTable[7][16] = { { 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1},
		{ 11 , 0 , 43 , 0 , 171 , 0 , 683 , 0 , 2731 , 0 , 10923 , 0 , 43691 , 0 , 174763 , 0},
		{ 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1},
		{ 7 , 13 , 27 , 0 , 103 , 205 , 411 , 0 , 1639 , 3277 , 6555 , 0 , 26215 , 52429 , 104859 , 0},
		{ 5 , 0 , 55 , 37 , 0 , 439 , 293 , 0 , 3511 , 2341 , 0 , 28087 , 18725 , 0 , 224695 , 149797},
		{ 3 , 23 , 13 , 29 , 93 , 0 , 189 , 373 , 745 , 1501 , 2979 , 5981 , 11917 , 23837 , 47709 , 0},
		{ 3 , 5 , 11 , 59 , 59 , 79 , 315 , 0 , 635 , 1339 , 2523 , 5051 , 10083 , 20165 , 40331 , 80699} }; // Table of x such that mx = 2^i 
	int m[7] = { 2,3,4,5,7,11,13 };
	uint64_t *U;												//U_i = S_i + S_(i+m) + S_(i+1) + S_(i+m+1) 
	uint64_t temp[3 * MAXDEGREE / 64 + 1] = { 0 };				//decimation sequence
	uint64_t temp2[3 * MAXDEGREE / 64 + 1] = { 0 };				
	uint64_t decimationPoly, poly;								
	uint64_t state = 1;											//state for LFSR
	uint64_t s, tempS, degreemask;
	int i, j, k, x, n, DSindex;

	for (i = 0; i < 7; i++)
	{
		if (m[i] == DSlength) DSindex = i;
	}

	U = (uint64_t *)calloc(streamlength / 64 + 1, sizeof(uint64_t));
	for (i = 0; i < streamlength / 64 + 1; i++)
	{
		U[i] = S[i] ^ ((S[i] << DSlength) ^ (S[i + 1] >> (64 - DSlength)));
	}
	printf("\nU = S + (S<<m) : \n");
	printBIT(U, 0, 0, 512);

	for (i = 0; i < streamlength / 64 + 1; i++)
	{
		U[i] = U[i] ^ ((U[i] << 1) ^ (U[i + 1] >> (63)));
	}
	printf("\nDelete BARKER & Message information T = U + (U<<1) : \n");
	for (i = 0; i < 256; )
	{
		printBIT(U, i, DSlength, (64 / DSlength)*DSlength);
		i += (64 / DSlength)*DSlength;
	}

	for (i = 0; i < 3 * MAXDEGREE; i++)
	{
		InsertBIT(temp, BIT(U, DSlength * i), i);	// Decimation for each DSlength bit of U  
	}
	decimationPoly = BM(temp, 3 * MAXDEGREE, deg);					// primitive polynomial of (U_0,U_DSlength,...) 

	degreemask = 0; n = *deg;
	for (i = 0; i < n; i++)
	{
		degreemask ^= ((uint64_t)1 << i);
	}
	
	x = xTable[DSindex][n - 5];
	if (x == 0) {
		printf("\nSince gcd of length of spreadinig code %d and PN-sequence period 2^%d -1 (%d) is not 1, descrambling is impossible in this case\n", DSlength, n, (1 << n) - 1);
		return FAIL;
	}

	printf("\ndecimation Polynomial : \n"); printpoly(decimationPoly, n);
	printf("\nLength of Spreading Code m : %d , Period of PN sequence 2^n -1 : %d, gcd(m,2^n -1) : %d \n", DSlength, ((1 << n) - 1), gcd(DSlength, (1 << n) - 1));
	printf("x (that satisfies mx = 2^i (mod 2^n -1)) : %d\n", x);

	// x bits decimation stream from PN sequence of DecimationPoly  
	for (i = 0; i < 3 * MAXDEGREE; i++)
	{
		for (j = 0; j < x; j++)
		{
			s = 0;
			tempS = state & (decimationPoly >> 1);

			for (k = 0; k < n; k++)
			{
				s ^= (tempS >> k) & 1;
			}
			state = ((state << 1)&degreemask) ^ (s & 1);
		}
		InsertBIT(temp2, s, i);
	}

	poly = BM(temp2, 3 * MAXDEGREE, deg);

	return poly;
}

// Input  >> S : Stream after spreading and scrambling, streamlength : bit length of stream S, DSlength : length of spreading code(estimated)
// Output >> Primitive polynomial for scrambling
// Estimate primitive polynomial and degree of primitive polynomial by using Method1 from length of spreading code S that is scrambled and spreaded by PN code.
uint64_t FindSyncScramblerPolynomial_Method1_PN(uint64_t *S, int streamlength, int DSlength, int *deg)
{
	uint32_t xTable[8][16] = { {11, 0, 43, 0, 171, 0, 683, 0, 2731, 0, 10923, 0, 43691, 0, 174763, 0},
	{5, 0, 55, 37, 0, 439, 293, 0, 3511, 2341, 0, 28087, 18725, 0, 224695, 149797},
	{15, 0, 9, 0, 239, 0, 137, 0, 3823, 0, 2185, 0, 61167, 0, 34953, 0},
	{0, 31, 21, 91, 17, 0, 991, 661, 2907, 529, 0, 31711, 21141, 93019, 16913, 0},
	{1, 0, 63, 0, 0, 0, 33, 0, 4031, 0, 0, 0, 2081, 0, 257983, 0},
	{11, 1, 0, 127, 85, 73, 887, 1387, 65, 0, 16255, 10837, 9289, 113527, 177515, 8257},
	{5, 0, 1, 0, 255, 0, 731, 0, 1189, 0, 129, 0, 65279, 0, 187099, 0},
	{15, 0, 43, 1, 0, 511, 341, 0, 273, 7663, 0, 21931, 257, 0, 261631, 174421} }; // Table of x such that mx = 2^i 
	int m[8] = { 3, 7, 15, 31, 63, 127, 255, 511 };

	uint64_t *U;												//U_i = S_i + S_(i+m) + S_(i+1) + S_(i+m+1) 
	uint64_t temp[3 * MAXDEGREE / 64 + 1] = { 0 };				//decimation sequence
	uint64_t temp2[3 * MAXDEGREE / 64 + 1] = { 0 };				
	uint64_t decimationPoly, poly;							
	uint64_t state = 1;											//state for LFSR
	uint64_t s, tempS, degreemask;
	int i, j, k, x, n, DSindex;

	for (i = 0; i < 8; i++)
	{
		if (m[i] == DSlength) DSindex = i;
	}

	U = (uint64_t *)calloc(streamlength / 64 + 1, sizeof(uint64_t));

	for (i = 0; i < streamlength / 64 + 1; i++)
	{
		U[i] = S[i] ^ (S[i + (DSlength / 64)] << (DSlength % 64)) ^ (S[i + (DSlength / 64) + 1] >> (64 - (DSlength % 64)));
	}

	printf("\nU = S + (S<<m) : \n");
	printBIT(U, 0, 0, 512);

	for (i = 0; i < streamlength / 64 + 1; i++)
	{
		U[i] = U[i] ^ ((U[i] << 1) ^ (U[i + 1] >> (63)));
	}
	printf("\nDelete BARKER & Message information T = U + (U<<1) : \n");

	if (DSlength > 64) {
		for (i = 0; i < 1024 / DSlength; i++)
		{
			printBIT(U, DSlength*i, 0, DSlength);
		}
	}
	else {
		for (i = 0; i < 256; )
		{
			printBIT(U, i, DSlength, (64 / DSlength)*DSlength);
			i += (64 / DSlength)*DSlength;
		}
	}

	for (i = 0; i < 3 * MAXDEGREE; i++)
	{
		InsertBIT(temp, BIT(U, DSlength * i), i);	// Decimation for each DSlength bit of U  
	}

	decimationPoly = BM(temp, 3 * MAXDEGREE, deg);	// primitive polynomial of (U_0,U_DSlength,...) 				

	
	if (*deg == 0 || *deg == 1)
	{
		printf("\nSince degree of primitive poly in PN spreading code and degree of primitive poly in PN-sequence used in scrambling are same, descrambling is impossible in this case\n"); 
		return FAIL;
	}

	degreemask = 0; n = *deg;
	for (i = 0; i < n; i++)
	{
		degreemask ^= ((uint64_t)1 << i);
	}

	x = xTable[DSindex][n - 5];

	if (x == 0) {
		printf("\nSince gcd of length of spreadinig code %d and PN-sequence period 2^%d -1 (%d) is not 1, descrambling is impossible in this case\n", DSlength, n, (1 << n) - 1);
		return FAIL;
	}

	printf("\ndecimation Polynomial : \n"); printpoly(decimationPoly, n);
	printf("\nLength of Spreading Code m : %d , Period of PN sequence 2^n -1 : %d, gcd(m,2^n -1) : %d \n", DSlength, ((1 << n) - 1), gcd(DSlength, (1 << n) - 1));
	printf("x (that satisfies mx = 2^i (mod 2^n -1)) : %d\n", x);
	

	// x bits decimation stream from PN sequence of DecimationPoly 
	for (i = 0; i < 3 * MAXDEGREE; i++)
	{
		for (j = 0; j < x; j++)
		{
			s = 0;
			tempS = state & (decimationPoly >> 1);

			for (k = 0; k < n; k++)
			{
				s ^= (tempS >> k) & 1;
			}
			state = ((state << 1)&degreemask) ^ (s & 1);
		}
		InsertBIT(temp2, s, i);
	}
	poly = BM(temp2, 3 * MAXDEGREE, deg);

	return poly;
}

