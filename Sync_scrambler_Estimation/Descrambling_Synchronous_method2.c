#include "Descrambling_Synchronous_Method2.h"

// Input  >> S : Stream after spreading and scrambling, streamlength : bit length of stream S
// Output >> Primitive polynomial for scrambling(stored in findpoly), Length of Spreading code(stored in DSlength), and descrambled stream S
// Descramble S after estimating primitive polynomial and initial value.
int descramblingSync_Method2(uint64_t *S, int streamlength, int DSlength, int WALSHorBARKERorPN)
{
	uint64_t poly;
	uint64_t IV;
	int deg;

	if (WALSHorBARKERorPN == WALSH)
		poly = FindSyncScramblerPolynomial_Method2_WALSH(S, streamlength, &deg);
	else if (WALSHorBARKERorPN == BARKER)
		poly = FindSyncScramblerPolynomial_Method2_BARKER(S, streamlength, DSlength, &deg);
	else if (WALSHorBARKERorPN == PN)
		poly = FindSyncScramblerPolynomial_Method2_PN(S, streamlength, DSlength, &deg);

	if (poly != FAIL) printf("\nEstimated primitive polynomial : "); printpoly(poly, deg);

	if (poly == FAIL) return FAIL;

	IV = FindInitialvlaue(S, poly, DSlength, &deg);
	printf("\nEstimated IV : %llx\n", IV);

	descramblingSync_poly(S, streamlength, poly, deg, IV);
	
	return poly;
}

// Input  >> S : Stream after spreading and scrambling, streamlength : bit length of stream S
// Output >> Primitive polynomial for scrambling & degree of primitive polynomial
// Estimate primitive polynomial and degree of primitive polynomial by using Method2 from S that is scrambled and spreaded by WALSH code.
uint64_t FindSyncScramblerPolynomial_Method2_WALSH(uint64_t *S, int streamlength, int *deg)
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

	// calculate minimal polynomial of (U_0,U_2,...)		
	for (i = 0; i < 3 * MAXDEGREE; i++)
	{
		InsertBIT(temp, BIT(U, 2 * i), i);
	}
	poly[0] = BM(temp, 3 * MAXDEGREE, &deg_temp[0]);

	// calculate minimal polynomial of (U_1,U_3,...) 
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

		// calculate primitive polynomial of (U_even number,...)
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
// Estimate primitive polynomial and degree of primitive polynomial by using Method2 from length of spreading code and S that is scrambled and spreaded by BARKER code.
uint64_t FindSyncScramblerPolynomial_Method2_BARKER(uint64_t *S, int streamlength, int DSlength, int *deg)
{
	int i, j, cnt = 0;
	int index;
	int tempdeg;
	uint64_t *U;											//U_i = S_i + S_(i+m) + S_(i+1) + S_(i+m+1) 
	uint64_t Z[100] = { 0 };
	uint64_t Ztemp[100] = { 0 };
	uint64_t b[2] = { 0 };
	uint64_t mask;
	uint64_t poly;
	uint64_t *temp_poly;									
	uint64_t temp, s;
	uint64_t temp2[2];

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

	temp_poly = (uint64_t *)calloc(DSlength, sizeof(uint64_t));			
	for (i = 0; i < DSlength; i++)
	{
		temp2[0] = 0;
		temp2[1] = 0;
		for (j = 0; j < 128; j++)
		{
			InsertBIT(temp2, BIT(U, DSlength * j + i), j);	 // Decimation for each DSlength bit of U 
		}
		temp_poly[i] = BM(temp2, 128, &tempdeg);
	}

	// Find a part that contains message information
	for (i = 0; i < DSlength; i++)
	{
		if (temp_poly[i] == temp_poly[(i + 1) % DSlength]) {
			index = i; break; 
		}
	}
	if (i == DSlength) {
		free(temp_poly);
		free(U);
		return FAIL;
	}
	cnt = 0;
	for (i = 0; i < DSlength; i++)
	{
		if (temp_poly[index] == temp_poly[i]) {
			cnt++;
		}
		else {
			randomIndex = i;
		}
	}
	if (cnt != DSlength - 1) { 
		free(temp_poly);
		free(U);
		return FAIL; 
	}

	// Save (randomIndex - 2 mod DSlength) to plainIndex => By using linear equation, restore randomIndex part
	plainIndex = (randomIndex - 2 + DSlength) % DSlength;		
	
	// Restore PN sequence by using linear equation   
	// Z[i] : Value in plainIndex  
	// b : Value of plainIndex + 1 
	for (i = 0; i < 63; i++)
	{
		Z[0] ^= BIT(U, (i*DSlength + plainIndex)) << (63 - i);
	}
	b[0] = BIT(U, plainIndex + 1) << 63;
	for (i = 1; i < 100; i++)
	{
		Z[i] = (Z[i - 1] << 1) ^ (BIT(U, ((i + 62)*DSlength + plainIndex)) << 1);
		b[i / 64] ^= BIT(U, (i*DSlength + plainIndex + 1)) << (63 - i % 64);
	}
	
	for (i = 0; i < 100; i++)
	{
		Ztemp[i] = Z[i];
	}
	for (i = 0; i < 100; i++)
	{
		Ztemp[i] = Ztemp[i] ^ ((b[i / 64] >> (63 - i % 64)) & 1);
	}
	//	 print Matrix (Z|b) 
	printf("\n\n");
	for (i = 0; i < 9; i ++ )
	{
		printf("%d     - ", i + 1);
		printBIT(Ztemp, 64 * i, 63, 64);
	}
	for (i = 9; i < 99; i++)
	{
		printf("%d    - ", i + 1);
		printBIT(Ztemp, 64 * i, 63, 64);
	}
	printf("%d   - ", 100);
	printBIT(Ztemp, 64 * 99, 63, 64);
	printf("\n\n");
	/////////////////////////////

	mask = GE2(Z, b, 100, 63, &tempdeg);		

	//	 print Matrix (Z|b) after GE
	printf("\n\n");
	for (i = 0; i < 9; i++)
	{
		printf("%d     - ", i + 1);
		printBIT(Z, 64 * i, 63, 64);
	}
	for (i = 9; i < 99; i++)
	{
		printf("%d    - ", i + 1);
		printBIT(Z, 64 * i, 63, 64);
	}
	printf("%d   - ", 100);
	printBIT(Z, 64 * 99, 63, 64);
	printf("\n\n");
	/////////////////////////////



	if (tempdeg > MAXDEGREE) {
		free(temp_poly);
		free(U);
		return FAIL;
	}
	for (i = 0; i < 3 * MAXDEGREE; i += DSlength)
	{
		temp = 0;
		s = 0;
		for (j = 0; j < tempdeg; j++)
		{
			temp ^= (U[(i + j * DSlength + randomIndex - 1) / 64] >> (63 - (i + j * DSlength + randomIndex - 1) % 64) & 1) << j;	
		}
		temp = temp & mask;

		for (j = 0; j < tempdeg; j++)
		{
			s ^= (temp >> j) & 1;						// s is restored bit of randomIndex	part		
		}

		InsertBIT(U, s, (uint64_t)(i + randomIndex));
	}

	printf("\nRestored Scrambling Sequence\n"); 
	for (i = 0; i < 256; )
	{
		printBIT(U, i, DSlength, (64 / DSlength)*DSlength);
		i += (64 / DSlength)*DSlength;
	}
	   
	poly = BM(U, 3 * MAXDEGREE, &tempdeg);								
	*deg = tempdeg;
	free(temp_poly);
	free(U);

	return poly;
}

// Input  >> S : Stream after spreading and scrambling, streamlength : bit length of stream S, DSlength : length of spreading code(estimated)
// Output >> Primitive polynomial for scrambling
// Estimate primitive polynomial and degree of primitive polynomial by using Method2 from length of spreading code S that is scrambled and spreaded by PN code.
uint64_t FindSyncScramblerPolynomial_Method2_PN(uint64_t *S, int streamlength, int DSlength, int *deg)
{
	int i, j, cnt = 0;
	int index, temp_degree;
	int tempdeg;
	uint64_t *U;											//U_i = S_i + S_(i+m) + S_(i+1) + S_(i+m+1) 
	uint64_t Z[100] = { 0 };
	uint64_t Ztemp[100] = { 0 };
	uint64_t b[2] = { 0 };
	uint64_t mask;
	uint64_t poly;
	uint64_t *temp_poly;									
	uint64_t temp, s;
	uint64_t temp2[2];
	int *check_deg;

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
	for (i = 0; i < 256; )
	{
		printBIT(U, i, DSlength, (64 / DSlength)*DSlength);
		i += (64 / DSlength)*DSlength;
	}

	temp_poly = (uint64_t *)calloc(DSlength, sizeof(uint64_t));		
	check_deg = (int *)calloc(DSlength, sizeof(int));	

	for (i = 0; i < DSlength; i++)
	{
		temp2[0] = 0;
		temp2[1] = 0;
		for (j = 0; j < 128; j++)
		{
			InsertBIT(temp2, BIT(U, DSlength * j + i), j);	// Decimation for each DSlength bit of U 		
		}
		temp_poly[i] = BM(temp2, 128, &temp_degree);
		check_deg[i] = temp_degree;
	}

	// Find a part that contains message information
	for (i = 0; i < DSlength; i++)
	{
		if (check_deg[i] == 0 || check_deg[i] == 1){
			cnt++;
		}
	}
	if (cnt == DSlength - 1 || cnt == DSlength) {
		free(temp_poly);
		free(U);
		free(check_deg);
		return FAIL;	
	}
	cnt = 0;	
	for (i = 0; i < DSlength; i++)
	{
		if (temp_poly[i] == temp_poly[(i + 1) % DSlength]) {
			index = i; break;
		}
	}
	if (i == DSlength) {
		free(temp_poly);
		free(U);
		free(check_deg);
		return FAIL;
	}
	for (i = 0; i < DSlength; i++)
	{
		if (temp_poly[index] == temp_poly[i]) cnt++;
		else randomIndex = i;
	}
	if (cnt != DSlength - 1) {
		free(temp_poly);
		free(U);
		free(check_deg);
		return FAIL;
	}
	// Save (randomIndex - 2 mod DSlength) to plainIndex => By using linear equation, restore randomIndex part
	plainIndex = (randomIndex - 2 + DSlength) % DSlength;		

	// Restore PN sequence by using linear equation   
	// Z[i] : Value in plainIndex 
	// b : Value of plainIndex + 1 
	for (i = 0; i < 63; i++)
	{
		Z[0] ^= BIT(U, (i*DSlength + plainIndex)) << (63 - i);
	}
	b[0] = BIT(U, plainIndex + 1) << 63;
	for (i = 1; i < 100; i++)
	{
		Z[i] = (Z[i - 1] << 1) ^ (BIT(U, ((i + 62) * DSlength + plainIndex)) << 1);
		b[i / 64] ^= BIT(U, (i*DSlength + plainIndex + 1)) << (63 - i % 64);
	}
		
	for (i = 0; i < 100; i++)
	{
		Ztemp[i] = Z[i];
	}
	for (i = 0; i < 100; i++)
	{
		Ztemp[i] = Ztemp[i] ^ ((b[i / 64] >> (63 - i % 64)) & 1);
	}
	//	 print Matrix (Z|b) 
	printf("\n\n");
	for (i = 0; i < 9; i++)
	{
		printf("%d     - ", i + 1);
		printBIT(Ztemp, 64 * i, 63, 64);
	}
	for (i = 9; i < 99; i++)
	{
		printf("%d    - ", i + 1);
		printBIT(Ztemp, 64 * i, 63, 64);
	}
	printf("%d   - ", 100);
	printBIT(Ztemp, 64 * 99, 63, 64);
	printf("\n\n");
	/////////////////////////////

	mask = GE2(Z, b, 100, 63, &tempdeg);	

	//	 print Matrix (Z|b) after GE
	printf("\n\n");
	for (i = 0; i < 9; i++)
	{
		printf("%d     - ", i + 1);
		printBIT(Z, 64 * i, 63, 64);
	}
	for (i = 9; i < 99; i++)
	{
		printf("%d    - ", i + 1);
		printBIT(Z, 64 * i, 63, 64);
	}
	printf("%d   - ", 100);
	printBIT(Z, 64 * 99, 63, 64);
	printf("\n\n");
	/////////////////////////////

	if (tempdeg > MAXDEGREE) {
		free(temp_poly);
		free(U);
		free(check_deg);
		return FAIL;
	}
	for (i = 0; i < 3 * MAXDEGREE; i += DSlength)
	{
		temp = 0;
		s = 0;
		for (j = 0; j < tempdeg; j++)
		{
			temp ^= (U[(i + j * DSlength + randomIndex - 1) / 64] >> (63 - (i + j * DSlength + randomIndex - 1) % 64) & 1) << j;	
		}
		temp = temp & mask;

		for (j = 0; j < tempdeg; j++)
		{
			s ^= (temp >> j) & 1;						// s is restored bit of randomIndex	part		
		}

		InsertBIT(U, s, (i + randomIndex)); 
	}

	printf("\nRestored Scrambling Sequence\n"); 
	for (i = 0; i < 256; )
	{
		printBIT(U, i, DSlength, (64 / DSlength)*DSlength);
		i += (64 / DSlength)*DSlength;
	}

	poly = BM(U, 3 * MAXDEGREE, &temp_degree);						
	free(temp_poly);
	free(U);
	free(check_deg);

	return poly;
}

