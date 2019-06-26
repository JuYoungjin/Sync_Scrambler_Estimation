#include "MakeStream.h"
#include "general.h"
#include <stdlib.h>


// Input  >> DS : Spreading Code, Dslength : length of DS, m : message, mlength : bit length of message
// Output >> T : Stream that is made from message by Spreading
// By using spreading Code DS and message m, Save spreaded stream to S (input S must be filled with 0)
void Spreading(uint64_t *T, uint64_t *DS, int DSlength, uint64_t *m, int mlength)
{
	int i, j;

	for (i = 0; i < mlength; i++)
	{
		for (j = 0; j < DSlength; j++)
		{
			T[(i*DSlength + j) / 64] ^= ((~((m[i / 64] >> (63 - (i % 64))) ^ (DS[j / 64] >> (63 - (j % 64)))) & 1) << (63 - (i*DSlength + j) % 64));
		}
	}
}

// Input  >> DS : Spreading Code, Dslength : length of DS, m : message, mlength : bit length of message, poly : Scrambling polynomial, degree : degree of poly, IV : initial value for Scrmabling
// Output >> S : Stream that is made from message by Spreading and Synch Scrambling
// Save stream S that is made from message by spreading and scrambling
void makeSyncStr(uint64_t *S, uint64_t *DS, int DSlength, uint64_t *m, int mlength, uint64_t poly, int degree, uint64_t IV)
{
	int i, j;
	uint64_t *temp;
	uint64_t tempS;
	uint64_t state, degreemask;
	uint64_t s;

	temp = (uint64_t *)calloc(mlength*DSlength / sizeof(uint64_t) + 1, sizeof(uint64_t));

	Spreading(temp, DS, DSlength, m, mlength);				//Spreading	

	printf("\nStream after Spreaded : \n");
	if (DSlength > 64) {
		for (i = 0; i < 1024 / DSlength; i++)
		{
			printBIT(temp, DSlength*i, 0, DSlength);
		}
	}
	else {
		for (i = 0; i < 256; )
		{
			printBIT(temp, i, DSlength, (64 / DSlength)*DSlength);
			i += (64 / DSlength)*DSlength;
		}
	}

	degreemask = 0;
	for (i = 0; i < degree; i++)
	{
		degreemask ^= ((uint64_t)1 << i);
	}
	state = IV & degreemask;									//Save the IV value in state

	printf("\nScrambler information"); 
	printf("\n - Polynomial : "); printpoly(poly, degree);
	printf(" - Initial Value for Scrambling : %llx \n", state);

	for (i = 0; i < mlength*DSlength; i++)
	{
		s = 0;
		tempS = state & (poly >> 1);

		for (j = 0; j < degree; j++)
		{
			s ^= (tempS >> j) & 1;							// Calculate LFSR value 
		}
		state = ((state << 1)&degreemask) ^ (s & 1);
		s ^= (temp[i / 64] >> (63 - (i % 64)) & 1);
		S[i / 64] ^= (s << (63 - (i % 64)));
	}

	printf("\nStream after Scrambled S : \n");
	printBIT(S, 0, 0, 512);

	free(temp);
}

