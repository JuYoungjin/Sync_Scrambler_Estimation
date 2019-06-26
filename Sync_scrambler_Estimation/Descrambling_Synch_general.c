#include "Descrambling_Synch_general.h"


// Input  >> S : Stream after spreading and scrambling, poly : primitive polynomial for scramble, DSlength : length of Spreading code
// Output >> Initial Value of Scrambling & degree of primitive polynomial
// Algorithm for finding IV by using primitive polynomial and stream S
uint64_t FindInitialvlaue(uint64_t *S, uint64_t poly, int DSlength, int *deg)
{
	uint64_t *p;
	uint64_t b[2] = { 0 };
	uint64_t Z[100] = { 0 };
	uint64_t temp;
	uint64_t IVtemp;
	uint64_t IV = 0;
	int i, j, size;
	int check_deg = 0;

	for (i = 63; i > 0; i--)
	{
		if ((poly >> i) & 1 != 0)
		{
			check_deg = i;
			break;
		}
	}
	*deg = check_deg;

	size = 100 * DSlength + plainIndex + 2;

	p = (uint64_t *)calloc(size, sizeof(uint64_t));

	poly = poly >> (uint64_t)1;

	p[0] = poly;
	for (i = 1; i < check_deg; i++)
	{
		p[i] = poly >> i;
		for (j = 0; j < i; j++)
		{
			if ((poly >> j) & 1 == 1)
			{
				p[i] ^= p[i - 1 - j];
			}
		}
	}
	for (i = check_deg; i < 100 * DSlength + plainIndex + 2; i++)
	{
		for (j = 0; j < check_deg; j++)
		{
			if ((poly >> (uint64_t)j) & 1 == 1)
			{
				p[i] ^= p[i - j - 1];
			}
		}
	}

	for (i = 0; i < 100; i++)
	{
		temp = BIT(S, i*DSlength + plainIndex) ^ BIT(S, (i + 1)*DSlength + plainIndex) ^ BIT(S, i*DSlength + plainIndex + 1) ^ BIT(S, (i + 1)*DSlength + plainIndex + 1);
		b[i / 64] ^= (temp << (uint64_t)(63 - (i % 64)));

		Z[i] = (p[i*DSlength + plainIndex] ^ p[(i + 1)*DSlength + plainIndex] ^ p[i*DSlength + plainIndex + 1] ^ p[(i + 1)*DSlength + plainIndex + 1]) << (64 - check_deg);
	}

	IVtemp = GE2(Z, b, 100, check_deg, &i);
	for (i = 0; i < check_deg; i++)
	{
		IV ^= ((IVtemp >> (check_deg - i - 1)) & 1) << i;
	}

	free(p);
	return IV;
}

// Input  >> S : Stream after spreading and scrambling, streamlength : bit length of stream S, poly : primitive polynomial for scramble, deg : degree of primitive polynomial
// Output >> Stream that is descrambled.
// After spreading and scrambling stream S is descrambled by the estimated polynomial and initial value. And store it in S
void descramblingSync_poly(uint64_t *S, int streamlength, uint64_t poly, int deg, uint64_t IV)
{
	uint64_t *temp;
	uint64_t state;
	uint64_t s, tempS, degreemask;
	int i, j;

	temp = (uint64_t *)calloc((streamlength + 63) / 64, sizeof(uint64_t));

	degreemask = 0;
	for (i = 0; i < deg; i++)
	{
		degreemask ^= ((uint64_t)1 << i);
	}
	state = IV & degreemask;

	for (i = 0; i < streamlength; i++)
	{
		s = 0;
		tempS = state & (poly >> 1);

		for (j = 0; j < deg; j++)
		{
			s ^= (tempS >> j) & 1;
		}
		temp[i / 64] ^= (BIT(S, i) ^ s) << (63 - i % 64);
		state = ((state << 1) & degreemask) ^ s;
	}

	for (i = 0; i < streamlength / 64 + 1; i++)
	{
		S[i] = temp[i];
	}

	free(temp);
}
