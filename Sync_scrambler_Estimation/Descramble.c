#include "Descramble.h"

// Input  >> S : Stream after spreading and scrambling, streamlength : bit length of S, use : use method1 then use =1, don't use method1 then use = 0
// Output >> findpoly : estimated polynomial, S : Descrambled Stream, DSlength : length of Spreading Code
// Use Method 1 and Method 2
int DESCRAMBLE(uint64_t *S, int streamlength, int DSlength, int WALSHorBARKERorPN, uint64_t *findpoly, int use)
{
	uint64_t poly;
	if (use == 1)
	{
		printf("\nMETHOD1\n");
		poly = DescramblingSync_Method1(S, streamlength, DSlength, WALSHorBARKERorPN);
		if (poly != FAIL) {
			*findpoly = poly;
			return SUCCESS;
		}
		
		printf("\nMETHOD2\n");
		poly = descramblingSync_Method2(S, streamlength, DSlength, WALSHorBARKERorPN);
		if (poly != FAIL) {
			*findpoly = poly;
			return SUCCESS;
		}
		else if (poly == FAIL) {
			printf("\nPolynomial Estimation Fail\n");
			return FAIL;
		}
	}
	else if (use == 2)
	{
		printf("\nMETHOD2\n");
		poly = descramblingSync_Method2(S, streamlength, DSlength, WALSHorBARKERorPN);
		if (poly != FAIL) {
			*findpoly = poly;
			return SUCCESS;
		}
		else if (poly == FAIL) {
			printf("\nPolynomial Estimation Fail\n");
			return FAIL;
		}
	}
}