#ifndef _DESCRAMBLINGSYNC_H_
#define _DESCRAMBLINGSYNC_H_
#include "general.h"
#include "Descrambling_Synch_general.h"


// Input  >> S : Stream after spreading and scrambling, streamlength : bit length of stream S, DSlength : length of Spreading code, WALSHorBARKERorPN : Config that indicate Spreading Code
// Output >> Primitive polynomial for scrambling(stored in findpoly), Length of Spreading code(stored in DSlength), and descrambled stream S
// Descramble S after estimating primitive polynomial and initial value.
int DescramblingSync_Method1(uint64_t *S, int streamlength, int DSlength, int WALSHorBARKERorPN);


// Input  >> S : Stream after spreading and scrambling, streamlength : bit length of stream S
// Output >> Primitive polynomial for scrambling & degree of primitive polynomial
// Estimate primitive polynomial and degree of primitive polynomial by using Method1 from S that is scrambled and spreaded by WALSH code.
uint64_t FindSyncScramblerPolynomial_Method1_WALSH(uint64_t *S, int streamlength, int *deg);


// Input  >> S : Stream after spreading and scrambling, streamlength : bit length of stream S, DSlength : length of spreading code
// Output >> Primitive polynomial for scrambling & degree of primitive polynomial
// Estimate primitive polynomial and degree of primitive polynomial by using Method1 from length of spreading code and S that is scrambled and spreaded by BARKER code.
uint64_t FindSyncScramblerPolynomial_Method1_BARKER(uint64_t *S, int streamlength, int DSlength, int *deg);

// Input  >> S : Stream after spreading and scrambling, streamlength : bit length of stream S, DSlength : length of spreading code(estimated)
// Output >> Primitive polynomial for scrambling
// Estimate primitive polynomial and degree of primitive polynomial by using Method1 from length of spreading code S that is scrambled and spreaded by PN code.
uint64_t FindSyncScramblerPolynomial_Method1_PN(uint64_t *S, int streamlength, int DSlength, int *deg);


#endif