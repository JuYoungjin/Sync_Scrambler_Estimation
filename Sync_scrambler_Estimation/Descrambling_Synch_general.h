#ifndef _DESCRAMBLING_SYNCHRONOUS_GENERAL_H_
#define _DESCRAMBLING_SYNCHRONOUS_GENERAL_H_
#include "general.h"

// plainIndex = randomIndex-2
// randomIndex : bit position where message was not deleted
// plainIndex = randomIndex - 2
int randomIndex, plainIndex;


// Input  >> S : Stream after spreading and scrambling, poly : primitive polynomial for scramble, DSlength : length of Spreading code
// Output >> Initial Value of Scrambling & degree of primitive polynomial
// Algorithm for finding IV by using primitive polynomial and stream S
uint64_t FindInitialvlaue(uint64_t *S, uint64_t poly, int DSlength, int *deg);

// Input  >> S : Stream after spreading and scrambling, streamlength : bit length of stream S, poly : primitive polynomial for scramble, deg : degree of primitive polynomial
// Output >> Stream that is descrambled.
// After spreading and scrambling stream S is descrambled by the estimated polynomial and initial value. And store it in S
void descramblingSync_poly(uint64_t *S, int streamlength, uint64_t poly, int deg, uint64_t IV);

#endif