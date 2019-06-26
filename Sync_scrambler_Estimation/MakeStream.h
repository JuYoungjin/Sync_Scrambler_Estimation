#ifndef _MAKESTREAM_H_
#define _MAKESTREAM_H_

#include <stdio.h>
#include <stdint.h>


// Input  >> DS : Spreading Code, Dslength : length of DS, m : message, mlength : bit length of message
// Output >> T : Stream that is made from message by Spreading
// By using spreading Code DS and message m, Save spreaded stream to S (input S must be filled with 0)
void Spreading(uint64_t *T, uint64_t *DS, int DSlength, uint64_t *m, int mlength);

// Input  >> DS : Spreading Code, Dslength : length of DS, m : message, mlength : bit length of message, poly : Scrambling polynomial, degree : degree of poly, IV : initial value for Scrmabling
// Output >> S : Stream that is made from message by Spreading and Synch Scrambling
// Save stream S that is made from message by spreading and scrambling
void makeSyncStr(uint64_t *S, uint64_t *DS, int DSlength, uint64_t *m, int mlength, uint64_t poly, int degree, uint64_t IV);


#endif