#ifndef _DESCRAMBLE_H_
#define _DESCRAMBLE_H_

#include "general.h"
#include "Descrambling_Synchronous_Method1.h"
#include "Descrambling_Synchronous_method2.h"
#include "Descrambling_Synch_general.h"



// Input  >> S : Stream after spreading and scrambling, streamlength : bit length of S, use : use method1 then use =1, don't use method1 then use = 0
// Output >> findpoly : estimated polynomial, S : Descrambled Stream, DSlength : length of Spreading Code
// Use Method 1 and Method 2
int DESCRAMBLE(uint64_t *S, int streamlength, int DSlength, int WALSHorBARKERorPN, uint64_t *findpoly, int use);


#endif