/*
Descramble when we know the length of spreading code
*/
#include "general.h"
#include <math.h>
#include <time.h>

uint64_t primitivePolys[51] = { 0,0, 0x7, 0xB, 0x13, 0x25, 0x43, 0x83, 0x187, 0x211, 0x409, 0x805, 0x1407, 0x2129, 0x5803, 0x8003, 0x19401, 0x20009, 0x40081, 0x80609, 0x100009, 0x200005, 0x400003, 0x800021, 0x1000825, 0x2000009, 0x480A001, 0x8C20001, 0x10000009, 0x20000005, 0x48000601, 0x80000009, 0x100010085, 0x200002001, 0x400021101, 0x800000005, 0x1000000801, 0x2000404005, 0x4008000061, 0x8000000011, 0x10028800001, 0x20000000009, 0x404C0000001, 0x80008400021, 0x108800040001, 0x208010000011, 0x410080040001, 0x800000000021, 0x1000000080203, 0x2000000000201, 0x4000480020001 };
uint64_t barker[7] = { 0x8000000000000000,0xC000000000000000,0xD000000000000000,0xE800000000000000, 0xE400000000000000,0xE240000000000000,0xF9A8000000000000 };
uint64_t WALSH8 = 0x9600000000000000;
uint64_t WALSH16 = 0xC33C000000000000;
uint64_t WALSH32 = 0xF0F00F0F00000000;
uint64_t WALSH64 = 0xA55A5AA55AA5A55A;
uint64_t WALSH128[2] = { 0x9966669966999966,0x6699996699666699 };
uint64_t WALSH256[4] = { 0xCC33CC3333CC33CC,0x33CC33CCCC33CC33,0x33CC33CCCC33CC33,0xCC33CC3333CC33CC };
uint64_t WALSH512[8] = { 0xA5A55A5AA5A55A5A ,0xA5A55A5AA5A55A5A ,0x5A5AA5A55A5AA5A5,0x5A5AA5A55A5AA5A5,0x5A5AA5A55A5AA5A5,0x5A5AA5A55A5AA5A5,0xA5A55A5AA5A55A5A,0xA5A55A5AA5A55A5A };
uint64_t PN2 = 0x6000000000000000;
uint64_t PN3 = 0x7400000000000000;
uint64_t PN4 = 0x7ACB000000000000;
uint64_t PN5 = 0x7CD215D800000000;
uint64_t PN6 = 0x7EACDDA4E2F28C20;
uint64_t PN7[2] = { 0x7F54CEE9637B5B24, 0x70BE57344F143040 };
uint64_t PN8[4] = { 0x7FB678D721EEF410, 0x2DF37174CAA48A13, 0x47D694FC06D44BCB, 0x1198E183B0AC9CEA };
uint64_t PN9[8] = { 0x7FC3DC2CDBD0E612, 0x2BAF25CE0774F528, 0x155F5A0DDB582EF8, 0xF34D71A2FE96298C, 0x066564FDA49BF2D4, 0x289D97B0D5390C42, 0x011191D5B1C4A8D9, 0xF3C5B94826747DE0 };

int main(void)
{
	int DSlength, n, n_PN;
	int WALSHorBARKERorPN;
	int i, j, k, cnt, temp;
	int X, Y;
	int success;
	int notsynch;		
	int use, method;			
	uint64_t S[2048] = { 0 };
	uint64_t RT[8] = { 0 };
	uint64_t R[2048] = { 0 };
	uint64_t m[16] = { 0 };
	uint64_t tempSpreading[2048] = { 0 };
	uint64_t checkpoly, findpoly;
			
	
	printf("-----------------------------------------------------------------------------\n");
	printf("                                                                             \n");
	printf("-----------------------------------------------------------------------------\n");
	srand(time(NULL));

	///// Select Spreading Code /////
	while (1)
	{
		printf("\nSelect Spreading Code for making Stream. (1. WALSH , 2. Barker, 3. PNcode)  >>>  ");
		scanf_s("%d", &WALSHorBARKERorPN);
		if (WALSHorBARKERorPN == WALSH || WALSHorBARKERorPN == BARKER || WALSHorBARKERorPN == PN)
			break;
		printf("Invalid Input.\n");
	}

	if (WALSHorBARKERorPN == WALSH)
	{
		while (1)
		{
			printf("\nSelect the length of WALSH spreading code. (8,16,32,64,128,256,512)  >>>  ");
			scanf_s("%d", &DSlength);
			if (DSlength == 8 || DSlength == 16 || DSlength == 32 || DSlength == 64 || DSlength == 128 || DSlength == 256 || DSlength == 512)
				break;
			printf("Invalid Input.\n");
		}
	}
	else if (WALSHorBARKERorPN == BARKER)
	{
		while (1)
		{
			printf("\nSelect the length of BARKER spreading code. (2,3,4,5,7,11,13)  >>>  "); 
			scanf_s("%d", &DSlength);
			if (DSlength == 2 || DSlength == 3 || DSlength == 4 || DSlength == 5 || DSlength == 7 || DSlength == 11 || DSlength == 13)
				break;
			printf("Invalid Input.\n");
		}
	}
	else if (WALSHorBARKERorPN == PN)
	{
		while (1)
		{
			printf("\nSelect the length of PN spreading code (2, 3, 4, 5, 6, 7, 8, 9)  >>>  ");  .
			scanf_s("%d", &n_PN);
			if (n_PN >= 2 && n_PN <= 9 )
			{
				DSlength = pow(2, n_PN) - 1;
				break;
			}
			printf("Invalid Input.\n");
		}
	}

	while (1)
	{
		printf("\nSelect the degree of primitive polynomial for Synchronizing Scrambler (integer from 5 to 20) >>>  "); 
		scanf_s("%d", &n);
		if (n >= 5 && n <= 20)
			break;
		printf("Invalid Input.\n");
	}

	///// Make Random Message ///// 
	for (i = 0; i < 16; i++)
	{
		m[i] = (((uint64_t)rand() & 0xFFFF) << 48) ^ (((uint64_t)rand() & 0xFFFF) << 32) ^ (((uint64_t)rand() & 0xFFFF) << 16) ^ (((uint64_t)rand() & 0xFFFF));
	}

	printf("256 bits random message used for spreading : \n"); 
	printBIT(m, 0, 0, 256);

	///// Spreadinig ////// 
	// 1. WALSH case
	if (WALSHorBARKERorPN == WALSH && DSlength == 8)
	{
		printf("\n%dbits WALSH code used for spreading : \n", DSlength);
		printBIT(&WALSH8, 0, 0, DSlength);
		makeSyncStr(S, &WALSH8, DSlength, m, 2048, primitivePolys[n], n, 0x123456789);
		Spreading(tempSpreading, &WALSH8, DSlength, m, 2048);
	}
	else if (WALSHorBARKERorPN == WALSH && DSlength == 16)
	{
		printf("\n%dbits WALSH code used for spreading : \n", DSlength);
		printBIT(&WALSH16, 0, 0, DSlength);
		makeSyncStr(S, &WALSH16, DSlength, m, 2048, primitivePolys[n], n, 0x123456789);
		Spreading(tempSpreading, &WALSH16, DSlength, m, 2048);
	}
	else if (WALSHorBARKERorPN == WALSH && DSlength == 32)
	{
		printf("\n%dbits WALSH code used for spreading : \n", DSlength);
		printBIT(&WALSH32, 0, 0, DSlength);
		makeSyncStr(S, &WALSH32, DSlength, m, 2048, primitivePolys[n], n, 0x123456789);
		Spreading(tempSpreading, &WALSH32, DSlength, m, 2048);
	}
	else if (WALSHorBARKERorPN == WALSH && DSlength == 64)
	{
		printf("\n%dbits WALSH code used for spreading : \n", DSlength);
		printBIT(&WALSH64, 0, 0, DSlength);
		makeSyncStr(S, &WALSH64, DSlength, m, 2048, primitivePolys[n], n, 0x123456789);
		Spreading(tempSpreading, &WALSH64, DSlength, m, 2048);
	}
	else if (WALSHorBARKERorPN == WALSH && DSlength == 128)
	{
		printf("\n%dbits WALSH code used for spreading : \n", DSlength);
		printBIT(WALSH128, 0, 0, DSlength);
		makeSyncStr(S, WALSH128, DSlength, m, 1024, primitivePolys[n], n, 0x123456789);
		Spreading(tempSpreading, WALSH128, DSlength, m, 1024);
	}
	else if (WALSHorBARKERorPN == WALSH && DSlength == 256)
	{
		printf("\n%dbits WALSH code used for spreading : \n", DSlength);
		printBIT(WALSH256, 0, 0, DSlength);
		makeSyncStr(S, WALSH256, DSlength, m, 512, primitivePolys[n], n, 0x123456789);
		Spreading(tempSpreading, WALSH256, DSlength, m, 512);
	}
	else if (WALSHorBARKERorPN == WALSH && DSlength == 512)
	{
		printf("\n%dbits WALSH code used for spreading : \n", DSlength);
		printBIT(WALSH512, 0, 0, DSlength);
		makeSyncStr(S, WALSH512, DSlength, m, 256, primitivePolys[n], n, 0x123456789);
		Spreading(tempSpreading, WALSH512, DSlength, m, 256);
	}
	// 2. BARKER case
	else if (WALSHorBARKERorPN == BARKER && DSlength == 2)
	{
		printf("\n%dbits BARKER code used for spreading : \n", DSlength);
		//printf("\n확산에 사용된 %dbit barker code : \n", DSlength);
		printBIT(&barker[0], 0, 0, DSlength);
		makeSyncStr(S, &barker[0], DSlength, m, 2048, primitivePolys[n], n, 0x123456789);
		Spreading(tempSpreading, &barker[0], DSlength, m, 2048);
	}
	else if (WALSHorBARKERorPN == BARKER && DSlength == 3)
	{
		printf("\n%dbits BARKER code used for spreading : \n", DSlength);
		printBIT(&barker[1], 0, 0, DSlength);
		makeSyncStr(S, &barker[1], DSlength, m, 2048, primitivePolys[n], n, 0x123456789);
		Spreading(tempSpreading, &barker[1], DSlength, m, 2048);
	}
	else if (WALSHorBARKERorPN == BARKER && DSlength == 4)
	{
		printf("\n%dbits BARKER code used for spreading : \n", DSlength);
		printBIT(&barker[2], 0, 0, DSlength);
		makeSyncStr(S, &barker[2], DSlength, m, 2048, primitivePolys[n], n, 0x123456789);
		Spreading(tempSpreading, &barker[2], DSlength, m, 2048);
	}
	else if (WALSHorBARKERorPN == BARKER && DSlength == 5)
	{
		printf("\n%dbits BARKER code used for spreading : \n", DSlength);
		printBIT(&barker[3], 0, 0, DSlength);
		makeSyncStr(S, &barker[3], DSlength, m, 2048, primitivePolys[n], n, 0x123456789);
		Spreading(tempSpreading, &barker[3], DSlength, m, 2048);
	}
	else if (WALSHorBARKERorPN == BARKER && DSlength == 7)
	{
		printf("\n%dbits BARKER code used for spreading : \n", DSlength);
		printBIT(&barker[4], 0, 0, DSlength);
		makeSyncStr(S, &barker[4], DSlength, m, 2048, primitivePolys[n], n, 0x123456789);
		Spreading(tempSpreading, &barker[4], DSlength, m, 2048);
	}
	else if (WALSHorBARKERorPN == BARKER && DSlength == 11)
	{
		printf("\n%dbits BARKER code used for spreading : \n", DSlength);
		printBIT(&barker[5], 0, 0, DSlength);
		makeSyncStr(S, &barker[5], DSlength, m, 2048, primitivePolys[n], n, 0x123456789);
		Spreading(tempSpreading, &barker[5], DSlength, m, 2048);
	}
	else if (WALSHorBARKERorPN == BARKER && DSlength == 13)
	{
		printf("\n%dbits BARKER code used for spreading : \n", DSlength);
		printBIT(&barker[6], 0, 0, DSlength);
		makeSyncStr(S, &barker[6], DSlength, m, 2048, primitivePolys[n], n, 0x123456789);
		Spreading(tempSpreading, &barker[6], DSlength, m, 2048);
	}
	// 3. PN case
	else if (WALSHorBARKERorPN == PN && n_PN == 2)
	{
		printf("\n%dbits PN code used for spreading : \n", DSlength);
		printBIT(&PN2, 0, 0, DSlength);
		makeSyncStr(S, &PN2, DSlength, m, 2048, primitivePolys[n], n, 0x123456789);
		Spreading(tempSpreading, &PN2, DSlength, m, 2048);
	}
	else if (WALSHorBARKERorPN == PN && n_PN == 3)
	{
		printf("\n%dbits PN code used for spreading : \n", DSlength);
		printBIT(&PN3, 0, 0, DSlength);
		makeSyncStr(S, &PN3, DSlength, m, 2048, primitivePolys[n], n, 0x123456789);
		Spreading(tempSpreading, &PN3, DSlength, m, 2048);
	}
	else if (WALSHorBARKERorPN == PN && n_PN == 4)
	{
		printf("\n%dbits PN code used for spreading : \n", DSlength);
		printBIT(&PN4, 0, 0, DSlength);
		makeSyncStr(S, &PN4, DSlength, m, 2048, primitivePolys[n], n, 0x123456789);
		Spreading(tempSpreading, &PN4, DSlength, m, 2048);
	}
	else if (WALSHorBARKERorPN == PN && n_PN == 5)
	{
		printf("\n%dbits PN code used for spreading : \n", DSlength);
		printBIT(&PN5, 0, 0, DSlength);
		makeSyncStr(S, &PN5, DSlength, m, 2048, primitivePolys[n], n, 0x123456789);
		Spreading(tempSpreading, &PN5, DSlength, m, 2048);
	}
	else if (WALSHorBARKERorPN == PN && n_PN == 6)
	{
		printf("\n%dbits PN code used for spreading : \n", DSlength);
		printBIT(&PN6, 0, 0, DSlength);
		makeSyncStr(S, &PN6, DSlength, m, 2048, primitivePolys[n], n, 0x123456789);
		Spreading(tempSpreading, &PN6, DSlength, m, 2048);
	}
	else if (WALSHorBARKERorPN == PN && n_PN == 7)
	{
		printf("\n%dbits PN code used for spreading : \n", DSlength);
		printBIT(&PN7, 0, 0, DSlength);
		makeSyncStr(S, &PN7, DSlength, m, 512, primitivePolys[n], n, 0x123456789);
		Spreading(tempSpreading, &PN7, DSlength, m, 512);
	}
	else if (WALSHorBARKERorPN == PN && n_PN == 8)
	{
		printf("\n%dbits PN code used for spreading : \n", DSlength);
		printBIT(&PN8, 0, 0, DSlength);
		makeSyncStr(S, &PN8, DSlength, m, 256, primitivePolys[n], n, 0x123456789);
		Spreading(tempSpreading, &PN8, DSlength, m, 256);
	}
	else if (WALSHorBARKERorPN == PN && n_PN == 9)
	{
		printf("\n%dbits PN code used for spreading : \n", DSlength);
		printBIT(&PN9, 0, 0, DSlength);
		makeSyncStr(S, &PN9, DSlength, m, 128, primitivePolys[n], n, 0x123456789);
		Spreading(tempSpreading, &PN9, DSlength, m, 128);
	}

	////////////////////////////////////////// Descramble /////////////////////////////////////////////////
	printf("\nDescramble\n"); 
	while (1)
	{
		printf("\n\nIf you want to use method1 with method2, then select 1. Otherwise, select 2. \n(1. Method1 & Method2,  2. Method2 only)  >>>  ");
		scanf_s("%d", &use);
		if (use == 1)
		{
			printf("\nUse both Method1 & Method2.\n");
			break;
		}
		else if (use == 2)
		{
			printf("\nUse Method2 only.\n");
			break;
		}
		printf("Invalid Input.\n");
	}
	
	success = DESCRAMBLE(S, 1500 * 64, DSlength, WALSHorBARKERorPN, &findpoly, use);

	checkpoly = primitivePolys[n];
	

	if (success != FAIL && DSlength <= 64)
	{
		printf("\n\nThe result of descrambling with the estimated polynomial and IV : \n");
		for (i = 0; i < 256; )
		{
			printBIT(S, i, DSlength, (64 / DSlength)*DSlength);
			i += (64 / DSlength)*DSlength;
		}

		printf("\n\nThe Stream after Spreaded : \n");
		for (i = 0; i < 256; )
		{
			printBIT(tempSpreading, i, DSlength, (64 / DSlength)*DSlength);
			i += (64 / DSlength)*DSlength;
		}
	}
	else if (success != FAIL)
	{
		printf("\n\nThe result of descrambling with the estimated polynomial and IV : \n");
		for (i = 0; i < 512; i += DSlength)
		{
			printBIT(S, i, 0, DSlength);
		}

		printf("\n\nThe Stream after Spreaded : \n");
		for (i = 0; i < 512; i += DSlength)
		{
			printBIT(tempSpreading, i, 0, DSlength);
		}
	}

	////////////////////////////////////// Estimate Message ///////////////////////////////////////////
	if (success != FAIL)
	{
		// For delete trash values, delete some bits(multiple of length of spreading code) from MSB
		X = success;
		j = -1;
		while (X > 0) {
			X >>= 1;
			j += 1;
		}

		X = DSlength * ((j / DSlength) + 1); // deleting bits length : 1. longer than primitive polynomial's degree & 2. multiple of length of spreading code 
		Y = X / DSlength;

		if ((X % 64) == 0) {
			for (i = 0; i < 255; i++) {
				S[i] = S[i + (X / 64)];
			}
		}
		else {
			for (i = 0; i < 255; i++) {
				S[i] = (S[i + (X / 64)] << (X % 64)) | (S[i + (X / 64) + 1] >> (64 - (X % 64)));
			}
		}
		
		// Restore message from length of spreading code  
		// BARKER case
		if (WALSHorBARKERorPN == BARKER) 
		{
			printf("\n%d length BARKER code for using spreading(estimated) : \n", DSlength);
			if ((S[0] >> 63) & 1 == 1) 
			{
				RT[0] = S[0];
				printBIT(RT, 0, 0, DSlength);
			}
			else 
			{
				RT[0] = ~S[0];
				printBIT(RT, 0, 0, DSlength);
			}
			RT[0] >>= (64 - DSlength);
			
			printf("\n%d length BARKER code for using spreading(Actually used) : \n", DSlength); 
			if (DSlength == 2) printBIT(&barker[0], 0, 0, DSlength);
			if (DSlength == 3) printBIT(&barker[1], 0, 0, DSlength);
			if (DSlength == 4) printBIT(&barker[2], 0, 0, DSlength);
			if (DSlength == 5) printBIT(&barker[3], 0, 0, DSlength);
			if (DSlength == 7) printBIT(&barker[4], 0, 0, DSlength);
			if (DSlength == 11) printBIT(&barker[5], 0, 0, DSlength);
			if (DSlength == 13) printBIT(&barker[6], 0, 0, DSlength);

			j = 64;
			for (i = 0; i < 64; i++)
			{
				if (j >= DSlength)
				{
					R[(i*DSlength) / 64] = (R[(i*DSlength) / 64] << DSlength) ^ RT[0];
				}
				else
				{
					R[(i * DSlength) / 64] = (R[(i * DSlength) / 64] << j) ^ (RT[0] >> (DSlength - j));
					R[((i * DSlength) + j) / 64] = RT[0] & (((uint64_t)0x01 << (DSlength)) - 1);
				}
				j -= DSlength;
				if (j <= 0) { j += 64; }
			}
			R[(i*DSlength) / 64] <<= j;

			for (i = 0; i < 32; i++)
			{
				R[i] ^= S[i];
			}
			
			printf("\nRestored random message :\n"); 
			for (i = 0; i < Y; i++) printf(" ");
			for (i = 0; i < 32 - Y; i++)
			{
				k = 0;
				for (j = 0; j < DSlength; j++)
				{
					k += (R[(i*DSlength + j) / 64] >> (63 - ((i*DSlength + j) % 64))) & (uint64_t)0x01;
				}
				if (k == 0) printf("%d", 1);
				else if (k == DSlength) printf("%d", 0);
				else printf("Wrong spreading code\n"); 
			}
			printf("\n\nActually used message : \n"); 
			printBIT(m, 0, 0, 32);
		}
		// WALSH case
		else if (WALSHorBARKERorPN == WALSH) 
		{
			printf("\n%d length WALSH code for using spreading(Actually used) :\n", DSlength); 
			
			
			if ((S[0] >> 63) & 1 == 1)
			{
				for (i = 0; i < (DSlength + 63) / 64; i++)
				{
					RT[i] ^= S[i];
				}
				printBIT(RT, 0, 0, DSlength);
			}
			else
			{
				for (i = 0; i < (DSlength + 63) / 64; i++)
				{
					RT[i] ^= ~S[i];
				}
				printBIT(RT, 0, 0, DSlength);
			}
			if (DSlength < 64) RT[0] >>= (64 - DSlength);
																							
			printf("\n\n%d length WALSH code for using spreading(Actually used) :\n", DSlength);
			if (DSlength == 8) printBIT(&WALSH8, 0, 0, DSlength);
			if (DSlength == 16) printBIT(&WALSH16, 0, 0, DSlength);
			if (DSlength == 32) printBIT(&WALSH32, 0, 0, DSlength);
			if (DSlength == 64) printBIT(&WALSH64, 0, 0, DSlength);
			if (DSlength == 128) printBIT(WALSH128, 0, 0, DSlength);
			if (DSlength == 256) printBIT(WALSH256, 0, 0, DSlength);
			if (DSlength == 512) printBIT(WALSH512, 0, 0, DSlength);
			
			// Length of WALSH code is shorter than 64 
			if (DSlength < 64) 
			{
				j = 64;
				for (i = 0; i < 64; i++)
				{
					if (j >= DSlength)
					{
						R[(i * DSlength) / 64] = (R[(i * DSlength) / 64] << DSlength) ^ RT[0];
					}
					else
					{
						R[(i * DSlength) / 64] = (R[(i * DSlength) / 64] << j) ^ (RT[0] >> (DSlength - j));
						R[((i * DSlength) + j) / 64] = RT[0] & (((uint64_t)0x01 << (DSlength - 1)) - 1);
					}
					j -= DSlength;
					if (j <= 0) j += 64;
				}
				R[(i * DSlength) / 64] <<= j;



				for (i = 0; i < 32; i++)
				{
					R[i] ^= S[i];
				}

				printf("\nRestored random message :\n"); 
				for (i = 0; i < Y; i++) printf(" ");	
				for (i = 0; i < 32 - Y; i++)
				{
					for (j = 0, k = 0; j < DSlength; j++)
					{
						k += (R[(i*DSlength + j) / 64] >> (63 - ((i*DSlength + j) % 64))) & (uint64_t)0x01;
					}
					if (k == 0) printf("%d", 1);
					else if (k == DSlength) printf("%d", 0);
					else printf("Wrong spreading code\n");

				}
				printf("\n\nActually used message : \n"); 
				printBIT(m, 0, 0, 32);
			}
			else // Length of WALSH code is longer than 64
			{
				for (i = 0; i < 32; i++)
				{
					for (j = 0; j < DSlength / 64; j++)
					{
						R[(DSlength / 64) * i + j] = S[(DSlength / 64) * i + j] ^ RT[j];
					}
				}
				printf("\nRestored random message :\n");
				for (i = 0; i < Y; i++) printf(" ");	
				for (i = 0; i < 32 - Y; i++)
				{
					for (j = 0, k = 0; j < DSlength; j++)
					{
						k += BIT(R, (i*DSlength + j));
					}
					if (k == 0) printf("%d", 1);
					else if (k == DSlength) printf("%d", 0);
					else printf("Wrong spreading code\n");

				}
				printf("\n\nActually used message : \n");
				printBIT(m, 0, 0, 32);
			}
		}
		// PN case
		else if (WALSHorBARKERorPN == PN)
		{
			printf("\n%d length BARKER code for using spreading(estimated) :\n", DSlength);

			cnt = 0;
			for (i = 0; i < DSlength; i++) //1 in PN code = 0 in PN code + 1 
			{
				if (BIT(S, i) == 0)
				{
					cnt += 1;
				}
			}
			if (cnt == (DSlength / 2)) 
			{
				for (i = 0; i < (DSlength + 63) / 64; i++)
				{
					RT[i] ^= S[i];
				}
				printBIT(RT, 0, 0, DSlength);
			}
			else
			{
				for (i = 0; i < (DSlength + 63) / 64; i++)
				{
					RT[i] ^= ~S[i];
				}
				printBIT(RT, 0, 0, DSlength);
			}
			if (DSlength < 64) RT[0] >>= (64 - DSlength);

			printf("\n\n%d length BARKER code for using spreading(Actually used) :\n");
			if (DSlength == 3) printBIT(&PN2, 0, 0, DSlength);
			if (DSlength == 7) printBIT(&PN3, 0, 0, DSlength);
			if (DSlength == 15) printBIT(&PN4, 0, 0, DSlength);
			if (DSlength == 31) printBIT(&PN5, 0, 0, DSlength);
			if (DSlength == 63) printBIT(&PN6, 0, 0, DSlength);
			if (DSlength == 127) printBIT(&PN7, 0, 0, DSlength);
			if (DSlength == 255) printBIT(&PN8, 0, 0, DSlength);
			if (DSlength == 511) printBIT(&PN8, 0, 0, DSlength);

			// Length of PN code is shorter than 64 
			if (DSlength < 64)
			{
				j = 64;
				for (i = 0; i < 64; i++)
				{
					if (j >= DSlength)
					{
						R[(i * DSlength) / 64] = (R[(i * DSlength) / 64] << DSlength) ^ RT[0];
					}
					else
					{
						R[(i * DSlength) / 64] = (R[(i * DSlength) / 64] << j) ^ (RT[0] >> (DSlength - j));
						R[((i * DSlength) + j) / 64] = RT[0] & (((uint64_t)0x01 << (DSlength - 1)) - 1);
					}
					j -= DSlength;
					if (j <= 0) j += 64;
				}
				R[(i * DSlength) / 64] <<= j;

				for (i = 0; i < 32; i++)
				{
					R[i] ^= S[i];
				}

				printf("\nRestored random message:\n");
				for (i = 0; i < Y; i++) printf(" ");	
				for (i = 0; i < 32 - Y; i++)
				{
					for (j = 0, k = 0; j < DSlength; j++)
					{
						k += (R[(i*DSlength + j) / 64] >> (63 - ((i*DSlength + j) % 64))) & (uint64_t)0x01;
					}
					if (k == 0) printf("%d", 1);
					else if (k == DSlength) printf("%d", 0);
					else printf("Wrong spreading code\n");

				}
				printf("\n\nActually used message : \n");
				printBIT(m, 0, 0, 32);
			}
			// Length of PN code is longer than 64 
			else 
			{
				j = 64;
				for (i = 0; i < 64 * 512; i++)
				{
					j = i % DSlength;
					temp = BIT(RT, j);
					InsertBIT(R, temp, i);
				}

				for (i = 0; i < 512; i++)
				{
					R[i] ^= S[i];
				}

				printf("\nRestored random message:\n");
				for (i = 0; i < Y; i++) printf(" ");	
				for (i = 0; i < 32 - Y; i++)
				{
					for (j = 0, k = 0; j < DSlength; j++)
					{
						k += BIT(R, (i*DSlength + j)); 
					}
					if (k == 0) printf("%d", 1);
					else if (k == DSlength) printf("%d", 0);
					else printf("Wrong spreading code\n");

				}
				printf("\n\nActually used message : \n");
				printBIT(m, 0, 0, 32);
			}
		}
	}

	system("pause");
	
	return 0;

}
