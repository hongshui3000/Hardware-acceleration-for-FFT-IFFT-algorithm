/*
 * fft.c
 *
 *  Created on: 08.03.2019
 *      Author: yukai.shen
 */

#include "fft.h"

extern short SampleArray_Real[FFTSIZE];
extern short SampleArray_Imag[FFTSIZE];
extern short W_Real[FFTSIZE/2];
extern short W_Imag[FFTSIZE/2];
short Sinewave[1024];


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 														FFT without TIE														  	//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//===
//DIT
//===
#if FFT==0
void DIT()
{
	int temp0_real, temp0_imag, temp1_real, temp1_imag;
	short NumberOfStage, StageIndex;
	short counter, step, node;

	NumberOfStage = POWER;

	Bitreverse();

	for(StageIndex=0; StageIndex<NumberOfStage; StageIndex++){
		counter = 0;	//counter for the Twiddle Factors
		step = 1<<StageIndex;	//Step length of the two sample nodes in each butterfly
		GetTwiddleFactors(StageIndex);

		for(node=0; node<FFTSIZE; node++){
			//distinguish by choosing the bottom node in each butterfly
			if((node>>StageIndex)%2 == 1){
				//Scaling up the input sample
				temp0_real = (int)SampleArray_Real[node-step] << 15;
				temp0_imag = (int)SampleArray_Imag[node-step] << 15;

				temp1_real = (int)SampleArray_Real[node] * W_Real[counter] - (int)SampleArray_Imag[node] * W_Imag[counter];
				temp1_imag = (int)SampleArray_Real[node] * W_Imag[counter] + (int)SampleArray_Imag[node] * W_Real[counter];

				//Scaling down the results and place them back into the SampleArray ready for the next stage
				SampleArray_Real[node-step] = (short)((temp0_real + temp1_real) >> 15);
				SampleArray_Imag[node-step] = (short)((temp0_imag + temp1_imag) >> 15);
				SampleArray_Real[node] = (short)((temp0_real - temp1_real) >> 15);
				SampleArray_Imag[node] = (short)((temp0_imag - temp1_imag) >> 15);

				counter += 1;
			}
		}
	}
	if(DO_INVERSE) IFFT();
}
#endif

//===
//DIF
//===
#if FFT==1
void DIF()
{
	short temp0_real, temp0_imag, temp1_real, temp1_imag;
	short NumberOfStage, StageIndex;
	short counter, step, node;

	NumberOfStage = POWER;

	for(StageIndex=NumberOfStage-1; StageIndex>=0; StageIndex--){
		counter = 0;
		step = 1 << StageIndex;
		GetTwiddleFactors(StageIndex);

		for(node=0; node<FFTSIZE; node++){
			if((node>>StageIndex)%2==1){
				temp0_real = SampleArray_Real[node-step] + SampleArray_Real[node];
				temp0_imag = SampleArray_Imag[node-step] + SampleArray_Imag[node];
				temp1_real = SampleArray_Real[node-step] - SampleArray_Real[node];
				temp1_imag = SampleArray_Imag[node-step] - SampleArray_Imag[node];

				SampleArray_Real[node-step] = temp0_real;
				SampleArray_Imag[node-step] = temp0_imag;
				SampleArray_Real[node] = (int)((int)temp1_real * W_Real[counter] - (int)temp1_imag * W_Imag[counter]) >> 15;
				SampleArray_Imag[node] = (int)((int)temp1_real * W_Imag[counter] + (int)temp1_imag * W_Real[counter]) >> 15;

				counter += 1;
			}
		}
	}
	Bitreverse();
	if(DO_INVERSE) IFFT();
}
#endif

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 														FFT with MAC														  	//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#if FFT==2
//===================
//DIT_macTIE_no_table
//===================

void DIT_macTIE()
{
	fft_16 *Sample0Real, *Sample0Imag, *Sample1Real, *Sample1Imag;
	short NumberOfStage, StageIndex;
	short counter, step, node;

	NumberOfStage = POWER;
	Bitreverse();

	for(StageIndex=0; StageIndex<NumberOfStage; StageIndex++){
		counter = 0;
		step = 1<<StageIndex;
		GetTwiddleFactors(StageIndex);

		for(node=0; node<FFTSIZE; node++){
			if((node>>StageIndex)%2 == 1){
				Sample0Real = (fft_16*)&SampleArray_Real[node-step];
				Sample0Imag = (fft_16*)&SampleArray_Imag[node-step];
				Sample1Real = (fft_16*)&SampleArray_Real[node];
				Sample1Imag = (fft_16*)&SampleArray_Imag[node];

				DIT_one_BUTTERFLY(*Sample0Real, *Sample0Imag, *Sample1Real, *Sample1Imag, W_Real[counter], W_Imag[counter]);
				counter += 1;
			}
		}
	}
	if(DO_INVERSE) IFFT();
}


//=====================
//DIT_macTIE_with_table
//=====================
/*
void DIT_macTIE()
{
	fft_16 *Sample0Real, *Sample0Imag, *Sample1Real, *Sample1Imag;
	short NumberOfStage, StageIndex;
	short step, node;
	short NumberOfGroups;
	short Rightangel = 256;	//N_WAVE >> 2
	short imag_angel, real_angel;
	short GroupIndex;
	short j;

	NumberOfStage = POWER;

	Bitreverse();

	for(StageIndex=0; StageIndex<NumberOfStage; StageIndex++){
		NumberOfGroups = FFTSIZE >> (StageIndex+1);
		step = 1<<StageIndex;
			for (GroupIndex=0; GroupIndex<NumberOfGroups; GroupIndex++){
				for(j=0; j<step; j++){
					imag_angel = j<<(9-StageIndex);
					real_angel = imag_angel + Rightangel;
					node = (step<<1)*GroupIndex + j;
					Sample0Real = (fft_16*)&SampleArray_Real[node];
					Sample0Imag = (fft_16*)&SampleArray_Imag[node];
					Sample1Real = (fft_16*)&SampleArray_Real[node+step];
					Sample1Imag = (fft_16*)&SampleArray_Imag[node+step];
					DIT_one_BUTTERFLY(*Sample0Real, *Sample0Imag, *Sample1Real, *Sample1Imag, real_angel, imag_angel);
				}
			}
		}
}
*/
#endif


#if FFT==3
//===================
//DIF_macTIE_no_table
//===================

void DIF_macTIE()
{
	fft_16 *Sample0Real, *Sample0Imag, *Sample1Real, *Sample1Imag;
	short NumberOfStage, StageIndex;
	short counter, step, node;

	NumberOfStage = POWER;

	for(StageIndex=NumberOfStage-1; StageIndex>=0; StageIndex--){
		counter = 0;
		step = 1<<StageIndex;
		GetTwiddleFactors(StageIndex);

		for(node=0; node<FFTSIZE; node++){
			if((node>>StageIndex)%2 == 1){
				Sample0Real = (fft_16*)&SampleArray_Real[node-step];
				Sample0Imag = (fft_16*)&SampleArray_Imag[node-step];
				Sample1Real = (fft_16*)&SampleArray_Real[node];
				Sample1Imag = (fft_16*)&SampleArray_Imag[node];

				DIF_one_BUTTERFLY(*Sample0Real, *Sample0Imag, *Sample1Real, *Sample1Imag, W_Real[counter], W_Imag[counter]);
				counter += 1;
			}
		}
	}
	Bitreverse();
	if(DO_INVERSE) IFFT();
}


//=====================
//DIF_macTIE_with_table
//=====================
/*
void DIF_macTIE()
{
	fft_16 *Sample0Real, *Sample0Imag, *Sample1Real, *Sample1Imag;
	short NumberOfStage, StageIndex;
	short step, node;
	short NumberOfGroups;
	short Rightangel = 256;	//N_WAVE >> 2;
	//short Straightangel = 512;
	short imag_angel, real_angel;
	short GroupIndex;
	short NumberOfSamplesInGroup;
	short j;

	NumberOfStage = POWER;
	for(StageIndex=0; StageIndex<NumberOfStage; StageIndex++){
		NumberOfGroups = 1 << StageIndex;
		step = FFTSIZE  >> (StageIndex+1);
		NumberOfSamplesInGroup = FFTSIZE >> StageIndex;
		for(GroupIndex=0; GroupIndex<NumberOfGroups; GroupIndex++){
			for(j=0; j<step; j++){
				node = NumberOfSamplesInGroup*GroupIndex + j;
				Sample0Real = (fft_16*)&SampleArray_Real[node];
				Sample0Imag = (fft_16*)&SampleArray_Imag[node];
				Sample1Real = (fft_16*)&SampleArray_Real[node+step];
				Sample1Imag = (fft_16*)&SampleArray_Imag[node+step];
				imag_angel = j<<(10+StageIndex-NumberOfStage);
				real_angel = imag_angel + Rightangel;
				DIF_one_BUTTERFLY(*Sample0Real, *Sample0Imag, *Sample1Real, *Sample1Imag, real_angel, imag_angel);
			}
		}
	}
	Bitreverse();
}
*/

#endif



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 														FFT with SIMD														  	//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//===========
//DIT_simdTIE
//===========
#if FFT==4
void DIT_simdTIE()
{
	short StageIndex, counter, NumberOfStage;
	short sub_stage;	//Every 8 input samples as a sub stage in each stage;
	fft_128 *SampleReal, *SampleImag;
	fft_64 *WReal, *WImag;

	counter = 0;	//counter of sub_stage
	StageIndex = 0;
	NumberOfStage = POWER;
	Bitreverse();
	GetTwiddleFactors(StageIndex);
	//Calculate the StageIndex=0 at first, every time 4 butterflies parallel calculation
	for(sub_stage=0; sub_stage<(FFTSIZE/8); sub_stage++){
		if(counter >= FFTSIZE) break;
		SampleReal = (fft_128*)(SampleArray_Real+8*sub_stage);
		SampleImag = (fft_128*)(SampleArray_Imag+8*sub_stage);
		WReal = (fft_64*)(W_Real+4*sub_stage);
		WImag = (fft_64*)(W_Imag+4*sub_stage);
		DIT_four_BUTTERFLY(*SampleReal, *SampleImag, *WReal, *WImag);
		counter += 8;
	}

	//Calculate from StageIndex=1
	for(StageIndex=1; StageIndex<NumberOfStage; StageIndex++){
		counter = 0;
		Remap_1(StageIndex);
		GetTwiddleFactors(StageIndex);

		for(sub_stage=0; sub_stage<(FFTSIZE/8); sub_stage++){
			if(counter >= FFTSIZE) break;
			SampleReal = (fft_128*)(SampleArray_Real+8*sub_stage);
			SampleImag = (fft_128*)(SampleArray_Imag+8*sub_stage);
			WReal = (fft_64*)(W_Real+4*sub_stage);
			WImag = (fft_64*)(W_Imag+4*sub_stage);

			DIT_four_BUTTERFLY(*SampleReal, *SampleImag, *WReal, *WImag);
			counter += 8;
		}

		Remap_2(StageIndex);
	}
	if(DO_INVERSE) IFFT();
}
#endif

//===========
//DIF_simdTIE
//===========
#if FFT==5
void DIF_simdTIE()
{
	short StageIndex, counter, NumberOfStage;
	short sub_stage;
	fft_128 *SampleReal, *SampleImag;
	fft_64 *WReal, *WImag;

	counter = 0;
	NumberOfStage = POWER;

	for(StageIndex=NumberOfStage-1; StageIndex>0; StageIndex--){
		counter = 0;
		Remap_1(StageIndex);
		GetTwiddleFactors(StageIndex);
		for(sub_stage=0; sub_stage<(FFTSIZE/8); sub_stage++){
			if(counter >= FFTSIZE) break;
			SampleReal = (fft_128*)(SampleArray_Real+8*sub_stage);
			SampleImag = (fft_128*)(SampleArray_Imag+8*sub_stage);
			WReal = (fft_64*)(W_Real+4*sub_stage);
			WImag = (fft_64*)(W_Imag+4*sub_stage);
			DIF_four_BUTTERFLY(*SampleReal, *SampleImag, *WReal, *WImag);
			counter += 8;
		}
		Remap_2(StageIndex);
	}

	//Calculate StageIndex = 0
	counter = 0;
	StageIndex = 0;
	GetTwiddleFactors(StageIndex);
	for(sub_stage=0; sub_stage<(FFTSIZE/8); sub_stage++){
		if(counter >= FFTSIZE) break;
		SampleReal = (fft_128*)(SampleArray_Real+8*sub_stage);
		SampleImag = (fft_128*)(SampleArray_Imag+8*sub_stage);
		WReal = (fft_64*)(W_Real+4*sub_stage);
		WImag = (fft_64*)(W_Imag+4*sub_stage);
		DIF_four_BUTTERFLY(*SampleReal, *SampleImag, *WReal, *WImag);
		counter += 8;
	}
	Bitreverse();
	if(DO_INVERSE) IFFT();
}
#endif

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 													 	Subfunctions 															//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Bitreverse()
{
	short temp_real;
	short temp_imag;
	int i=0, j=0, k=0;
	int t;

	for(i=0; i<FFTSIZE; i++){
		k = i;							//k: SampleIndex before reversal
		j = 0;							//j: SampleIndex after reversal
		t = POWER;						//t: How many times need to do the reversal

		while( (t--)>0 ){
			j = j<<1;
			j |= (k&1);					//get the lowest bit of k and add it to j
			k = k>>1;
		}

		if(j>i){
			temp_real = SampleArray_Real[i];
			temp_imag = SampleArray_Imag[i];
			SampleArray_Real[i] = SampleArray_Real[j];
			SampleArray_Imag[i] = SampleArray_Imag[j];
			SampleArray_Real[j] = temp_real;
			SampleArray_Imag[j] = temp_imag;

		}
	}
}

void Output()
{
	int i;

	for(i=0; i<FFTSIZE; i++){
		printf("\n%d: real:%d\timag:%d\n", i, SampleArray_Real[i], SampleArray_Imag[i]);
	}
}


void GetTwiddleFactors(short StageIndex)
{
	short i,j,k, PowerStage, RightAngel, W_inStage;
	short angel,counter;

	RightAngel = N_WAVE >> 2;
	PowerStage = 1<<StageIndex;
	W_inStage = (FFTSIZE>>(StageIndex+1));
	k = 0;

	for(j=0; j<PowerStage; j++){
		angel = j<<(9-StageIndex);
		W_Real[k] = Sinewave[RightAngel+angel];
		if(DO_INVERSE)
			W_Imag[k] = Sinewave[angel];
		else
			W_Imag[k] = -Sinewave[angel];
		k = k + 1;
	}

	//Copy the twiddle factors for each butterfly
	for(i=1; i<W_inStage; i++){
		for(j=0; j<PowerStage; j++){
			W_Real[k] = W_Real[j];
			W_Imag[k] = W_Imag[j];
			k = k+1;
		}
	}
}

void IFFT()
{
	short i;
	for(i=0; i<FFTSIZE; i++){
		SampleArray_Real[i] = SampleArray_Real[i] >> POWER;
		SampleArray_Imag[i] = SampleArray_Imag[i] >> POWER;
	}

}

/*
 * The Operation that four Butterflies parallel calculate is designed based on DIT for StageIndex = 0 (or DIF for last stage).
 * Remap_1(): From the StageIndex = 1, need a rempping function to reorder the sample array
 * to make sure that TIE operation (four butterflies parallel calculate)can be applied.
 * Remap_2(): Inverse function of Remap_1(), reorder the results after the TIE operation calculation
 * and get ready for the next stage.
 */
void Remap_1(short StageIndex)
{
	short NumberOfGroups;	//Number of groups in each stage
	short RemapTime;		//How many times needed for remapping in each group
	short GroupIndex;		//Group index in each stage
	short temp_real[FFTSIZE], temp_imag[FFTSIZE];
	short node;				//Index for Sample array
	short counter;			//Index for input sample array while remapping
	short counter_temp;		//Index for temp array while remapping
	short i;

	NumberOfGroups = FFTSIZE >> (StageIndex+1);
	RemapTime = 1 << StageIndex;
	counter_temp = 0;
	for(GroupIndex=0; GroupIndex<NumberOfGroups; GroupIndex++){
		counter = counter_temp;
		for(i=0; i<RemapTime; i++){
			temp_real[counter_temp] = SampleArray_Real[counter];
			temp_imag[counter_temp] = SampleArray_Imag[counter];
			temp_real[counter_temp+1] = SampleArray_Real[counter+RemapTime];
			temp_imag[counter_temp+1] = SampleArray_Imag[counter+RemapTime];
			counter += 1;
			counter_temp += 2;
		}
	}

	//Put the after-remapping array from temp array to the input sample array
	for(node=0; node<FFTSIZE; node++){
		SampleArray_Real[node] = temp_real[node];
		SampleArray_Imag[node] = temp_imag[node];
	}
}

void Remap_2(short StageIndex)
{
	short NumberOfGroups;
	short RemapTime;
	short GroupIndex;
	short temp_real[FFTSIZE], temp_imag[FFTSIZE];
	short node;
	short counter;
	short counter_temp;
	short i;

	NumberOfGroups = FFTSIZE >> (StageIndex+1);
	RemapTime = 1 << StageIndex;
	counter = 0;
	for(GroupIndex=0; GroupIndex<NumberOfGroups; GroupIndex++){
		counter_temp = counter;
		for(i=0; i<RemapTime; i++){
			temp_real[counter_temp] = SampleArray_Real[counter];
			temp_imag[counter_temp] = SampleArray_Imag[counter];
			temp_real[counter_temp+RemapTime] = SampleArray_Real[counter+1];
			temp_imag[counter_temp+RemapTime] = SampleArray_Imag[counter+1];
			counter += 2;
			counter_temp += 1;
		}
	}
	//Put the after-remapping array from temp array to the input sample array
	for(node=0; node<FFTSIZE; node++){
		SampleArray_Real[node] = temp_real[node];
		SampleArray_Imag[node] = temp_imag[node];
	}
}

//Lookup tables for Twiddle Factors
//Stored in System RAM, Starting_Address: 0x60009d80
short Sinewave[1024] = {
      0,    201,    402,    603,    804,   1005,   1206,   1406,
   1607,   1808,   2009,   2209,   2410,   2610,   2811,   3011,
   3211,   3411,   3611,   3811,   4011,   4210,   4409,   4608,
   4807,   5006,   5205,   5403,   5601,   5799,   5997,   6195,
   6392,   6589,   6786,   6982,   7179,   7375,   7571,   7766,
   7961,   8156,   8351,   8545,   8739,   8932,   9126,   9319,
   9511,   9703,   9895,  10087,  10278,  10469,  10659,  10849,
  11038,  11227,  11416,  11604,  11792,  11980,  12166,  12353,
  12539,  12724,  12909,  13094,  13278,  13462,  13645,  13827,
  14009,  14191,  14372,  14552,  14732,  14911,  15090,  15268,
  15446,  15623,  15799,  15975,  16150,  16325,  16499,  16672,
  16845,  17017,  17189,  17360,  17530,  17699,  17868,  18036,
  18204,  18371,  18537,  18702,  18867,  19031,  19194,  19357,
  19519,  19680,  19840,  20000,  20159,  20317,  20474,  20631,
  20787,  20942,  21096,  21249,  21402,  21554,  21705,  21855,
  22004,  22153,  22301,  22448,  22594,  22739,  22883,  23027,
  23169,  23311,  23452,  23592,  23731,  23869,  24006,  24143,
  24278,  24413,  24546,  24679,  24811,  24942,  25072,  25201,
  25329,  25456,  25582,  25707,  25831,  25954,  26077,  26198,
  26318,  26437,  26556,  26673,  26789,  26905,  27019,  27132,
  27244,  27355,  27466,  27575,  27683,  27790,  27896,  28001,
  28105,  28208,  28309,  28410,  28510,  28608,  28706,  28802,
  28897,  28992,  29085,  29177,  29268,  29358,  29446,  29534,
  29621,  29706,  29790,  29873,  29955,  30036,  30116,  30195,
  30272,  30349,  30424,  30498,  30571,  30643,  30713,  30783,
  30851,  30918,  30984,  31049,
  31113,  31175,  31236,  31297,
  31356,  31413,  31470,  31525,  31580,  31633,  31684,  31735,
  31785,  31833,  31880,  31926,  31970,  32014,  32056,  32097,
  32137,  32176,  32213,  32249,  32284,  32318,  32350,  32382,
  32412,  32441,  32468,  32495,  32520,  32544,  32567,  32588,
  32609,  32628,  32646,  32662,  32678,  32692,  32705,  32717,
  32727,  32736,  32744,  32751,  32757,  32761,  32764,  32766,
  32767,  32766,  32764,  32761,  32757,  32751,  32744,  32736,
  32727,  32717,  32705,  32692,  32678,  32662,  32646,  32628,
  32609,  32588,  32567,  32544,  32520,  32495,  32468,  32441,
  32412,  32382,  32350,  32318,  32284,  32249,  32213,  32176,
  32137,  32097,  32056,  32014,  31970,  31926,  31880,  31833,
  31785,  31735,  31684,  31633,  31580,  31525,  31470,  31413,
  31356,  31297,  31236,  31175,  31113,  31049,  30984,  30918,
  30851,  30783,  30713,  30643,  30571,  30498,  30424,  30349,
  30272,  30195,  30116,  30036,  29955,  29873,  29790,  29706,
  29621,  29534,  29446,  29358,  29268,  29177,  29085,  28992,
  28897,  28802,  28706,  28608,  28510,  28410,  28309,  28208,
  28105,  28001,  27896,  27790,  27683,  27575,  27466,  27355,
  27244,  27132,  27019,  26905,  26789,  26673,  26556,  26437,
  26318,  26198,  26077,  25954,  25831,  25707,  25582,  25456,
  25329,  25201,  25072,  24942,  24811,  24679,  24546,  24413,
  24278,  24143,  24006,  23869,  23731,  23592,  23452,  23311,
  23169,  23027,  22883,  22739,  22594,  22448,  22301,  22153,
  22004,  21855,  21705,  21554,  21402,  21249,  21096,  20942,
  20787,  20631,  20474,  20317,  20159,  20000,  19840,  19680,
  19519,  19357,  19194,  19031,  18867,  18702,  18537,  18371,
  18204,  18036,  17868,  17699,  17530,  17360,  17189,  17017,
  16845,  16672,  16499,  16325,  16150,  15975,  15799,  15623,
  15446,  15268,  15090,  14911,  14732,  14552,  14372,  14191,
  14009,  13827,  13645,  13462,  13278,  13094,  12909,  12724,
  12539,  12353,  12166,  11980,  11792,  11604,  11416,  11227,
  11038,  10849,  10659,  10469,  10278,  10087,   9895,   9703,
   9511,   9319,   9126,   8932,   8739,   8545,   8351,   8156,
   7961,   7766,   7571,   7375,   7179,   6982,   6786,   6589,
   6392,   6195,   5997,   5799,   5601,   5403,   5205,   5006,
   4807,   4608,   4409,   4210,   4011,   3811,   3611,   3411,
   3211,   3011,   2811,   2610,   2410,   2209,   2009,   1808,
   1607,   1406,   1206,   1005,    804,    603,    402,    201,
      0,   -201,   -402,   -603,   -804,  -1005,  -1206,  -1406,
  -1607,  -1808,  -2009,  -2209,  -2410,  -2610,  -2811,  -3011,
  -3211,  -3411,  -3611,  -3811,  -4011,  -4210,  -4409,  -4608,
  -4807,  -5006,  -5205,  -5403,  -5601,  -5799,  -5997,  -6195,
  -6392,  -6589,  -6786,  -6982,  -7179,  -7375,  -7571,  -7766,
  -7961,  -8156,  -8351,  -8545,  -8739,  -8932,  -9126,  -9319,
  -9511,  -9703,  -9895, -10087, -10278, -10469, -10659, -10849,
 -11038, -11227, -11416, -11604, -11792, -11980, -12166, -12353,
 -12539, -12724, -12909, -13094, -13278, -13462, -13645, -13827,
 -14009, -14191, -14372, -14552, -14732, -14911, -15090, -15268,
 -15446, -15623, -15799, -15975, -16150, -16325, -16499, -16672,
 -16845, -17017, -17189, -17360, -17530, -17699, -17868, -18036,
 -18204, -18371, -18537, -18702, -18867, -19031, -19194, -19357,
 -19519, -19680, -19840, -20000, -20159, -20317, -20474, -20631,
 -20787, -20942, -21096, -21249, -21402, -21554, -21705, -21855,
 -22004, -22153, -22301, -22448, -22594, -22739, -22883, -23027,
 -23169, -23311, -23452, -23592, -23731, -23869, -24006, -24143,
 -24278, -24413, -24546, -24679, -24811, -24942, -25072, -25201,
 -25329, -25456, -25582, -25707, -25831, -25954, -26077, -26198,
 -26318, -26437, -26556, -26673, -26789, -26905, -27019, -27132,
 -27244, -27355, -27466, -27575, -27683, -27790, -27896, -28001,
 -28105, -28208, -28309, -28410, -28510, -28608, -28706, -28802,
 -28897, -28992, -29085, -29177, -29268, -29358, -29446, -29534,
 -29621, -29706, -29790, -29873, -29955, -30036, -30116, -30195,
 -30272, -30349, -30424, -30498, -30571, -30643, -30713, -30783,
 -30851, -30918, -30984, -31049, -31113, -31175, -31236, -31297,
 -31356, -31413, -31470, -31525, -31580, -31633, -31684, -31735,
 -31785, -31833, -31880, -31926, -31970, -32014, -32056, -32097,
 -32137, -32176, -32213, -32249, -32284, -32318, -32350, -32382,
 -32412, -32441, -32468, -32495, -32520, -32544, -32567, -32588,
 -32609, -32628, -32646, -32662, -32678, -32692, -32705, -32717,
 -32727, -32736, -32744, -32751, -32757, -32761, -32764, -32766,
 -32767, -32766, -32764, -32761, -32757, -32751, -32744, -32736,
 -32727, -32717, -32705, -32692, -32678, -32662, -32646, -32628,
 -32609, -32588, -32567, -32544, -32520, -32495, -32468, -32441,
 -32412, -32382, -32350, -32318, -32284, -32249, -32213, -32176,
 -32137, -32097, -32056, -32014, -31970, -31926, -31880, -31833,
 -31785, -31735, -31684, -31633, -31580, -31525, -31470, -31413,
 -31356, -31297, -31236, -31175, -31113, -31049, -30984, -30918,
 -30851, -30783, -30713, -30643, -30571, -30498, -30424, -30349,
 -30272, -30195, -30116, -30036, -29955, -29873, -29790, -29706,
 -29621, -29534, -29446, -29358, -29268, -29177, -29085, -28992,
 -28897, -28802, -28706, -28608, -28510, -28410, -28309, -28208,
 -28105, -28001, -27896, -27790, -27683, -27575, -27466, -27355,
 -27244, -27132, -27019, -26905, -26789, -26673, -26556, -26437,
 -26318, -26198, -26077, -25954, -25831, -25707, -25582, -25456,
 -25329, -25201, -25072, -24942, -24811, -24679, -24546, -24413,
 -24278, -24143, -24006, -23869, -23731, -23592, -23452, -23311,
 -23169, -23027, -22883, -22739, -22594, -22448, -22301, -22153,
 -22004, -21855, -21705, -21554, -21402, -21249, -21096, -20942,
 -20787, -20631, -20474, -20317, -20159, -20000, -19840, -19680,
 -19519, -19357, -19194, -19031, -18867, -18702, -18537, -18371,
 -18204, -18036, -17868, -17699, -17530, -17360, -17189, -17017,
 -16845, -16672, -16499, -16325, -16150, -15975, -15799, -15623,
 -15446, -15268, -15090, -14911, -14732, -14552, -14372, -14191,
 -14009, -13827, -13645, -13462, -13278, -13094, -12909, -12724,
 -12539, -12353, -12166, -11980, -11792, -11604, -11416, -11227,
 -11038, -10849, -10659, -10469, -10278, -10087,  -9895,  -9703,
  -9511,  -9319,  -9126,  -8932,  -8739,  -8545,  -8351,  -8156,
  -7961,  -7766,  -7571,  -7375,  -7179,  -6982,  -6786,  -6589,
  -6392,  -6195,  -5997,  -5799,  -5601,  -5403,  -5205,  -5006,
  -4807,  -4608,  -4409,  -4210,  -4011,  -3811,  -3611,  -3411,
  -3211,  -3011,  -2811,  -2610,  -2410,  -2209,  -2009,  -1808,
  -1607,  -1406,  -1206,  -1005,   -804,   -603,   -402,   -201,
};












