/*
 * main.c
 *
 *  Created on: 08.03.2019
 *      Author: yukai.shen
 */


#include "fft.h"

short SampleArray_Real[FFTSIZE];
short SampleArray_Imag[FFTSIZE];
short W_Real[FFTSIZE/2];
short W_Imag[FFTSIZE/2];

int main()
{
	int i;
	FILE *inpFile;
	FILE *outFile;

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 														Input Generate														  	//
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	for(i=0; i<FFTSIZE; i++){
		SampleArray_Imag[i] = 0;
		SampleArray_Real[i] = 1000*cos(2*PI*i/FFTSIZE);
	}


	inpFile = fopen("input", "w");
	if(inpFile != NULL){
		//Write DIT Results in a file
		for(i=0; i<FFTSIZE; i++)
			fprintf(inpFile,"Input[%d]: real:%d\timag:%d\n", i, SampleArray_Real[i], SampleArray_Imag[i]);
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 														FFT without TIE														  	//
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//====
// DIT
//====
#if FFT==0
	DIT();
	outFile = fopen("outputDIT", "w");
	if(outFile != NULL){
		//Write DIT Results in a file
		for(i=0; i<FFTSIZE; i++)
			fprintf(outFile,"Output[%d]: real:%d\timag:%d\n", i, SampleArray_Real[i], SampleArray_Imag[i]);
	}else
	{
		printf("\nERROR: Couldn't create file for outputDIT!\n");
		return 1;
	}
#endif

//====
// DIF
//====
#if FFT==1
	DIF();
	outFile = fopen("outputDIF", "w");
	if(outFile != NULL){
		for(i=0; i<FFTSIZE; i++)
			fprintf(outFile,"Output[%d]: real:%d\timag:%d\n", i, SampleArray_Real[i], SampleArray_Imag[i]);
	}else
	{
		printf("\nERROR: Couldn't create file for outputDIF!\n");
		return 1;
	}
#endif

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 														FFT with MAC														  	//
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//===========
// DIT_macTIE
//===========
#if FFT==2
	DIT_macTIE();
	outFile = fopen("outputDIT_macTie", "w");
	if(outFile != NULL){
		for(i=0; i<FFTSIZE; i++)
			fprintf(outFile,"Output[%d]: real:%d\timag:%d\n", i, SampleArray_Real[i], SampleArray_Imag[i]);
	}else
	{
		printf("\nERROR: Couldn't create file for outputDIT_macTIE!\n");
		return 1;
	}
#endif

//===========
// DIF_macTIE
//===========
#if FFT==3
	DIF_macTIE();
	outFile = fopen("outputDIF_macTie", "w");
	if(outFile != NULL){
		for(i=0; i<FFTSIZE; i++)
			fprintf(outFile,"Output[%d]: real:%d\timag:%d\n", i, SampleArray_Real[i], SampleArray_Imag[i]);
	}else
	{
		printf("\nERROR: Couldn't create file for outputDIF_macTIE!\n");
		return 1;
	}
#endif
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 														FFT with SIMD														  	//
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//============
// DIT_simdTIE
//============
#if FFT==4
	DIT_simdTIE();
	outFile = fopen("outputDIT_simdTIE", "w");
	if(outFile != NULL){
		for(i=0; i<FFTSIZE; i++)
			fprintf(outFile,"Output[%d]: real:%d\timag:%d\n", i, SampleArray_Real[i], SampleArray_Imag[i]);
	}else
		{
			printf("\nERROR: Couldn't create file for outputDIT_simdTIE!\n");
			return 1;
		}
#endif

//============
// DIF_simdTIE
//============
#if FFT==5
	DIF_simdTIE();
	outFile = fopen("outputDIF_simdTIE", "w");
	if(outFile != NULL){
		for(i=0; i<FFTSIZE; i++)
			fprintf(outFile,"Output[%d]: real:%d\timag:%d\n", i, SampleArray_Real[i], SampleArray_Imag[i]);
	}else
		{
			printf("\nERROR: Couldn't create file for outputDIF_simdTIE!\n");
			return 1;
		}
#endif

	return 0;
}
