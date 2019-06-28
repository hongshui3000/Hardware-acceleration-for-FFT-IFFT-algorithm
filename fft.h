/*
 * fft.h
 *
 *  Created on: 08.03.2019
 *      Author: yukai.shen
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <xtensa/tie/fft_tie.h>

#ifndef FFT_H
#define FFT_H

#define PI 3.1415926535
#define POWER 4
#define FFTSIZE (1<<POWER)
#define N_WAVE 1024
#define DO_INVERSE 0
#define FFT 4 //0:DIT(); 1:DIF(); 2:DIT_macTIE(); 3:DIF_macTIE(); 4:DIT_simdTIE(); 5:DIF_simdTIE();

void DIT();
void DIF();
void IFFT();
void DIT_macTIE();
void DIF_macTIE();
void DIT_simdTIE();
void DIF_simdTIE();
void Bitreverse();
void Remap_1(short);
void Remap_2(short);
void GetTwiddleFactors(short);
void Output();




#endif /* FFT_H  */
