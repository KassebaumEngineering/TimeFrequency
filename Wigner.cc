//
// Wigner.cc
//
// C++ Implementation of Wigner Transform Wigner Object
//
//  $Id: Wigner.cc,v 1.2 1994/10/07 06:55:30 jak Exp $
//
//  Author: John Kassebaum
//
// $Log: Wigner.cc,v $
// Revision 1.2  1994/10/07 06:55:30  jak
// Wigner now works!  Bug fixes to the Spectrogram also.  Stride can now
// be set from the command line!  -jak
//
// Revision 1.1  1994/10/06  17:51:59  jak
// Made Fixes to several of the spectrum programs - have a preliminary version
// of the Wigner distribution. -jak
//
//

static char rcsid_Wigner_cc[] = "$Id: Wigner.cc,v 1.2 1994/10/07 06:55:30 jak Exp $";

#include "Wigner.h"
#include <math.h>
#include <Complex.h>
#include <stdlib.h>
#include <stdio.h>

#define DEBUG
#define POSITIVE_F_ONLY

#define BLOCKSIZE    1024
#define Null(A)             ((A *) 0)
#define New(A)              ((A *) malloc( sizeof(A) ) )
#define NewBlock(A,N)       ((A *) malloc( sizeof(A) * (N)) )
#define BiggerBlock(A,B,N)  ((A *) realloc( (void *)(B), sizeof(A) * (N)))

//
// Constructor
//

Wigner:: Wigner(): Spectrogram()
{
    ;
};

//
// Destructor
//

Wigner:: ~Wigner()
{
};


void Wigner::compute()
{
    int i,j;
	
	if( getSpectrumDataLength() == 0 ) {
	    setWindowStride( getWindowSize() / 2 );
	}
	
    if( !isAnalytic ) makeAnalytic();
	
#ifdef DEBUG
    fprintf(stderr, "Computing Wigner Transform of Analytic Signal\n");
#endif
	for( i=0; i< time_slots; i++){
		int tau, t1, t2, half, s_offset;
        Complex temp[ getWindowSize() ];

	  //
	  // Wigner Weighting of Signal
	  //
#ifdef DEBUG
    fprintf(stderr, "W");
#endif
		half = ( getWindowSize() >> 1 );
		s_offset = i*getWindowStride(); // start of window!

		t1 = i*getWindowStride() + half;
		for( tau = 0; (tau < getWindowSize()) && (t1 - tau >= 0) && (t1 + tau < signal_length) ; tau++ ){
			temp[ tau ] = conj( signal[ t1 - tau ] ) * signal[ t1 + tau ] ;
		}
		bcopy((char *)temp, (char *)spectrogram[i], (getWindowSize() * sizeof( Complex )));

	  //
	  // Perform FFT
	  //
#ifdef DEBUG
    fprintf(stderr, "F");
#endif
		fft( spectrogram[i], getWindowSize() );
#ifdef DEBUG
    fprintf(stderr, "-");
#endif
	}
#ifdef DEBUG
    fprintf(stderr, "\nDone With STFTs!\n");
#endif
	
	isComputed = 1;
};

void Wigner:: print_Gnuplot()
{
  register int i,j, status;
  float length, cnt, temp;

#ifdef POSITIVE_F_ONLY
  length = (float) getWindowSize() / 2.0;
#else
  length = (float) getWindowSize(); // 2.0;
#endif

  if((status = fwrite((char*)&length, sizeof(float), 1, stdout))!=1){
    fprintf(stderr,"Spectrogram:: print_Gnuplot() - fwrite 1 returned %d\n",status);
    abort();
  }
  
  for( cnt = 0.0; cnt < getBandWidth()/4.0; cnt += getFrequencyResolution()/2.0){
      fwrite((char*)&cnt, sizeof(float), 1, stdout);
  }

#ifndef POSITIVE_F_ONLY
  for( cnt = -getBandWidth()/4.0; cnt < 0.0; cnt += getFrequencyResolution()/2.0){
      fwrite((char*)&cnt, sizeof(float), 1, stdout);
  }
#endif
    
  for(j=0, cnt=(getWindowSize()/getBandWidth()); j < time_slots; j++, cnt += (getWindowStride() / getBandWidth())) {
      fwrite((char*)(&cnt), sizeof(float), 1, stdout);
#ifdef POSITIVE_F_ONLY
      for(i=0; i< getWindowSize()/2; i++) {
#else
      for(i=0; i< getWindowSize(); i++) {
#endif
          temp = norm( spectrogram[j][i] );
          fwrite((char*)(&temp), sizeof(float), 1, stdout);
      }
  }
};

