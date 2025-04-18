//
// Wigner.cc
//
// C++ Implementation of Wigner Transform Wigner Object
//
//  $Id: Wigner.cc,v 1.3 1994/10/27 09:11:31 jak Exp $
//
//  Author: John Kassebaum
//
// $Log: Wigner.cc,v $
// Revision 1.3  1994/10/27 09:11:31  jak
// Fixes, including anti-aliasing additions. -jak
//
// Revision 1.2  1994/10/07  06:55:30  jak
// Wigner now works!  Bug fixes to the Spectrogram also.  Stride can now
// be set from the command line!  -jak
//
// Revision 1.1  1994/10/06  17:51:59  jak
// Made Fixes to several of the spectrum programs - have a preliminary version
// of the Wigner distribution. -jak
//
//

[[maybe_unused]] static char rcsid_Wigner_cc[] = "$Id: Wigner.cc,v 1.3 1994/10/27 09:11:31 jak Exp $";

#include "Wigner.h"
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <cstring> // For memcpy
using Complex = std::complex<double>;

#define DEBUG
#define POSITIVE_F_ONLY
#define ALIASFREE

#define BLOCKSIZE    1024
#define Null(A)             (static_cast<A *>(nullptr))
#define New(A)              (static_cast<A *>(malloc(sizeof(A))))
#define NewBlock(A,N)       (static_cast<A *>(malloc(sizeof(A) * (N))))
#define BiggerBlock(A,B,N)  (static_cast<A *>(realloc(static_cast<void *>(B), sizeof(A) * (N))))

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
    int i; // Removed unused j
	
	if( getSpectrumDataLength() == 0 ) {
	    setWindowStride( getWindowSize() / 2 );
	}
	
    if( !isAnalytic ) makeAnalytic();
	
#ifdef DEBUG
    fprintf(stderr, "Computing Wigner Transform of Analytic Signal\n");
#endif
	for( i=0; i< time_slots; i++){
		int tau, t1, half; // Removed unused t2
		      std::vector<Complex> temp(getWindowSize());

	  //
	  // Wigner Weighting of Signal
	  //
#ifdef DEBUG
    fprintf(stderr, "W");
#endif
		half = ( getWindowSize() >> 1 );
		// int s_offset = i*getWindowStride(); // start of window! (Unused)

		t1 = i*getWindowStride() + half;
		
#ifdef ALIASFREE
		for( tau = 0; (tau < getWindowSize()) && (t1 - tau/2 >= 0) && (t1 + tau/2 + 1 < signal_length) ; tau += 2 ){
			temp[ tau ] = std::conj( signal[ t1 - tau/2 ] ) * signal[ t1 + tau/2 ] ;
			temp[ tau+1 ] = std::conj( signal[ t1 - tau/2 ] ) * signal[ t1 + tau/2 + 1] ;
		}
#else
		for( tau = 0; (tau < getWindowSize()) && (t1 - tau >= 0) && (t1 + tau < signal_length) ; tau++ ){
			temp[ tau ] = std::conj( signal[ t1 - tau ] ) * signal[ t1 + tau ] ;
		}
#endif

		memcpy(spectrogram[i], temp.data(), getWindowSize() * sizeof(Complex));

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
    fprintf(stderr, "\nDone With Wigner Distributions!\n");
#endif
	
	isComputed = 1;
};

void Wigner:: print_Gnuplot()
{
  int i,j, status;
  float length, cnt, temp;

#ifdef POSITIVE_F_ONLY
  length = static_cast<float>(getWindowSize()) / 2.0;
#else
  length = static_cast<float>(getWindowSize()); // 2.0;
#endif // POSITIVE_F_ONLY

  if((status = fwrite((char*)&length, sizeof(float), 1, stdout))!=1){
    fprintf(stderr,"Spectrogram:: print_Gnuplot() - fwrite 1 returned %d\n",status);
    abort();
  }
  
#ifdef ALIASFREE
  for( cnt = 0.0; cnt < getBandWidth()/2.0; cnt += getFrequencyResolution()){
      fwrite((char*)&cnt, sizeof(float), 1, stdout);
  }
#else
  for( cnt = 0.0; cnt < getBandWidth()/4.0; cnt += getFrequencyResolution()/2.0){
      fwrite((char*)&cnt, sizeof(float), 1, stdout);
  }
#endif

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

