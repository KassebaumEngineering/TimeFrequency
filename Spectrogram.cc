//
// Spectrogram.cc
//
// C++ Implementation of Short Time Fourier Transform Spectrogram Object
//
//  $Id: Spectrogram.cc,v 1.3 1994/10/07 06:55:28 jak Exp $
//
//  Author: John Kassebaum
//
/* $Log: Spectrogram.cc,v $
 * Revision 1.3  1994/10/07 06:55:28  jak
 * Wigner now works!  Bug fixes to the Spectrogram also.  Stride can now
 * be set from the command line!  -jak
 *
// Revision 1.2  1994/10/06  17:51:55  jak
// Made Fixes to several of the spectrum programs - have a preliminary version
// of the Wigner distribution. -jak
//
// Revision 1.1.1.1  1994/10/04  07:21:05  jak
// Placing Time/Frequency Code under CVS control.  Only Spectrogram
// works currently.  -jak
//*/

[[maybe_unused]] static char rcsid_Spectrogram_cc[] = "$Id: Spectrogram.cc,v 1.3 1994/10/07 06:55:28 jak Exp $";

#include "Spectrogram.h"
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstdio>

using Complex = std::complex<double>;

#define DEBUG
#define POSITIVE_F_ONLY

#define BLOCKSIZE    1024
#define Null(A)             (static_cast<A *>(nullptr))
#define New(A)              (static_cast<A *>(malloc(sizeof(A))))
#define NewBlock(A,N)       (static_cast<A *>(malloc(sizeof(A) * (N))))
#define BiggerBlock(A,B,N)  (static_cast<A *>(realloc(static_cast<void *>(B), sizeof(A) * (N))))

//
// Constructor
//

Spectrogram:: Spectrogram(): TimeFrequency(), time_slots(0), spectrogram(0), isComputed(0)
{
    ;
};

//
// Destructor
//

Spectrogram:: ~Spectrogram()
{
    if( spectrogram ) free( spectrogram );
};

void Spectrogram:: setWindowStride( unsigned short astride )
{
    int i;
	double ratio;
	
	TimeFrequency::setWindowStride( astride );

	ratio = getWindowSize() / getWindowStride();
	time_slots = static_cast<int>(floor(((static_cast<double>(signal_length) / static_cast<double>(getWindowStride())) - ratio + 1.0)));
	
    if( spectrogram ) {
	    if( (spectrogram = BiggerBlock( Complex *, spectrogram, time_slots )) == Null( Complex * )){
	     perror("Spectrogram:: setWindowStride - Can't reallocate Memory for spectrogram ");
	  abort();
		}
	    if( (spectrogram[0] = BiggerBlock( Complex, spectrogram[0], (time_slots * getWindowSize()) )) == Null( Complex )){
		    perror("Spectrogram:: setWindowStride - Can't reallocate Memory for spectrogram data");
			abort();
		}
		spectrum_data_length = (time_slots * getWindowSize());
		bzero( (char *)spectrogram[0], spectrum_data_length * sizeof( Complex ) );
	} else {
	    if( (spectrogram = NewBlock( Complex *, time_slots )) == Null( Complex * )){
	     perror("Spectrogram:: setWindowStride - Can't allocate Memory for spectrogram ");
	  abort();
		}
	    if( (spectrogram[0] = NewBlock( Complex, (time_slots * getWindowSize()) )) == Null( Complex )){
		    perror("Spectrogram:: setWindowStride - Can't allocate Memory for spectrogram data");
			abort();
		}
		spectrum_data_length = (time_slots * getWindowSize());
		bzero( (char *)spectrogram[0], spectrum_data_length * sizeof( Complex ) );
	}
	
	for( i=0; i< time_slots; i++){
	    spectrogram[i] = &( spectrogram[0][ i*getWindowSize() ] );
	}
	isComputed = 0;
};

void Spectrogram::compute()
{
    int i,j, half;
	
	if( !spectrogram ) {
	    setWindowStride( getWindowSize() / 2 );
	}
	
    if( !isAnalytic ) makeAnalytic();
	
#ifdef DEBUG
    fprintf(stderr, "Computing STFT of Analytic Signal\n");
#endif
	for( i=0; i< time_slots; i++){
	  //
	  // Copy Signal to Spectrogram time slot
	  //
#ifdef DEBUG
    fprintf(stderr, "C");
#endif
	    bcopy((char *)&(signal[ i*getWindowStride() ]), (char *)spectrogram[i], (getWindowSize() * sizeof( Complex )));
	  //
	  // Window Signal
	  //
#ifdef DEBUG
    fprintf(stderr, "W");
#endif
	    half = ( getWindowSize() >> 1 );
	    for(j=0; j < half; j++){
		    double hamming;
			hamming = 0.54 + 0.46*cos( M_PI * static_cast<double>(j) / static_cast<double>(half) );
			spectrogram[i][half + j] *= hamming;
			spectrogram[i][half - j - 1] *= hamming;
		}
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

void Spectrogram:: print_Mathematica()
{
    int x,y;

    fprintf(stdout, "\n{");
    for (y=0; y < time_slots; y++) {
      fprintf(stdout,"{");
      for (x=0; x< getWindowSize()-1; x++) {
          fprintf(stdout,"%lf,", norm( spectrogram[y][x] ) );
      }
      /* last number */
	  if( y == time_slots-1 ){
          fprintf(stdout,"%lf}", norm( spectrogram[y][x] ) );
	  } else {
          fprintf(stdout,"%lf},\n", norm( spectrogram[y][x] ) );
	  }
    }
    fprintf(stdout,"}\n");
};

void Spectrogram:: print_Gnuplot()
{
  int i,j, status;
  float length, cnt, temp;

#ifdef POSITIVE_F_ONLY
  length = static_cast<float>(getWindowSize()) / 2.0;
#else
  length = static_cast<float>(getWindowSize()); // 2.0;
#endif

  if((status = fwrite((char*)&length, sizeof(float), 1, stdout))!=1){
    fprintf(stderr,"Spectrogram:: print_Gnuplot() - fwrite 1 returned %d\n",status);
    abort();
  }
  
  for( cnt = 0.0; cnt < getBandWidth()/2.0; cnt += getFrequencyResolution()){
      fwrite((char*)&cnt, sizeof(float), 1, stdout);
  }

#ifndef POSITIVE_F_ONLY
  for( cnt = -getBandWidth()/2.0; cnt < 0.0; cnt += getFrequencyResolution()){
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

