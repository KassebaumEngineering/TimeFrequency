//
// Spectrogram.cc
//
// C++ Implementation of Short Time Fourier Transform Spectrogram Object
//
//  $Id: Spectrogram.cc,v 1.1 1994/10/04 07:21:04 jak Exp $
//
//  Author: John Kassebaum
//
/* $Log: Spectrogram.cc,v $
/* Revision 1.1  1994/10/04 07:21:04  jak
/* Initial revision
/**/

static char rcsid_Spectrogram_cc[] = "$Id: Spectrogram.cc,v 1.1 1994/10/04 07:21:04 jak Exp $";

#include "Spectrogram.h"
#include <math.h>
#include <Complex.h>
#include <stdlib.h>
#include <stdio.h>

#define DEBUG

#define BLOCKSIZE    1024
#define Null(A)             ((A *) 0)
#define New(A)              ((A *) malloc( sizeof(A) ) )
#define NewBlock(A,N)       ((A *) malloc( sizeof(A) * (N)) )
#define BiggerBlock(A,B,N)  ((A *) realloc( (void *)(B), sizeof(A) * (N)))

//
// Constructor
//

Spectrogram:: Spectrogram(): TimeFrequency(), spectrogram(0), time_slots(0),  isComputed(0)
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
    int spectrogram_size, i; 
	double ratio;
	
	TimeFrequency::setWindowStride( astride );

	ratio = getWindowSize() / getWindowStride();
	time_slots = (int) floor( ( (double)signal_length / (double)getWindowStride() ) - ratio + 1.0);
	
    if( spectrogram ) {
	    if( (spectrogram = BiggerBlock( Complex *, spectrogram, time_slots )) == Null( Complex )){
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
	    if( (spectrogram = NewBlock( Complex *, time_slots )) == Null( Complex )){
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
	    spectrogram[i] = &( spectrogram[0][ i*getWindowStride() ] );
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
			hamming = 0.54 + 0.46*cos( M_PI * (double)j / (double)(half) );
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
  register int i,j, status;
  float length, cnt, temp;

  length = (float) getWindowSize() / 2.0;

  if((status = fwrite((char*)&length, sizeof(float), 1, stdout))!=1){
    fprintf(stderr,"Spectrogram:: print_Gnuplot() - fwrite 1 returned %d\n",status);
    abort();
  }
  
  for( cnt = getBandWidth()/2.0 - getFrequencyResolution(); cnt >= 0.0; cnt -= getFrequencyResolution()){
      fwrite((char*)&cnt, sizeof(float), 1, stdout);
  }
    
  for(j=0, cnt=0.0; j < time_slots; j++, cnt += (getWindowStride() / getBandWidth())) {
      fwrite((char*)(&cnt), sizeof(float), 1, stdout);
      for(i=getWindowSize()/2; i< getWindowSize(); i++) {
          temp = norm( spectrogram[j][i] );
          fwrite((char*)(&temp), sizeof(float), 1, stdout);
      }
  }
};

