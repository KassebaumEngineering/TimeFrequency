//
// Wigner.cc
//
// C++ Implementation of Wigner Transform Wigner Object
//
//  $Id: Wigner.cc,v 1.1 1994/10/06 17:51:59 jak Exp $
//
//  Author: John Kassebaum
//
// $Log: Wigner.cc,v $
// Revision 1.1  1994/10/06 17:51:59  jak
// Made Fixes to several of the spectrum programs - have a preliminary version
// of the Wigner distribution. -jak
//
//

static char rcsid_Wigner_cc[] = "$Id: Wigner.cc,v 1.1 1994/10/06 17:51:59 jak Exp $";

#include "Wigner.h"
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
//	double  *WignerT;
	
//	if( !WignerT ){
//	    if((WignerT = NewBlock( double *, getWindowSize() )) == 0 ){
//		    perror("Wigner::compute() - Can't Allocate memory For Wigner Transform Row Pointers.");
//			abort();
//		}
//	    if((WignerT[0] = NewBlock( double, getWindowSize() * getWindowSize() )) == 0 ){
//		    perror("Wigner::compute() - Can't Allocate memory For Wigner Transform Matrix.");
//			abort();
//		}
//		bzero( WignerT[0], getWindowSize() * getWindowSize());
//		for(i=1; i< getWindowSize(); i++){
//		    WignerT[i] = &( WignerT[0][ i * getWindowSize() ] );
//			WignerT[i][i] = 1.0;
//		}
//#ifdef DEBUG
//    fprintf(stderr, "Wigner Transform Matrix Computed!\n");
//#endif
//	}
	
	if( getSpectrumDataLength() == 0 ) {
	    setWindowStride( getWindowSize() / 2 );
	}
	
    if( !isAnalytic ) makeAnalytic();
	
#ifdef DEBUG
    fprintf(stderr, "Computing Wigner Transform of Analytic Signal\n");
#endif
	for( i=0; i< time_slots; i++){
		int tau, t1, t2, half, s_offset;
        Complex temp[200];
	  //
	  // Copy Signal to Spectrum time slot
	  //
//#ifdef DEBUG
//    fprintf(stderr, "C");
//#endif
	  //bcopy((char *)&(signal[ i*getWindowStride() ]), (char *)spectrogram[i], (getWindowSize() * sizeof( Complex )));
	  //
	  // Wigner Weighting of Signal
	  //
#ifdef DEBUG
    fprintf(stderr, "W");
#endif
//        {
//		    Complex temp[100][100];
//			int tau, t1, t2, half, s_offset;
			
			half = ( getWindowSize() >> 1 );
			s_offset = i*getWindowStride(); // start of window!

            t1 = i*getWindowStride() + half;
            for( tau = 0; (tau < getWindowSize()) && (t1 - tau >= 0) && (t1 + tau < signal_length) ; tau++ ){
                temp[ tau ] = conj( signal[ t1 - tau ] ) * signal[ t1 + tau ] ;
            }
            bcopy((char *)temp, (char *)spectrogram[i], (getWindowSize() * sizeof( Complex )));
//=====
//            for( tau = 0; (tau < half) ; tau += 1 ){
//			    spectrogram[i][tau*2] = Complex( 0.0, 0.0 );
//			    for( t1 = 0; (t1< half) && (s_offset + t1 + tau >= 0) && (s_offset + t1 - tau < signal_length); t1 += 1){
//                    spectrogram[i][tau*2] += (signal[s_offset + t1 + tau] * conj(signal[s_offset + t1 - tau])) ;
//                }
//			}
//			temp = Complex( 0.0, 0.0 );
//			for( t1 = 0; t1 < half; t1++ ){
//			    for( t2=0; t2 < half; t2++){
//				    temp[t1][t2] = ( signal[s_offset - t1 ] * conj( signal[s_offset + t2]) );
//				}
//			}
//=====
			
//		}
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

