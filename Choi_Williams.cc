//
// Choi_Williams.cc
//
// C++ Implementation of Choi_Williams Transform Object
//
//  $Id: Choi_Williams.cc,v 1.1 1994/10/27 09:12:37 jak Exp $
//
//  Author: John Kassebaum
//
// $Log: Choi_Williams.cc,v $
// Revision 1.1  1994/10/27 09:12:37  jak
// Checking in first working anti-aliased Choi-Williams Distribution TFD -jak
//
//
//

static char rcsid_Choi_Williams_cc[] = "$Id: Choi_Williams.cc,v 1.1 1994/10/27 09:12:37 jak Exp $";

#include "Choi_Williams.h"
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

Choi_Williams:: Choi_Williams(): Wigner(), CW_Kernel(0), temp(0), sigma(10.0)
{
    ;
};

//
// Destructor
//

Choi_Williams:: ~Choi_Williams()
{
    if( temp[0] ) delete [] temp[0];
    if( temp ) free( temp );
    if( CW_Kernel[0] ) {
        free( CW_Kernel[0] );
        CW_Kernel[0] = Null( double );
    }
    if( CW_Kernel ) free( CW_Kernel );
};


void Choi_Williams::compute()
{
    int i,j;
    int NumberOfRows;
    int NumberOfCols;

    NumberOfRows = getWindowSize()/2 ;
	NumberOfCols = getWindowSize();

	if( !CW_Kernel ){
	    if((CW_Kernel = NewBlock( double *, NumberOfRows )) == 0 ){
		    perror("Choi_Williams::compute() - Can't Allocate memory For Choi_Williams Transform Row Pointers.");
			abort();
		}
	    if((CW_Kernel[0] = NewBlock( double, NumberOfRows * NumberOfCols)) == 0 ){
		    perror("Choi_Williams::compute() - Can't Allocate memory For Choi_Williams Transform Matrix.");
			abort();
		}
		bzero( (char *)CW_Kernel[0], NumberOfRows * NumberOfCols * sizeof(double)); //zero data
		for(i=1; i< NumberOfRows; i++){
		    CW_Kernel[i] = &( CW_Kernel[0][ i * NumberOfCols ] );
		}

      //--------------------
      // The Choi-Williams Kernel
      //
		for( j=0; j< NumberOfCols; j++){
            for( i=0; i< NumberOfRows; i++){
	            double variance, dist;
				
				dist = (double) (i - (NumberOfCols - 1 - j)/2 )*(i - (NumberOfCols - 1 - j)/2 );
				if (dist == 0.0) dist = 0.5;
				variance = 2.0 * (double) (( j - NumberOfCols/2 -1) * ( j - NumberOfCols/2 -1)) / sigma;
	            if( variance != 0.0 ) 
				    CW_Kernel[i][j] =  exp(-dist * variance);//erfc( dist / variance ); 
				else 
				    CW_Kernel[NumberOfRows/2][j] = 1.0;
                
                //if( i%2 ) {
                //    if( !(j%2) ) CW_Kernel[i][j] = 0.0;
                //} else {
                //    if( j%2 ) CW_Kernel[i][j] = 0.0;
                //}
			}
        }

//--------------------
// Test Case - The Wigner Kernel
//
//        for( j = 0, i = NumberOfRows-1; j < NumberOfCols-1; j += 2, i-- ){
//fprintf( stderr, "[%d,%d][%d,%d]",i,j,i,j+1);
//            CW_Kernel[i][j] = 1.0;
//            CW_Kernel[i][j + 1] = 1.0;
//        }
//--------------------

#ifdef DEBUG
/*
    fprintf(stderr, "Choi_Williams Transform Matrix Computed!\n");
	fprintf( stderr, "{\n");
	for( i=0; i< NumberOfRows; i++){
	    fprintf( stderr, "{ ", i);
		for( j=0; j< NumberOfCols; j++){
			fprintf( stderr, "%lf", CW_Kernel[i][j] ); 
			if( j < NumberOfCols - 1 ) fprintf( stderr, ", ");
		}
		fprintf( stderr, "}");
		if( i < NumberOfRows - 1) fprintf( stderr, ",");
		fprintf( stderr, "\n");
	}
	fprintf( stderr, "}\n");
*/
#endif
	}
	
	if( !temp ){
		if((temp = NewBlock( Complex*, NumberOfRows ) ) == 0 ){
			perror("Choi_Williams::compute() - Can't Allocate memory for temporary storage pointers.");
			abort();
		}
		if((temp[0] = new Complex[ NumberOfRows * NumberOfCols ] ) == 0 ){
			perror("Choi_Williams::compute() - Can't Allocate memory for temporary storage.");
			abort();
		}
		bzero( &(temp[0][0]), ( NumberOfRows * NumberOfCols ) * sizeof(Complex*) );
		for(i=1; i< NumberOfRows; i++){
			temp[i] = &( temp[0][ i * NumberOfCols ] );
		}
    }

	if( getSpectrumDataLength() == 0 ) {
	    setWindowStride( getWindowSize() / 2 );
	}
	
    if( !isAnalytic ) makeAnalytic();
	
#ifdef DEBUG
    fprintf(stderr, "Computing Choi_Williams Transform of Analytic Signal\n");
#endif
	for( i=0; i< time_slots; i++){
	    int tau, mu, t1, t2, t, half, s_offset;

		half = ( getWindowSize() >> 1 );
		s_offset = i*getWindowStride(); // start of window!
		t = i*getWindowStride() + half;

#ifdef DEBUG
    fprintf(stderr, "O");
#endif
      //
	  // Form the Current Outer Product Matrix in Tau-Column Form
	  //
        //for( mu = -half/2; mu < half/2 ; mu++ ){
        for( mu = 0; mu < half ; mu++ ){
            for( tau = 0; (tau < (half<<1)-1) ; tau+=2 ){
                register int tau_over_2 = tau >> 1;
                register int mu_plus_t = mu + s_offset + half/2;
                register int t_minus_mu =  t - mu;
                //register int row = mu + half/2;
                register int row = mu;
                if( mu_plus_t - tau_over_2 < 0 ) {
                //if( t_minus_mu - tau_over_2 < 0 ) {
                    temp[row][tau]   = 0;
                    temp[row][tau+1] = 0;
                } else if ( mu_plus_t + tau_over_2 + 1 >= signal_length) {
                    temp[row][tau+1] = 0;
                    if ( mu_plus_t + tau_over_2 >= signal_length)
                        temp[row][tau] = 0;
                    else
                        temp[row][tau] = conj( signal[ mu_plus_t - tau_over_2] ) * signal[ mu_plus_t + tau_over_2];
                        //temp[row][tau] = conj( signal[ t_minus_mu - tau_over_2] ) * signal[ mu_plus_t + tau_over_2];
                } else {
                    temp[row][tau]   = conj( signal[ mu_plus_t - tau_over_2] ) * signal[ mu_plus_t + tau_over_2];
                    temp[row][tau+1] = conj( signal[ mu_plus_t - tau_over_2] ) * signal[ mu_plus_t + tau_over_2 + 1];
                    //temp[row][tau]   = conj( signal[ t_minus_mu - tau_over_2] ) * signal[ mu_plus_t + tau_over_2];
                    //temp[row][tau+1] = conj( signal[ t_minus_mu - tau_over_2] ) * signal[ mu_plus_t + tau_over_2 + 1];
                }
            }
        }

#ifdef DEBUG
    fprintf(stderr, "W");
#endif
			
      //
	  // Weight the Outer Product Matrix With the Choi-Williams Convolution Kernel
	  // and sum along the columns to produce the Pre-FFT data vectors.
      //
        for( tau = 0; tau < NumberOfCols; tau++ ){
            spectrogram[i][tau] = 0.0;
            for( mu = 0; mu < NumberOfRows ; mu++ ){
                spectrogram[i][tau] += ( temp[ mu ][ tau ] * CW_Kernel[ mu ][ tau ] );
            }
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
    fprintf(stderr, "\nDone With Choi_Williams Distribution!\n");
#endif
	
	isComputed = 1;
};

