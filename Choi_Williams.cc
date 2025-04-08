//
// Choi_Williams.cc
//
// C++ Implementation of Choi_Williams Transform Object
//
//  $Id: Choi_Williams.cc,v 1.2 1994/11/18 05:52:47 jak Exp $
//
//  Author: John Kassebaum
//
// $Log: Choi_Williams.cc,v $
// Revision 1.2  1994/11/18 05:52:47  jak
// Small changes to improve Choi_Williams operations. -jak
//
// Revision 1.1  1994/10/27  09:12:37  jak
// Checking in first working anti-aliased Choi-Williams Distribution TFD -jak
//
//
//

static char rcsid_Choi_Williams_cc[] = "$Id: Choi_Williams.cc,v 1.2 1994/11/18 05:52:47 jak Exp $";

#include "Choi_Williams.h"
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstdio>
using Complex = std::complex<double>;


#define DEBUG

#define BLOCKSIZE    1024
#define Null(A)             (static_cast<A *>(nullptr))
#define New(A)              (static_cast<A *>(malloc(sizeof(A))))
#define NewBlock(A,N)       (static_cast<A *>(malloc(sizeof(A) * (N))))
#define BiggerBlock(A,B,N)  (static_cast<A *>(realloc(static_cast<void *>(B), sizeof(A) * (N))))

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
				    CW_Kernel[i][j] =  exp(-dist * variance); 
				else 
				    CW_Kernel[NumberOfRows/2][j] = 1.0;
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
                int tau_over_2 = tau >> 1;
                int mu_plus_t = mu + s_offset + half/2;
                int t_minus_mu =  t - mu;
                //int row = mu + half/2;
                int row = mu;
                if( mu_plus_t - tau_over_2 < 0 ) {
                //if( t_minus_mu - tau_over_2 < 0 ) {
                    temp[row][tau]   = 0;
                    temp[row][tau+1] = 0;
                } else if ( mu_plus_t + tau_over_2 + 1 >= signal_length) {
                    temp[row][tau+1] = 0;
                    if ( mu_plus_t + tau_over_2 >= signal_length)
                        temp[row][tau] = 0;
                    else
                        temp[row][tau] = std::conj( signal[ mu_plus_t - tau_over_2] ) * signal[ mu_plus_t + tau_over_2];
                        //temp[row][tau] = conj( signal[ t_minus_mu - tau_over_2] ) * signal[ mu_plus_t + tau_over_2];
                } else {
                    temp[row][tau]   = std::conj( signal[ mu_plus_t - tau_over_2] ) * signal[ mu_plus_t + tau_over_2];
                    temp[row][tau+1] = std::conj( signal[ mu_plus_t - tau_over_2] ) * signal[ mu_plus_t + tau_over_2 + 1];
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

