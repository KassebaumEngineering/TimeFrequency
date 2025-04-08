//
// TimeFrequency.cc
//
// C++ Implementation of TimeFrequency Object
//
//  $Id: TimeFrequency.cc,v 1.2 1994/10/07 06:55:29 jak Exp $
//
//  Author: John Kassebaum
//
/* $Log: TimeFrequency.cc,v $
 * Revision 1.2  1994/10/07 06:55:29  jak
 * Wigner now works!  Bug fixes to the Spectrogram also.  Stride can now
 * be set from the command line!  -jak
 *
// Revision 1.1.1.1  1994/10/04  07:21:05  jak
// Placing Time/Frequency Code under CVS control.  Only Spectrogram
// works currently.  -jak
//*/

[[maybe_unused]] static char rcsid_TimeFrequency_cc[] = "$Id: TimeFrequency.cc,v 1.2 1994/10/07 06:55:29 jak Exp $";

#include "TimeFrequency.h"
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>

using Complex = std::complex<double>;

#define DEBUG

#define DEFAULT_WSIZE  256   // Must be a power of 2
#define HALF_HILBERT   34    // half-length of Hilbert transform filter 
#define HILBERT_SIZE   2*HALF_HILBERT + 1   // Must be odd

#define BLOCKSIZE    1024
#define Null(A)             (static_cast<A *>(nullptr))
#define New(A)              (static_cast<A *>(malloc(sizeof(A))))
#define NewBlock(A,N)       (static_cast<A *>(malloc(sizeof(A) * (N))))
#define BiggerBlock(A,B,N)  (static_cast<A *>(realloc(static_cast<void *>(B), sizeof(A) * (N))))

//
// Constructor
//

TimeFrequency:: TimeFrequency(): myWindow( Rectangular ), window_size( DEFAULT_WSIZE ), stride( DEFAULT_WSIZE/2), sampling_frequency(1.0), signal(0), isAnalytic(0)
{
	frequency_band       = (sampling_frequency);              // ( in Hz ) see Nyquist
	   frequency_resolution = (frequency_band / static_cast<double>(window_size));  // ( in Hz per sample )
};

//
// Destructor
//

TimeFrequency:: ~TimeFrequency()
{
    if( signal ) free( signal );
};

void TimeFrequency::setDataSignal( float *data, int n_items )
{
    int i;
	
    if( (signal = NewBlock( Complex, n_items )) == Null( Complex )){
	     perror("TimeFrequency::setDataSignal() -> Can't alloc space for signal");
		 abort();
	}
	for(i=0; i< n_items; i++){
	    signal[i] = Complex( data[i], 0.0 );
	}
	signal_length = n_items;
	isAnalytic = 0;
};

void TimeFrequency::setDataSignal( double *data, int n_items )
{
    int i;
	
    if( (signal = NewBlock( Complex, n_items )) == Null( Complex )){
	     perror("TimeFrequency::setDataSignal() -> Can't alloc space for signal");
		 abort();
	}
	for(i=0; i< n_items; i++){
	    signal[i] = Complex( data[i], 0.0 );
	}
	signal_length = n_items;
	isAnalytic = 0;
};


int TimeFrequency::setWindowSize( unsigned short asize )
{
	short mask, bitcount, the_size = 0;
	   int cnt, rtn;

	mask = asize;
    bitcount = 0;
	cnt = 0;
	while (mask){
		if (mask & 0x01) bitcount++;
		mask = (mask >> 1);
		cnt++;
	}
	if( bitcount != 1 ){
		//
		// wsize must be a power of 2 in number of samples.
		//
		if ((1 << cnt) <= 0xffff ){
		    the_size = (1 << cnt);
			rtn = the_size;
		} else {
		    rtn = 0;
		}
	} else {
		the_size = asize;
		rtn = the_size;
    }	    

    window_size          = the_size;
    frequency_resolution = (frequency_band / static_cast<double>(window_size));  // ( in Hz per sample )
	
	return rtn;;
};

void TimeFrequency::setWindowStride( unsigned short asize )
{
	short mask, bitcount;
	
	mask = asize;
    bitcount = 0;
	while (mask){
		if (mask & 0x01) bitcount++;
		mask = (mask >> 1);
	}
	if( bitcount != 1 ){
		fprintf(stderr,"bitcount == %d , asize == %X\n", bitcount, asize);
		fprintf(stderr,"TimeFrequency::setWindowStride(): Stride of length %d is not a power of 2!\n", asize);
		abort();
	}

    stride = asize;
};

void TimeFrequency::setSamplingFrequency( double afreq )
{ 
	sampling_frequency   = afreq;
	frequency_band       = (sampling_frequency);              // ( in Hz ) see Nyquist
	   frequency_resolution = (frequency_band / static_cast<double>(window_size));  // ( in Hz per sample )
};

void TimeFrequency::makeAnalytic()
{
    int k, j, i;
	static int initialized = 0;
	static double hilbert_h[ HILBERT_SIZE ];

  //
  // Creat Hilbert Transform filter Using a Hamming Window.
  //
    if (!initialized){
	    double hamming;
		hilbert_h[HALF_HILBERT] = 0.0;
		for (i=1; i<=HALF_HILBERT; i++) {
			hamming = 0.54 + 0.46*cos( M_PI * static_cast<double>(i) / static_cast<double>(HALF_HILBERT) );
			hilbert_h[HALF_HILBERT-i] = hamming * ( -static_cast<double>(i%2) * 2.0 / ( M_PI * static_cast<double>(i) ) );
			hilbert_h[HALF_HILBERT+i] = -hilbert_h[HALF_HILBERT-i];
		}
		initialized = 1;
#ifdef DEBUG
    fprintf(stderr, "Hilbert Transform Filter Built!\n");
#endif
	}
//	isAnalytic = 1;
//  return;

#ifdef DEBUG
    fprintf(stderr, "Convolving Signal with Hilbert Filter\n");
#endif

  //
  // Convolve Hilbert Transform with Signal and store in imaginary part of signal.
  //
  // Analytic Signal = Signal  +  ( I * Hilbert Transform of Signal )
  //
	for( i=0 ; i< signal_length; i++ ){
	    double s_hilbert;
		int start_H, stop_H;
		start_H = ((i - HALF_HILBERT) <= 0) ? (HALF_HILBERT - i) : 0;
		stop_H = ((signal_length - i) >= HALF_HILBERT) ? HILBERT_SIZE : (HALF_HILBERT + (signal_length - i));
		k = ( (i-HALF_HILBERT) <= 0 ) ? 0 : i-HALF_HILBERT ;
	    for( s_hilbert = 0.0, j = start_H; j < stop_H; j++ , k++ ){
		    s_hilbert += ( hilbert_h[j] * signal[k].real() );
		}
		signal[i] += Complex( 0.0, s_hilbert );
#ifdef DEBUG
    fprintf(stderr, ".");
#endif
	}
	
#ifdef DEBUG
    fprintf(stderr, "\nDone! - Signal is Analytic!\n");
#endif
    isAnalytic = 1;
};

int    TimeFrequency:: setFrequencyResolution( double delta_freq )
{
    double wsize;
	int i, flag, rtn;;
	
    wsize = ( frequency_band / delta_freq );
  //
  // wsize must be a power of 2 in number of samples.
  //
    for( flag = 0, i = 2; i < 0xffff; i = i<<1 ){
	    if( static_cast<double>(i) >= wsize ){
	     wsize = static_cast<double>(i);
		    flag = 1;
		    break;
		}
	}
	
	if( flag ) {
		window_size = static_cast<unsigned short>(wsize);
		frequency_resolution = delta_freq;  // ( in Hz per sample )
		rtn = window_size;
	} else {
	    fprintf(stderr, "TimeFrequency::setFrequencyResolution( %lf ) -> Result in Window Size %lf greater than max %d samples!\n", delta_freq, wsize, 0xffff );
		rtn = 0;
	}
	
	return rtn;
};

double  TimeFrequency:: getFrequencyResolution()
{
    return  frequency_resolution;
};

double  TimeFrequency:: getBandWidth()
{
    return  frequency_band;
};

