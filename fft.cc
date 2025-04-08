//
//
//  fft.cc  a simple radix 2, in-place fft, and inv_fft
//  also includes simple data windows.
//
//  $Id: fft.cc,v 1.1 1994/10/04 07:21:04 jak Exp $
//
//
// History:
/* $Log: fft.cc,v $
 * Revision 1.1  1994/10/04 07:21:04  jak
 * Initial revision
 */
 
[[maybe_unused]] static char rcsid_fft_cc[] = "$Id: fft.cc,v 1.1 1994/10/04 07:21:04 jak Exp $";
 
#include <cmath>
#include <complex>
#include <cstdio>

using Complex = std::complex<double>;

#define DEBUG


void fft( Complex *a, int N ) // ----> N MUST be a power of 2 !
{
    int i, j, k, le, le1;
    int nv2, ip;
	Complex t, u, w;
	
#ifdef DEBUG
    {
	    int mask, bitcount;
		
		mask = N;
        bitcount = 0;
		while (mask){
			if (mask & 0x01) bitcount++;
		    mask = (mask >> 1);
		}
		if( bitcount != 1 ){
		    fprintf(stderr,"bitcount == %d , N == %X\n", bitcount, N);
		    fprintf(stderr,"fft(): sequence of length %d is not a power of 2!\n",N);
			abort();
		}
	}
#endif

  //
  // Place in Bit Reverse Order
  //
	nv2 = N >> 1; // N / 2
	for( j=1, i=1; ( i < N ); i++, j += k ){
	    if(i < j) {
		    t    = a[j];
			a[j] = a[i];
			a[i] = t;
		} 
        k = nv2;
		while ((k < j) && (k > 1)) {
		    j -= k;
            k =  k >> 1;
		}
	}

  //
  // Decimation in time, radix 2 fft
  //
	for( le = 0x02; le<= N; le = le << 1) {
		le1 = le >> 1;   // = le / 2 ;
		u = Complex( 1.0, 0.0 );
		w = Complex( cos( M_PI / static_cast<double>(le1) ), sin( M_PI / static_cast<double>(le1) ) );
		for( j = 0; j < le1; j++ ){
			for( i = j; i < N; i += le ){
				ip = i + le1;
				t = a[ ip ] * u;
				a[ ip ] = a[ i ] - t;
				a[ i ]  = a[ i ] + t;
			}
			u = u * w;
		}
	}
}


void inv_fft(Complex *a, int N) // ----> N MUST be a power of 2 !
{
    int i;

#ifdef DEBUG
    {
	    int mask, bitcount;
		
		mask = N;
        bitcount = 0;
		while (mask){
			if (mask & 0x01) bitcount++;
		    mask = (mask >> 1);
		}
		if( bitcount != 1 ){
		    fprintf(stderr,"bitcount == %d , N == %X\n", bitcount, N);
		    fprintf(stderr,"inv_fft(): sequence of length %d is not a power of 2!\n",N);
			abort();
		}
	}
#endif
//
// inv_fft = 1/N * Complex_conj( fft( Complex_conj( a ) ) )  
//
    for( i=0; i< N; i++){
	    a[i] = std::conj( a[i] );
	}
	fft( a, N );
    for( i=0; i< N; i++){
	    a[i] = std::conj( a[i] ) / static_cast<double>(N) ;
	}
	
}
