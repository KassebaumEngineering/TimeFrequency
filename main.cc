//
// main.cc
//
// C++ Main Program for Calculating Data Distributions 
//
//  $Id: main.cc,v 1.5 1994/11/18 05:52:48 jak Exp $
//
//  Author: John Kassebaum
//
// $Log: main.cc,v $
// Revision 1.5  1994/11/18 05:52:48  jak
// Small changes to improve Choi_Williams operations. -jak
//
// Revision 1.4  1994/10/27  09:11:34  jak
// Fixes, including anti-aliasing additions. -jak
//
// Revision 1.3  1994/10/07  06:55:33  jak
// Wigner now works!  Bug fixes to the Spectrogram also.  Stride can now
// be set from the command line!  -jak
//
// Revision 1.2  1994/10/06  17:52:02  jak
// Made Fixes to several of the spectrum programs - have a preliminary version
// of the Wigner distribution. -jak
//
// Revision 1.1.1.1  1994/10/04  07:21:05  jak
// Placing Time/Frequency Code under CVS control.  Only Spectrogram
// works currently.  -jak
//
//
//

static char rcsid_main_cc[] = "$Id: main.cc,v 1.5 1994/11/18 05:52:48 jak Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Spectrogram.h"
#include "Wigner.h"
#include "Choi_Williams.h"

//
// Usage Info -------------------------------------------------------------
//
#define USAGE    "Usage: %s [args] < input data_file > output data_file\n"
#define DESCR_0    "Where args are one or more of:\n"
#define DESCR_1    " -h            prints this help message and exits.\n"
#define DESCR_2    " -v            verbosely prints out data\n"
#define DESCR_3    " -s            prints statistics\n"
#define DESCR_4    " -w  <int>     set sampling window width (in samples)\n"
#define DESCR_4_1  "   (default %d samples)\n"
#define DESCR_5    " -ss <int>     set spectrum spacing (in samples)\n"
#define DESCR_6    " -r  <float>   set sampling frequency value (in Hz)\n"
#define DESCR_6_1  "   (default %d Hz)\n"
#define DESCR_7    " -bi           set input file type to binary (float)\n"
#define DESCR_8    " -bo           set output file type to binary (float)\n"
#define DESCR_9    " -no           request no output\n"
#define DESCR_10   " -d <string>   set distribution type (default FFT)\n"
#define DESCR_10_1 "    choices: FFT, Wigner, Choi_Williams, None\n"
#define DESCR_11   " -sigma <float>    set sigma value (useful for Choi_WIlliams only)\n"
#define DESCR_12   "\n data_file  -  A file of floating point numbers\n"
#define DESCR_13   "(default) in ascii format, with one number per line. \n "
#define DESCR_14   "example:\n"
#define DESCR_15   " -47.52\n -39.59\n -27.72\n -15.84\n"
//
// ------------------------------------------------------------------------
//

#define DEFAULT_WINDOW			128 	// samples
#define DEFAULT_SAMPLING_RATE	1.0 	// Hz

enum TFD { 
   STFT = 0,  // Short Time Fourier Transform
   WT   = 1,  // Wigner Transform
   CWD  = 2,  // Choi- Williams Distribution
   NONE = 3,  // None - request no distribution calculation (useful for binary to float conversions)
   NA         // None Applicable - requested unknown type
}; 

#define BLOCKSIZE    1024
#define Null(A)             ((A *) 0)
#define New(A)              ((A *) malloc( sizeof(A) ) )
#define NewBlock(A,N)       ((A *) malloc( sizeof(A) * (N)) )
#define BiggerBlock(A,B,N)  ((A *) realloc( (void *)(B), sizeof(A) * (N)))

main(int argc,char **argv)
{
    int c,i,j,k;
    int verbose, stats, flag;
	int window_width, stride, binary_input, binary_output, no_output;
	float sampling_rate, sigma;
	TFD  chosen_distribution;
    long int data_count;
    long int space_alloc;
    float min_value, max_value;
    double tempf;
    int step_count;
    float step;
    float *data;
    float *output_data;

    verbose = 0;
    stats = 0;
	window_width = DEFAULT_WINDOW;
    stride = DEFAULT_WINDOW/2;
	binary_input = 0;
	binary_output = 0;
	no_output = 0;
	sampling_rate = DEFAULT_SAMPLING_RATE;
	chosen_distribution = STFT;
    sigma = window_width*window_width;
	
  //
  // Special Check for program name - may assign default distribution type
  //
	if ( !strcmp( argv[0],"Spectrogram") ){
	    chosen_distribution = STFT;
	} else if ( !strcmp( argv[0],"Wigner") ){
	    chosen_distribution = WT;
	} else if ( !strcmp( argv[0],"Choi_Williams") ){
	    chosen_distribution = CWD;
	}
	
  //
  // Parse Command Line
  //
    for (c=1; c< argc; c++) {
        if ( !strcmp( argv[ c ],"-h") || !strcmp( argv[ c ],"-help")){
            fprintf(stderr,USAGE,argv[0]);
            fprintf(stderr,DESCR_0);
            fprintf(stderr,DESCR_1);
            fprintf(stderr,DESCR_2);
            fprintf(stderr,DESCR_3);
            fprintf(stderr,DESCR_4);
            fprintf(stderr,DESCR_4_1, DEFAULT_WINDOW);
            fprintf(stderr,DESCR_5);
            fprintf(stderr,DESCR_6);
            fprintf(stderr,DESCR_6_1, DEFAULT_SAMPLING_RATE);
            fprintf(stderr,DESCR_7);
            fprintf(stderr,DESCR_8);
            fprintf(stderr,DESCR_9);
            fprintf(stderr,DESCR_10);
            fprintf(stderr,DESCR_10_1);
            fprintf(stderr,DESCR_11);
            fprintf(stderr,DESCR_12);
            fprintf(stderr,DESCR_13);
            fprintf(stderr,DESCR_14);
            fprintf(stderr,DESCR_15);
            exit(0);
        } else if (!strcmp( argv[ c ],"-v")){
            verbose = 1;
        } else if (!strcmp( argv[ c ],"-s")){
            stats = 1;
        } else if (!strcmp( argv[ c ],"-w")){
		    c++;
			window_width = atoi( argv[ c ] );
        } else if (!strcmp( argv[ c ],"-ss")){
		    c++;
            stride = atoi( argv[ c ] );
        } else if (!strcmp( argv[ c ],"-r")){
		    c++;
			sampling_rate = atof( argv[ c ] );
        } else if (!strcmp( argv[ c ],"-sigma")){
		    c++;
			sigma = atof( argv[ c ] );
        } else if (!strcmp( argv[ c ],"-bi")){
			binary_input = 1;
        } else if (!strcmp( argv[ c ],"-bo")){
			binary_output = 1;
        } else if (!strcmp( argv[ c ],"-no")){
			no_output = 1;
        } else if (!strcmp( argv[ c ],"-d")){
		    c++;
			if( !strcmp( argv[ c ],"FFT") ){
			    chosen_distribution = STFT;
			} else if( !strcmp( argv[ c ],"Wigner") ){
			    chosen_distribution = WT;
			} else if( !strcmp( argv[ c ],"Choi_Williams") ){
			    chosen_distribution = CWD;
			} else if( !strcmp( argv[ c ],"None") ){
			    chosen_distribution = NONE;
			} else {
			    chosen_distribution = NA;
				fprintf(stderr,"%s: Unknown Distribution Type '%s'!\n",argv[0],argv[ c ]);
				abort();
			}
        }
    }
    
    data_count = 0;
    space_alloc = 0;
    data = (float *)0;
    
    if( data = NewBlock( float, BLOCKSIZE ) ){
        space_alloc = BLOCKSIZE;
    } else {
        perror(argv[0]);
        abort();
    };
        
	if (binary_input) {
	  //
	  // binary input data
	  //
	    int fd, bytes_read, block_size;
		fd = 0;    // note: stdin is always fd = 0
		
		block_size = BLOCKSIZE;
	    while ((bytes_read = read( fd, (char *) &data[ data_count ], (block_size*sizeof(float)) ))
		                   > 0 )
//	    while ((bytes_read = fread( (char *) &data[ data_count ], sizeof(float), block_size, stdin ))
//		                   > 0 )
		{
		    data_count += ( bytes_read / sizeof(float) );
			if( (block_size = (space_alloc - data_count)) < 1) {
				if( data = BiggerBlock( float, data, (space_alloc + BLOCKSIZE) ) ){
					space_alloc += BLOCKSIZE;
					block_size = space_alloc - data_count;
				} else {
					perror(argv[0]);
					abort();
				}
			}
		}
		
	    if( bytes_read == -1 ) {
			perror(argv[0]);
			abort();
		}
		
		if( verbose ){
		    for( i=0; i< data_count; i++){
		        fprintf(stderr, "%f\n", data[ i ]);
			}
		}
		
		if( stats ){
		    int i;
		    flag = 0;
		    for( i=0; i< data_count; i++){
				if (!flag) {
					min_value = max_value = data[ i ];
					flag = 1;
				} else {
					if( data[ i ] < min_value )
						min_value = data[ i ];
					if( data[ i ] > max_value )
						max_value = data[ i ];
				}
			}
		} 
		
 	} else { 
	  //
	  // ascii input data
	  //
		flag = 0;
		while ( scanf( "%f", &data[ data_count ]) != EOF ){
			if (!flag) {
				min_value = max_value = data[ data_count ];
				flag = 1;
			} else {
				if( data[ data_count ] < min_value )
					min_value = data[ data_count ];
				if( data[ data_count ] > max_value )
					max_value = data[ data_count ];
			}
			if (verbose) fprintf(stdout, "%f\n", data[ data_count ]);
			data_count++;
			if ( data_count >= space_alloc ){
				if( data = BiggerBlock( float, data, (space_alloc + BLOCKSIZE) ) ){
					space_alloc += BLOCKSIZE;
				} else {
					perror(argv[0]);
					abort();
				};
			};
		}
	}
    
    if (stats) {
        fprintf(stderr, "Space Allocated for %d floats.\n", space_alloc);
        fprintf(stderr, "Space Used for %d floats.\n", data_count);
        fprintf(stderr, "Minimum value = %f.\n", min_value);
        fprintf(stderr, "Maximum value = %f.\n", max_value);
    }
        
	//
	// Call subprogram to calculate requested time/frequency distribution --------------
	//
        switch( chosen_distribution ) {
		    case STFT: {
			    Spectrogram mySpectrogram;
				
			    //fprintf(stderr,"Call to STFT\n");
				no_output = 1;
				
				mySpectrogram.setDataSignal( data, data_count );
				mySpectrogram.setSamplingFrequency( sampling_rate );
				mySpectrogram.setWindowSize( window_width );
				mySpectrogram.setWindowStride( stride );
				if (stats) {
					fprintf(stderr, "Sampling Frequency %f => %f Hz.\n", sampling_rate,  mySpectrogram.getSamplingFrequency());
					fprintf(stderr, "Window Width %d => %d samples.\n", window_width, mySpectrogram.getWindowSize());
					fprintf(stderr, "Window Stride %d => %d samples.\n", stride, mySpectrogram.getWindowStride());
					fprintf(stderr, "Time Slots = %d \n", mySpectrogram.getTimeSlots());
					fprintf(stderr, "Spectrum Data length = %d \n", mySpectrogram.getSpectrumDataLength());
 				}
				mySpectrogram.compute();
                //mySpectrogram.print_Mathematica();
                mySpectrogram.print_Gnuplot();
				
				break;
			}
		    case WT: {
			    Wigner myWigner;
			
			    fprintf(stderr,"Call to WT\n");
				no_output = 1;

				myWigner.setDataSignal( data, data_count );
				myWigner.setSamplingFrequency( sampling_rate );
				myWigner.setWindowSize( window_width );
				myWigner.setWindowStride( stride );
				if (stats) {
					fprintf(stderr, "Sampling Frequency %f => %f Hz.\n", sampling_rate,  myWigner.getSamplingFrequency());
					fprintf(stderr, "Window Width %d => %d samples.\n", window_width, myWigner.getWindowSize());
					fprintf(stderr, "Window Stride %d => %d samples.\n", stride, myWigner.getWindowStride());
					fprintf(stderr, "Time Slots = %d \n", myWigner.getTimeSlots());
					fprintf(stderr, "Spectrum Data length = %d \n", myWigner.getSpectrumDataLength());
 				}
				myWigner.compute();
                //myWigner.print_Mathematica();
                myWigner.print_Gnuplot();
				
				break;
			}
		    case CWD:{
			    Choi_Williams myCW;

			    fprintf(stderr,"Call to CWD\n");
				no_output = 1;

				myCW.setDataSignal( data, data_count );
				myCW.setSamplingFrequency( sampling_rate );
				myCW.setWindowSize( window_width );
				myCW.setWindowStride( stride );
				myCW.setSigma( sigma );
				if (stats) {
					fprintf(stderr, "Sampling Frequency %f => %f Hz.\n", sampling_rate,  myCW.getSamplingFrequency());
					fprintf(stderr, "Window Width %d => %d samples.\n", window_width, myCW.getWindowSize());
					fprintf(stderr, "Window Stride %d => %d samples.\n", stride, myCW.getWindowStride());
					fprintf(stderr, "Sigma => %f \n", myCW.getSigma());
					fprintf(stderr, "Time Slots = %d \n", myCW.getTimeSlots());
					fprintf(stderr, "Spectrum Data length = %d \n", myCW.getSpectrumDataLength());
 				}
				myCW.compute();
                //myCW.print_Mathematica();
                myCW.print_Gnuplot();

				break;
            }
		    case NONE:
			    fprintf(stderr,"No Distribution Calculated.\n");
		        output_data = data;
		        break;
		}
	//
	// ---------------------------------------------------------------------------------
	//
	
	if( !no_output ) {
		if( binary_output ) {
			//
			// Binary Output Data
			//
			int fd;
			fd = 1; // stdout is always file descriptor 1
			if( write( fd, output_data, data_count*sizeof(float) ) == -1){
				perror(argv[0]);
				abort();
			}
		} else {
			//
			// Ascii Output Data
			//
			for(i=0; i< data_count; i++){
				fprintf(stdout,"%f\n",output_data[i]);
			}
		}
	}
	
    free( data );
}
