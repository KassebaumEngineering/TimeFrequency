//
// Spectrogram.h
//
// C++ interface to Short Time Fourier Transform Spectrogram Object
//
//  $Id: Spectrogram.h,v 1.1 1994/10/04 07:21:04 jak Exp $
//
//  Author: John Kassebaum
//
/* $Log: Spectrogram.h,v $
/* Revision 1.1  1994/10/04 07:21:04  jak
/* Initial revision
/**/

#ifndef _Spectrogram_h
#define _Spectrogram_h

static char rcsid_Spectrogram_h[] = "$Id: Spectrogram.h,v 1.1 1994/10/04 07:21:04 jak Exp $";

#include "TimeFrequency.h"

class Complex;

class Spectrogram : public TimeFrequency {
public:
    int              time_slots;
	int              spectrum_data_length;
	Complex        **spectrogram;
	int              isComputed;

    Spectrogram( void );
	~Spectrogram( void );
	
	void setWindowStride( unsigned short );
	inline int getTimeSlots( void ) const { return time_slots; };
	inline int getSpectrumDataLength( void ) const { return spectrum_data_length; };
	
	void compute( void );
	
	inline Complex **getSpectrogram( void )
	{ 
	    if (!isComputed) compute();
		return spectrogram;
	}
	
	void print_Mathematica( void );
	void print_Gnuplot( void );
};


#endif