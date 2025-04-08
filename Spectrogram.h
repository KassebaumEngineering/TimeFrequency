//
// Spectrogram.h
//
// C++ interface to Short Time Fourier Transform Spectrogram Object
//
//  $Id: Spectrogram.h,v 1.3 1994/10/27 09:11:30 jak Exp $
//
//  Author: John Kassebaum
//
/* $Log: Spectrogram.h,v $
 * Revision 1.3  1994/10/27 09:11:30  jak
 * Fixes, including anti-aliasing additions. -jak
 *
 * Revision 1.2  1994/10/06  17:51:56  jak
 * Made Fixes to several of the spectrum programs - have a preliminary version
 * of the Wigner distribution. -jak
 *
 * Revision 1.1.1.1  1994/10/04  07:21:05  jak
 * Placing Time/Frequency Code under CVS control.  Only Spectrogram
 * works currently.  -jak
 **/

#ifndef _Spectrogram_h
#define _Spectrogram_h

static char rcsid_Spectrogram_h[] = "$Id: Spectrogram.h,v 1.3 1994/10/27 09:11:30 jak Exp $";

#include "TimeFrequency.h"


class Spectrogram : public TimeFrequency {
public:
    int              time_slots;
	int              spectrum_data_length;
	Complex        **spectrogram;
	int              isComputed;

    Spectrogram( void );
	virtual ~Spectrogram( void );
	
	void setWindowStride( unsigned short );
	inline int getTimeSlots( void ) const { return time_slots; };
	inline int getSpectrumDataLength( void ) const { return spectrum_data_length; };
	
	virtual void compute( void );
	
	inline Complex **getSpectrogram( void )
	{ 
	    if (!isComputed) compute();
		return spectrogram;
	}
	
	virtual void print_Mathematica( void );
	virtual void print_Gnuplot( void );
};


#endif