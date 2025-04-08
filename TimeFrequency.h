//
// TimeFrequency.h
//
// C++ interface to TimeFrequency Object
//
//  $Id: TimeFrequency.h,v 1.2 1994/10/06 17:51:58 jak Exp $
//
//  Author: John Kassebaum
//
/* $Log: TimeFrequency.h,v $
/* Revision 1.2  1994/10/06 17:51:58  jak
/* Made Fixes to several of the spectrum programs - have a preliminary version
/* of the Wigner distribution. -jak
/*
 * Revision 1.1.1.1  1994/10/04  07:21:05  jak
 * Placing Time/Frequency Code under CVS control.  Only Spectrogram
 * works currently.  -jak
 **/

#ifndef _TimeFrequency_h
#define _TimeFrequency_h

static char rcsid_TimeFrequency_h[] = "$Id: TimeFrequency.h,v 1.2 1994/10/06 17:51:58 jak Exp $";

#include <complex>
using Complex = std::complex<double>;

extern void fft( std::complex<double> *, int );
extern void inv_fft( std::complex<double> *, int );

enum Data_Window {
    Hanning = 0,
	Hamming = 1,
	Rectangular = 2
};

class TimeFrequency {
public:
    Data_Window      myWindow;
    unsigned short   window_size;
	unsigned short   stride;
	double           sampling_frequency;
	double           frequency_resolution;
	double           frequency_band;
	std::complex<double> *signal;
	int              signal_length;
	int              isAnalytic;

    TimeFrequency( void );
	virtual ~TimeFrequency( void );
	
	inline void setDataWindow( Data_Window awin ) { myWindow = awin; };
	int setWindowSize( unsigned short );
	void setWindowStride( unsigned short );
	inline int getWindowStride( void ) const {  return stride; };
	void setDataSignal( float *, int );
	void setDataSignal( double *, int );
	inline int  getSignalLength( void ) const { return signal_length; };
	inline std::complex<double> * getDataSignal( void ) const { return signal; };
	inline Data_Window getDataWindow( void ) const { return myWindow; };
	inline int getWindowSize( void ) const { return window_size; };
	
	void makeAnalytic( void ); // Will make the current signal analytic
	                           // via a Hilbert Transform Filter
	
	void setSamplingFrequency( double );
	inline double getSamplingFrequency( void ) const 
	{ 
	    return sampling_frequency; 
	};
	
	int     setFrequencyResolution( double );
	double  getFrequencyResolution( void );
    double  getBandWidth( void );
};


#endif