//
// TimeFrequency.h
//
// C++ interface to TimeFrequency Object
//
//  $Id: TimeFrequency.h,v 1.1 1994/10/04 07:21:04 jak Exp $
//
//  Author: John Kassebaum
//
/* $Log: TimeFrequency.h,v $
/* Revision 1.1  1994/10/04 07:21:04  jak
/* Initial revision
/**/

#ifndef _TimeFrequency_h
#define _TimeFrequency_h

static char rcsid_TimeFrequency_h[] = "$Id: TimeFrequency.h,v 1.1 1994/10/04 07:21:04 jak Exp $";

class Complex;

extern void fft( Complex *, int );
extern void inv_fft( Complex *, int );

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
	Complex         *signal;
	int              signal_length;
	int              isAnalytic;

    TimeFrequency( void );
	~TimeFrequency( void );
	
	inline void setDataWindow( Data_Window awin ) { myWindow = awin; };
	int setWindowSize( unsigned short );
	void setWindowStride( unsigned short );
	inline int getWindowStride( void ) const {  return stride; };
	void setDataSignal( float *, int );
	void setDataSignal( double *, int );
	inline int  getSignalLength( void ) const { return signal_length; };
	inline Complex * getDataSignal( void ) const { return signal; };
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