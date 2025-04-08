//
// Wigner.h
//
// C++ interface to Short Time Fourier Transform Wigner Object
//
//  $Id: Wigner.h,v 1.3 1994/10/27 09:11:33 jak Exp $
//
//  Author: John Kassebaum
//
/* $Log: Wigner.h,v $
/* Revision 1.3  1994/10/27 09:11:33  jak
/* Fixes, including anti-aliasing additions. -jak
/*
 * Revision 1.2  1994/10/07  06:55:32  jak
 * Wigner now works!  Bug fixes to the Spectrogram also.  Stride can now
 * be set from the command line!  -jak
 *
 * Revision 1.1  1994/10/06  17:52:01  jak
 * Made Fixes to several of the spectrum programs - have a preliminary version
 * of the Wigner distribution. -jak
 **/

#ifndef _Wigner_h
#define _Wigner_h

static char rcsid_Wigner_h[] = "$Id: Wigner.h,v 1.3 1994/10/27 09:11:33 jak Exp $";

#include "Spectrogram.h"


class Wigner : public Spectrogram {
public:

    Wigner( void );
	virtual ~Wigner( void );
		
	virtual void compute( void );
    virtual void print_Gnuplot( void );
};


#endif