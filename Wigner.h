//
// Wigner.h
//
// C++ interface to Short Time Fourier Transform Wigner Object
//
//  $Id: Wigner.h,v 1.1 1994/10/06 17:52:01 jak Exp $
//
//  Author: John Kassebaum
//
/* $Log: Wigner.h,v $
/* Revision 1.1  1994/10/06 17:52:01  jak
/* Made Fixes to several of the spectrum programs - have a preliminary version
/* of the Wigner distribution. -jak
/**/

#ifndef _Wigner_h
#define _Wigner_h

static char rcsid_Wigner_h[] = "$Id: Wigner.h,v 1.1 1994/10/06 17:52:01 jak Exp $";

#include "Spectrogram.h"

class Complex;

class Wigner : public Spectrogram {
public:

    Wigner( void );
	~Wigner( void );
		
	void compute( void );
		
};


#endif