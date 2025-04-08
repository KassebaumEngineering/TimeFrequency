//
// Choi_Williams.h
//
// C++ interface to Choi_Williams Transform Object
//
//  $Id: Choi_Williams.h,v 1.1 1994/10/27 09:12:39 jak Exp $
//
//  Author: John Kassebaum
//
/* $Log: Choi_Williams.h,v $
/* Revision 1.1  1994/10/27 09:12:39  jak
/* Checking in first working anti-aliased Choi-Williams Distribution TFD -jak
/**/

#ifndef _Choi_Williams_h
#define _Choi_Williams_h

static char rcsid_Choi_Williams_h[] = "$Id: Choi_Williams.h,v 1.1 1994/10/27 09:12:39 jak Exp $";

#include "Wigner.h"


class Choi_Williams : public Wigner {
public:
	double  **CW_Kernel;
    Complex **temp;
    double    sigma;

    Choi_Williams( void );
	virtual ~Choi_Williams( void );
		
	virtual void compute( void );
    inline void setSigma( double a_val ){ sigma = a_val; };
    inline double getSigma( void ){ return sigma; };
		
};


#endif