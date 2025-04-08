#
# Makefile for Spectral Analysis Tools
# 
# $Id: Makefile,v 1.5 1994/12/06 05:28:15 jak Exp $
#
# History:
# $Log: Makefile,v $
# Revision 1.5  1994/12/06 05:28:15  jak
# Removed error condition in makefile install target. -jak
#
# Revision 1.4  1994/10/27  09:11:28  jak
# Fixes, including anti-aliasing additions. -jak
#
# Revision 1.3  1994/10/07  06:55:26  jak
# Wigner now works!  Bug fixes to the Spectrogram also.  Stride can now
# be set from the command line!  -jak
#
# Revision 1.2  1994/10/06  17:51:53  jak
# Made Fixes to several of the spectrum programs - have a preliminary version
# of the Wigner distribution. -jak
#
# Revision 1.1.1.1  1994/10/04  07:21:05  jak
# Placing Time/Frequency Code under CVS control.  Only Spectrogram
# works currently.  -jak
#
#
#

CXX=clang++
CC=clang
CXXFLAGS=-std=c++20 -Wall -Wextra -Wpedantic -ggdb -O2
CFLAGS=-Wall -Wextra -pedantic -ggdb -O2
LDLIBS=-lm

SOURCES= Makefile main.cc fft.cc TimeFrequency.h TimeFrequency.cc Spectrogram.cc Spectrogram.h sine.c Wigner.cc Wigner.h Choi_Williams.cc Choi_Williams.h
OBJECTS= fft.o Spectrogram.o TimeFrequency.o Wigner.o Choi_Williams.o
EXES= main
LINKS= Spectrogram Wigner Choi_Williams

DEPENDOBJ = $(OBJECTS:.o=.d) 

all: $(EXES) $(LINKS)

$(LINKS): %: main 
	-ln -s $< $@ 

%.o: %.cc
	$(CXX) -c $(CXXFLAGS) $< -o $@

sine: sine.o
	$(CC) $(CFLAGS) $< -o $@

$(EXES): %: %.o $(OBJECTS)
	$(CXX) $(CXXFLAGS) $< $(OBJECTS) $(LDLIBS) -o $@

clean:
	rm -f *.o 
	
distclean:
	rm -f *.o *.d $(EXES) $(LINKS) sine

install:
	cp main $(HOME)/Apps/tfd
	-ln -s ./tfd $(HOME)/Apps/Spectrogram
	-ln -s ./tfd $(HOME)/Apps/Wigner
	-ln -s ./tfd $(HOME)/Apps/Choi_Williams

$(DEPENDOBJ): %.d: %.cc
	$(SHELL) -ec '$(CC) -M -c $(CFLAGS) $< | sed '\''s/$*.o/& $@/g'\'' > $@'

# Add Dependancies
include $(DEPENDOBJ)
