#
# Makefile for Spectral Analysis Tools
# 
# $Id: Makefile,v 1.2 1994/10/06 17:51:53 jak Exp $
#
# History:
# $Log: Makefile,v $
# Revision 1.2  1994/10/06 17:51:53  jak
# Made Fixes to several of the spectrum programs - have a preliminary version
# of the Wigner distribution. -jak
#
# Revision 1.1.1.1  1994/10/04  07:21:05  jak
# Placing Time/Frequency Code under CVS control.  Only Spectrogram
# works currently.  -jak
#
#
#

CC=gcc
CFLAGS=-g -O2 -V2.5.8
LIBS=-lm -lg++

SOURCES= Makefile main.cc fft.cc TimeFrequency.h TimeFrequency.cc Spectrogram.cc Spectrogram.h sine.c Wigner.cc Wigner.h
OBJECTS= fft.o Spectrogram.o TimeFrequency.o Wigner.o
EXES= main
LINKS= Spectrogram Wigner Choi-Williams

DEPENDOBJ = $(OBJECTS:.o=.d) 

all: $(EXES) $(LINKS)

$(LINKS): %: main 
	-ln -s $< $@ 

%.o: %.cc
	$(CC) -c $(CFLAGS) $< -o $@

sine: sine.o
	$(CC) $< -o $@

$(EXES): %: %.o $(OBJECTS)
	$(CC) $(CFLAGS) $< $(OBJECTS) $(LIBS) -o $@ 

clean:
	rm -f *.o 
	
distclean:
	rm -f *.o *.d $(EXES) $(LINKS) sine

install:
	cp main $(HOME)/Apps/tfd
	ln -s ./tfd $(HOME)/Apps/Spectrogram
	ln -s ./tfd $(HOME)/Apps/Wigner
	ln -s ./tfd $(HOME)/Apps/Choi-Williams

$(DEPENDOBJ): %.d: %.cc
	$(SHELL) -ec '$(CC) -M -c $(CFLAGS) $< | sed '\''s/$*.o/& $@/g'\'' > $@'

# Add Dependancies
include $(DEPENDOBJ)
