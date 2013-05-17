# This material contains unpublished, proprietary software of 
# Entropic Research Laboratory, Inc. Any reproduction, distribution, 
# or publication of this work must be authorized in writing by Entropic 
# Research Laboratory, Inc., and must bear the notice: 
#
#    "Copyright (c) 1990-1993  Entropic Research Laboratory, Inc. 
#                   All rights reserved"
#
# The copyright notice above does not evidence any actual or intended 
# publication of this source code.     
#
# @(#)Makefile	1.3 1/22/97 ERL
# 
# Makefile for: 
#
# Written by:  Derek Lin
# Checked by:
# Revised by:

# dpwe
#
DPWEDIR=../dpwelib
DPWEINC=${DPWEDIR}/dpwelib

PROGCFLAGS = -I. -I${DPWEDIR}/dpwelib
LIBFLAGS = -lm -L${DPWEDIR} -ldpwe

CFLAGS = -g $(PROGCFLAGS)

OBJS =  get_f0_snd.o get_cands.o sigproc.o dp_f0.o dpwe_add.o
SRCS =  get_f0_snd.c get_cands.c sigproc.c dp_f0.c dpwe_add.c
PROGNAME = get_f0_snd
MANNAME = get_f0.1
DATA = mpgr1_sx419.wav mpgr1_sx419.txt-ref get_f0_params
DOCS = README.md LICENSE-ESPS.txt get_f0.1


$(PROGNAME): $(OBJS)
	cc  $(CFLAGS) $(OBJS) $(LIBFLAGS) -o $(PROGNAME)

clean:	
	-rm -f $(OBJS) $(PROGNAME) core 

test:	$(PROGNAME)
	./$(PROGNAME) mpgr1_sx419.wav tmp.out
	diff mpgr1_sx419.txt-ref tmp.out
	rm tmp.out


