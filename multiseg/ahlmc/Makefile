F95 = gfortran

FFLAGS	      = -g -w

CFLAGS	      = 

DEST	      = objects

EXTHDRS	      =

HDRS	      =

INSTALL	      = install

LD	      = gfortran

LDFLAGS	      = 

LIBS	      =

MAKEFILE      = Makefile

OBJS	      = objects/ahlmc.o \
 		objects/ahlmnewt.o \
 		objects/ahlmgam.o \
 		objects/ahlmpset.o \
 		objects/ahlmptim.o \
 		objects/ahlmflux.o \
 		objects/ahlmerr.o \
 		objects/ahlmresa.o \
 		objects/ahlmresb.o \
 		objects/ahlmresc.o \
 		objects/ahlmresd.o \
 		objects/ahlmresf.o \
 		objects/ahlmpick.o \
 		objects/ahlmreset.o \
 		objects/ahlcnewt.o \
 		objects/ahlcgam.o \
 		objects/ahlcpset.o \
 		objects/ahlcptim.o \
 		objects/ahlcflux.o \
 		objects/ahlcerr.o \
 		objects/ahlcresa.o \
 		objects/ahlcresb.o \
 		objects/ahlcresc.o \
 		objects/ahlcresd.o \
 		objects/ahlcresf.o \
 		objects/ahlcpick.o \
 		objects/ahlcreset.o \
		objects/hksub.o \
		objects/ae1.o \
		objects/nhe3.o \
		objects/tsc.o \
		objects/nkcc.o \
		objects/kcc.o \
		objects/les8.o \
		objects/eqlib.o

PRINT	      = pr

PROGRAM       = ahlmc

SHELL	      = /bin/sh

SRCS	      = ahlmc.f \
 		ahlm/ahlmnewt.f \
 		ahlm/ahlmgam.f \
 		ahlm/ahlmpset.f \
 		ahlm/ahlmptim.f \
 		ahlm/ahlmflux.f \
 		ahlm/ahlmerr.f \
 		ahlm/ahlmresa.f \
 		ahlm/ahlmresb.f \
 		ahlm/ahlmresc.f \
 		ahlm/ahlmresd.f \
 		ahlm/ahlmresf.f \
 		ahlm/ahlmpick.f \
 		ahlm/ahlmreset.f \
 		ahlc/ahlcnewt.f \
 		ahlc/ahlcgam.f \
 		ahlc/ahlcpset.f \
 		ahlc/ahlcptim.f \
 		ahlc/ahlcflux.f \
 		ahlc/ahlcerr.f \
 		ahlc/ahlcresa.f \
 		ahlc/ahlcresb.f \
 		ahlc/ahlcresc.f \
 		ahlc/ahlcresd.f \
 		ahlc/ahlcresf.f \
 		ahlc/ahlcpick.f \
 		ahlc/ahlcreset.f \
		hksub.f \
		ae1.f \
		nhe3.f \
		tsc.f \
		nkcc.f \
		kcc.f \
		les8.f \
		eqlib.f


all:		$(PROGRAM)

$(PROGRAM):     $(OBJS) $(LIBS)
		@echo "Linking $(PROGRAM) ..."
		@$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)
		@echo "done"

$(OBJS):	$(SRCS)
		gfortran -c $(FFLAGS) $(SRCS); mv *.o ./objects   

clean:;		@rm -f $(OBJS) core

clobber:;	@rm -f $(OBJS) $(PROGRAM) core tags

install:	$(PROGRAM)
		echo Installing $(PROGRAM) in $(DEST)
		if [ $(DEST) != . ]; then \
		(rm -f $(DEST)/$(PROGRAM); $(INSTALL) -f $(DEST) $(PROGRAM)); fi


