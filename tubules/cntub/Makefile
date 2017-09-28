F95 = gfortran

FFLAGS	      = -g

CFLAGS	      = 

DEST	      = .

EXTHDRS	      =

HDRS	      =

INSTALL	      = install

LD	      = gfortran

LDFLAGS	      = 

LIBS	      =

MAKEFILE      = Makefile

OBJS	      = cnt.o \
 		cntgam.o \
 		cntpset.o \
 		cntptim.o \
 		cntflux.o \
 		cnterr.o \
 		cntresa.o \
 		cntresb.o \
 		cntresc.o \
 		cntresd.o \
 		cntpick.o \
 		cnreset.o \
		hksub.o \
		ae1.o \
		les8.o

PRINT	      = pr

PROGRAM       = cnt

SHELL	      = /bin/sh

SRCS	      = cnt.f \
 		cntgam.f \
 		cntpset.f \
 		cntptim.f \
 		cntflux.f \
 		cnterr.f \
 		cntresa.f \
 		cntresb.f \
 		cntresc.f \
 		cntresd.f \
 		cntpick.f \
 		cnreset.f \
		hksub.f \
		ae1.f \
		les8.f

$(PROGRAM):     $(OBJS) $(LIBS)
		echo "Linking $(PROGRAM) ..."
		$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)
		echo "done"

$(OBJS):	$(SRCS)
		gfortran -c $(FFLAGS) $(SRCS)

install:	$(PROGRAM)
		echo Installing $(PROGRAM) in $(DEST)
		if [ $(DEST) != . ]; then \
		(rm -f $(DEST)/$(PROGRAM); $(INSTALL) -f $(DEST) $(PROGRAM)); fi

clean:;		rm -f $(OBJS) core

