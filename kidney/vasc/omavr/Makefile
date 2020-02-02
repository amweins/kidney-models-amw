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

OBJS	      = omavr.o \
 		omanewt.o \
 		omagam.o \
 		omapset.o \
 		omaflux.o \
 		omaerr.o \
 		omaresa.o \
 		omaresb.o \
 		omareset.o \
		les.o 

PRINT	      = pr

PROGRAM       = omavr

SHELL	      = /bin/sh

SRCS	      = omavr.f \
 		omanewt.f \
 		omagam.f \
 		omapset.f \
 		omaflux.f \
 		omaerr.f \
 		omaresa.f \
 		omaresb.f \
 		omareset.f \
		les.f 

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

