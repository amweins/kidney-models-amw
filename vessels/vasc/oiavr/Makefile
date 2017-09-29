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

OBJS	      = oiavr.o \
 		oianewt.o \
 		oiagam.o \
 		oiapset.o \
 		oiaflux.o \
 		oiaerr.o \
 		oiaresa.o \
 		oiaresb.o \
 		oiareset.o \
		les.o 

PRINT	      = pr

PROGRAM       = oiavr

SHELL	      = /bin/sh

SRCS	      = oiavr.f \
 		oianewt.f \
 		oiagam.f \
 		oiapset.f \
 		oiaflux.f \
 		oiaerr.f \
 		oiaresa.f \
 		oiaresb.f \
 		oiareset.f \
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

