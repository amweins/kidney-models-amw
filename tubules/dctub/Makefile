F95 = gfortran

FFLAGS	      = -g -w

CFLAGS	      = 

DEST	      = .

EXTHDRS	      =

HDRS	      =

INSTALL	      = install

LD	      = gfortran

LDFLAGS	      = 

LIBS	      =

MAKEFILE      = Makefile

OBJS	      = dct.o \
 		dctgam.o \
 		dctpset.o \
 		dctptim.o \
 		dctflux.o \
 		dcterr.o \
 		dctresa.o \
 		dctresb.o \
 		dctresc.o \
 		dctresd.o \
 		dctpick.o \
 		dcreset.o \
		hksub.o \
		tsc.o \
		nhe3.o \
		les8.o

PRINT	      = pr

PROGRAM       = dct

SHELL	      = /bin/sh

SRCS	      = dct.f \
 		dctgam.f \
 		dctpset.f \
 		dctptim.f \
 		dctflux.f \
 		dcterr.f \
 		dctresa.f \
 		dctresb.f \
 		dctresc.f \
 		dctresd.f \
 		dctpick.f \
 		dcreset.f \
		hksub.f \
		tsc.f \
		nhe3.f \
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

