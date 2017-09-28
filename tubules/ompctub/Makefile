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

OBJS	      = ompct.o \
 		ompctgam.o \
 		ompctpset.o \
 		ompctptim.o \
 		ompctflux.o \
 		ompcterr.o \
 		ompctresa.o \
 		ompctresb.o \
 		ompctresc.o \
 		ompctresd.o \
 		ompctpick.o \
 		ompctreset.o \
		hksub.o \
		ae1.o \
		les8.o

PRINT	      = pr

PROGRAM       = ompct

SHELL	      = /bin/sh

SRCS	      = ompct.f \
 		ompctgam.f \
 		ompctpset.f \
 		ompctptim.f \
 		ompctflux.f \
 		ompcterr.f \
 		ompctresa.f \
 		ompctresb.f \
 		ompctresc.f \
 		ompctresd.f \
 		ompctpick.f \
 		ompctreset.f \
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

