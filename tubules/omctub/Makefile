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

OBJS	      = omct.o \
 		omctgam.o \
 		omctpset.o \
 		omctptim.o \
 		omctflux.o \
 		omcterr.o \
 		omctresa.o \
 		omctresb.o \
 		omctresc.o \
 		omctresd.o \
 		omctpick.o \
 		omcreset.o \
		hksub.o \
		ae1.o \
		les8.o

PRINT	      = pr

PROGRAM       = omct

SHELL	      = /bin/sh

SRCS	      = omct.f \
 		omctgam.f \
 		omctpset.f \
 		omctptim.f \
 		omctflux.f \
 		omcterr.f \
 		omctresa.f \
 		omctresb.f \
 		omctresc.f \
 		omctresd.f \
 		omctpick.f \
 		omcreset.f \
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

