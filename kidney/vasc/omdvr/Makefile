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

OBJS	      = omdvr.o \
 		omdnewt.o \
 		omdgam.o \
 		omdpset.o \
 		omdflux.o \
 		omderr.o \
 		omdresa.o \
 		omdresb.o \
 		omdreset.o \
		les.o 

PRINT	      = pr

PROGRAM       = omdvr

SHELL	      = /bin/sh

SRCS	      = omdvr.f \
 		omdnewt.f \
 		omdgam.f \
 		omdpset.f \
 		omdflux.f \
 		omderr.f \
 		omdresa.f \
 		omdresb.f \
 		omdreset.f \
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

