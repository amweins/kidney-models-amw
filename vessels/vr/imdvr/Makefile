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

OBJS	      = imdvr.o \
 		imdnewt.o \
 		imdgam.o \
 		imdpset.o \
 		imdflux.o \
 		imderr.o \
 		imdresa.o \
 		imdresb.o \
 		imdreset.o \
		les.o 

PRINT	      = pr

PROGRAM       = imdvr

SHELL	      = /bin/sh

SRCS	      = imdvr.f \
 		imdnewt.f \
 		imdgam.f \
 		imdpset.f \
 		imdflux.f \
 		imderr.f \
 		imdresa.f \
 		imdresb.f \
 		imdreset.f \
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

