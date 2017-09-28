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

OBJS	      = ahl.o \
 		ahlgam.o \
 		ahlpset.o \
 		ahlptim.o \
 		ahlflux.o \
 		ahlerr.o \
 		ahlresa.o \
 		ahlresb.o \
 		ahlresc.o \
 		ahlresd.o \
 		ahlresf.o \
 		ahlpick.o \
 		ahlreset.o \
		hksub.o \
		tsc.o \
		nhe3.o \
		nkcc.o \
		kcc.o \
		les8.o


PRINT	      = pr

PROGRAM       = ahl

SHELL	      = /bin/sh

SRCS	      = ahl.f \
 		ahlgam.f \
 		ahlpset.f \
 		ahlptim.f \
 		ahlflux.f \
 		ahlerr.f \
 		ahlresa.f \
 		ahlresb.f \
 		ahlresc.f \
 		ahlresd.f \
 		ahlresf.f \
 		ahlpick.f \
 		ahlreset.f \
		hksub.f \
		tsc.f \
		nhe3.f \
		nkcc.f \
		kcc.f \
		les8.f


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


