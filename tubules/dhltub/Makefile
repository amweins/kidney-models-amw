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

OBJS	      = objects/dhl.o \
 		objects/dhlgam.o \
 		objects/dhlpset.o \
 		objects/dhlptim.o \
 		objects/dhlflux.o \
 		objects/dhlerr.o \
 		objects/dhlresa.o \
 		objects/dhlresb.o \
 		objects/dhlresc.o \
 		objects/dhlresd.o \
 		objects/dhlresf.o \
 		objects/dhlpick.o \
 		objects/dhlreset.o \
		objects/les8.o

PRINT	      = pr

PROGRAM       = dhl

SHELL	      = /bin/sh

SRCS	      = dhl.f \
 		dhlgam.f \
 		dhlpset.f \
 		dhlptim.f \
 		dhlflux.f \
 		dhlerr.f \
 		dhlresa.f \
 		dhlresb.f \
 		dhlresc.f \
 		dhlresd.f \
 		dhlresf.f \
 		dhlpick.f \
 		dhlreset.f \
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


