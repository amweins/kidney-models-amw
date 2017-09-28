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

OBJS	      = cct.o \
 		cctgam.o \
 		cctpset.o \
 		cctptim.o \
 		cctflux.o \
 		ccterr.o \
 		cctresa.o \
 		cctresb.o \
 		cctresc.o \
 		cctresd.o \
 		cctpick.o \
 		ccreset.o \
		hksub.o \
		ae1.o \
		les8.o

PRINT	      = pr

PROGRAM       = cct

SHELL	      = /bin/sh

SRCS	      = cct.f \
 		cctgam.f \
 		cctpset.f \
 		cctptim.f \
 		cctflux.f \
 		ccterr.f \
 		cctresa.f \
 		cctresb.f \
 		cctresc.f \
 		cctresd.f \
 		cctpick.f \
 		ccreset.f \
		hksub.f \
		ae1.f \
		les8.f

all:		$(PROGRAM)

$(PROGRAM):     $(OBJS) $(LIBS)
		@echo "Linking $(PROGRAM) ..."
		@$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)
		@echo "done"

$(OBJS):	$(SRCS)
		gfortran -c $(FFLAGS) $(SRCS)

clean:;		@rm -f $(OBJS) core

clobber:;	@rm -f $(OBJS) $(PROGRAM) core tags

install:	$(PROGRAM)
		echo Installing $(PROGRAM) in $(DEST)
		if [ $(DEST) != . ]; then \
		(rm -f $(DEST)/$(PROGRAM); $(INSTALL) -f $(DEST) $(PROGRAM)); fi


