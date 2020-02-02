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

OBJS	      = mrdvr.o \
 		mrdnewt.o \
 		mrdgam.o \
 		mrdpset.o \
 		mrdflux.o \
 		mrderr.o \
 		mrdresa.o \
 		mrdresb.o \
 		mrdreset.o \
		les.o 

PRINT	      = pr

PROGRAM       = mrdvr

SHELL	      = /bin/sh

SRCS	      = mrdvr.f \
 		mrdnewt.f \
 		mrdgam.f \
 		mrdpset.f \
 		mrdflux.f \
 		mrderr.f \
 		mrdresa.f \
 		mrdresb.f \
 		mrdreset.f \
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

