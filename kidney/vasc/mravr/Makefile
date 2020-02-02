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

OBJS	      = mravr.o \
 		mranewt.o \
 		mragam.o \
 		mrapset.o \
 		mraflux.o \
 		mraerr.o \
 		mraresa.o \
 		mraresb.o \
 		mrareset.o \
		les.o 

PRINT	      = pr

PROGRAM       = mravr

SHELL	      = /bin/sh

SRCS	      = mravr.f \
 		mranewt.f \
 		mragam.f \
 		mrapset.f \
 		mraflux.f \
 		mraerr.f \
 		mraresa.f \
 		mraresb.f \
 		mrareset.f \
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

