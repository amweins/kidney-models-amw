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

OBJS	      = imavr.o \
 		imanewt.o \
 		imagam.o \
 		imapset.o \
 		imaflux.o \
 		imaerr.o \
 		imaresa.o \
 		imaresb.o \
 		imareset.o \
		les.o 

PRINT	      = pr

PROGRAM       = imavr

SHELL	      = /bin/sh

SRCS	      = imavr.f \
 		imanewt.f \
 		imagam.f \
 		imapset.f \
 		imaflux.f \
 		imaerr.f \
 		imaresa.f \
 		imaresb.f \
 		imareset.f \
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

