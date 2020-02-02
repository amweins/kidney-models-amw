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

OBJS	      = oidvr.o \
 		oidnewt.o \
 		oidgam.o \
 		oidpset.o \
 		oidflux.o \
 		oiderr.o \
 		oidresa.o \
 		oidresb.o \
 		oidreset.o \
		les.o 

PRINT	      = pr

PROGRAM       = oidvr

SHELL	      = /bin/sh

SRCS	      = oidvr.f \
 		oidnewt.f \
 		oidgam.f \
 		oidpset.f \
 		oidflux.f \
 		oiderr.f \
 		oidresa.f \
 		oidresb.f \
 		oidreset.f \
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

