F90 = gfortran

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

OBJS	      = pct.o \
		pcnewton.o \
 		pctgam.o \
 		pctpset.o \
 		pctptim.o \
 		pctflux.o \
 		pcterr.o \
 		pctresa.o \
 		pctresb.o \
 		pctresc.o \
 		pctresd.o \
 		pctpick.o \
 		pcreset.o \
		pctorque.o \
		nhe3.o \
		les8.o

PRINT	      = pr

PROGRAM       = pct

SHELL	      = /bin/sh

SRCS	      = pct.f \
		pcnewton.f \
 		pctgam.f \
 		pctpset.f \
 		pctptim.f \
 		pctflux.f \
 		pcterr.f \
 		pctresa.f \
 		pctresb.f \
 		pctresc.f \
 		pctresd.f \
 		pctpick.f \
 		pcreset.f \
		pctorque.f \
		nhe3.f \
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

