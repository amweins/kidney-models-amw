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

OBJS	      = pct.o \
 		pctgam.o \
 		pctpset.o \
 		pctptim.o \
 		pctflux.o \
 		pcterr.o \
 		pctresa.o \
 		pctresb.o \
 		pctresc.o \
 		pctresd.o \
 		pctresf.o \
 		pctpick.o \
 		pcreset.o \
 		pctqo2.o \
		nhe3.o \
		les8.o

PRINT	      = pr

PROGRAM       = pct

SHELL	      = /bin/sh

SRCS	      = pct.f \
 		pctgam.f \
 		pctpset.f \
 		pctptim.f \
 		pctflux.f \
 		pcterr.f \
 		pctresa.f \
 		pctresb.f \
 		pctresc.f \
 		pctresd.f \
 		pctresf.f \
 		pctpick.f \
 		pcreset.f \
 		pctqo2.f \
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

