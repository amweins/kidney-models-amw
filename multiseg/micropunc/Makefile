F90 = f90

FFLAGS	      = -g 

CFLAGS	      = 

DEST	      =objects

EXTHDRS	      =

HDRS	      =

INSTALL	      = install

LD	      = f90

LDFLAGS	      = 

LIBS	      =

MAKEFILE      = Makefile

OBJS	      = objects/dn.o \
 		objects/dctnewt.o \
 		objects/dctgam.o \
 		objects/dctpset.o \
 		objects/dctptim.o \
 		objects/dctflux.o \
 		objects/dcterr.o \
 		objects/dctresa.o \
 		objects/dctresb.o \
 		objects/dctresc.o \
 		objects/dctresd.o \
 		objects/dctpick.o \
 		objects/dcreset.o \
 		objects/cntnewt.o \
 		objects/cntgam.o \
 		objects/cntpset.o \
 		objects/cntptim.o \
 		objects/cntflux.o \
 		objects/cnterr.o \
 		objects/cntresa.o \
 		objects/cntresb.o \
 		objects/cntresc.o \
 		objects/cntresd.o \
 		objects/cntpick.o \
 		objects/cnreset.o \
		objects/hksub.o \
		objects/ae1.o \
		objects/nhe3.o \
		objects/tsc.o \
		objects/les8.o \
		objects/eqlib.o

PRINT	      = pr

PROGRAM       = dn

SHELL	      = /usr/bin/sh

SRCS	      = dn.f \
 		dct/dctnewt.f \
 		dct/dctgam.f \
 		dct/dctpset.f \
 		dct/dctptim.f \
 		dct/dctflux.f \
 		dct/dcterr.f \
 		dct/dctresa.f \
 		dct/dctresb.f \
 		dct/dctresc.f \
 		dct/dctresd.f \
 		dct/dctpick.f \
 		dct/dcreset.f \
 		cnt/cntnewt.f \
 		cnt/cntgam.f \
 		cnt/cntpset.f \
 		cnt/cntptim.f \
 		cnt/cntflux.f \
 		cnt/cnterr.f \
 		cnt/cntresa.f \
 		cnt/cntresb.f \
 		cnt/cntresc.f \
 		cnt/cntresd.f \
 		cnt/cntpick.f \
 		cnt/cnreset.f \
		hksub.f \
		ae1.f \
		nhe3.f \
		tsc.f \
		les8.f \
		eqlib.f

all:		$(PROGRAM)

$(PROGRAM):     $(OBJS) $(LIBS)
		@echo "Linking $(PROGRAM) ..."
		@$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)
		@echo "done"

$(OBJS):	$(SRCS)
		f90 $(FFLAGS) $(SRCS)      

clean:;		@rm -f $(OBJS) core

clobber:;	@rm -f $(OBJS) $(PROGRAM) core tags

depend:;	@mkmf -f $(MAKEFILE) ROOT=$(ROOT)

echo:;		@echo $(HDRS) $(SRCS)

index:;		@ctags -wx $(HDRS) $(SRCS)

install:	$(PROGRAM)
		@echo Installing $(PROGRAM) in $(DEST)
		@-strip $(PROGRAM)
		@if [ $(DEST) != . ]; then \
		(rm -f $(DEST)/$(PROGRAM); $(INSTALL) -f $(DEST) $(PROGRAM)); fi

print:;		@$(PRINT) $(HDRS) $(SRCS)

tags:           $(HDRS) $(SRCS); @ctags $(HDRS) $(SRCS)

update:		$(DEST)/$(PROGRAM)

$(DEST)/$(PROGRAM): $(SRCS) $(LIBS) $(HDRS) $(EXTHDRS)
		@$(MAKE) -f $(MAKEFILE) ROOT=$(ROOT) DEST=$(DEST) install

.f.o:
	$(F90) $(FFLAGS) $*.f
