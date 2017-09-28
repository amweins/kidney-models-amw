F90 = f90

FFLAGS	      = -g

CFLAGS	      = 

DEST	      = .

EXTHDRS	      =

HDRS	      =

INSTALL	      = install

LD	      = f90 

LDFLAGS	      = 

LIBS	      =

MAKEFILE      = Makefile

OBJS	      = dct.o \
 		dctgam.o \
 		dctpset.o \
 		dctptim.o \
 		dctflux.o \
 		dcterr.o \
 		dctresa.o \
 		dctresb.o \
 		dctresc.o \
 		dctresd.o \
 		dctpick.o \
 		dcreset.o \
 		dctorque.o \
		hksub.o \
		tsc.o \
		nhe3.o \
		les8.o

PRINT	      = pr

PROGRAM       = dct

SHELL	      = /bin/sh

SRCS	      = dct.f \
 		dctgam.f \
 		dctpset.f \
 		dctptim.f \
 		dctflux.f \
 		dcterr.f \
 		dctresa.f \
 		dctresb.f \
 		dctresc.f \
 		dctresd.f \
 		dctpick.f \
 		dcreset.f \
 		dctorque.f \
		hksub.f \
		tsc.f \
		nhe3.f \
		les8.f

all:		$(PROGRAM)

$(PROGRAM):     $(OBJS) $(LIBS)
		@echo "Linking $(PROGRAM) ..."
		@$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)
		@echo "done"

$(OBJS):	$(SRCS)
		f90 +io77 $(FFLAGS) $(SRCS)      

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
