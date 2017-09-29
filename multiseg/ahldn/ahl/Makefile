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

OBJS	      = ahl.o \
 		ahlnewt.o \
 		ahlgam.o \
 		ahlpset.o \
 		ahlptim.o \
 		ahlflux.o \
 		ahlerr.o \
 		ahlresa.o \
 		ahlresb.o \
 		ahlresc.o \
 		ahlresd.o \
 		ahlresf.o \
 		ahlpick.o \
 		ahlreset.o \
		hksub.o \
		tsc.o \
		nhe3.o \
		nkcc.o \
		kcc.o \
		les8.o

PRINT	      = pr

PROGRAM       = ahl

SHELL	      = /bin/sh

SRCS	      = ahl.f \
 		ahlnewt.f \
 		ahlgam.f \
 		ahlpset.f \
 		ahlptim.f \
 		ahlflux.f \
 		ahlerr.f \
 		ahlresa.f \
 		ahlresb.f \
 		ahlresc.f \
 		ahlresd.f \
 		ahlresf.f \
 		ahlpick.f \
 		ahlreset.f \
		hksub.f \
		tsc.f \
		nhe3.f \
		nkcc.f \
		kcc.f \
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
