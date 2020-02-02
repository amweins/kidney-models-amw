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

OBJS	      = tahl.o \
 		tahlnewt.o \
 		tahlgam.o \
 		tahlpset.o \
 		tahlptim.o \
 		tahlflux.o \
 		tahlerr.o \
 		tahlresa.o \
 		tahlresb.o \
 		tahlresc.o \
 		tahlresd.o \
 		tahlresf.o \
 		tahlpick.o \
 		tahlreset.o \
		les8.o

PRINT	      = pr

PROGRAM       = tahl

SHELL	      = /bin/sh

SRCS	      = tahl.f \
 		tahlnewt.f \
 		tahlgam.f \
 		tahlpset.f \
 		tahlptim.f \
 		tahlflux.f \
 		tahlerr.f \
 		tahlresa.f \
 		tahlresb.f \
 		tahlresc.f \
 		tahlresd.f \
 		tahlresf.f \
 		tahlpick.f \
 		tahlreset.f \
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
