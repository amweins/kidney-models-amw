F90 = f90

FFLAGS	      = -g +io77


CFLAGS	      = 

DEST	      = .

EXTHDRS	      =

HDRS	      =

INSTALL	      = install

LD	      = f90 

LDFLAGS	      = 

LIBS	      =

MAKEFILE      = Makefile

OBJS	      = cnt.o \
 		cntgam.o \
 		cntpset.o \
 		cntptim.o \
 		cntflux.o \
 		cnterr.o \
 		cntresa.o \
 		cntresb.o \
 		cntresc.o \
 		cntresd.o \
 		cntpick.o \
 		cntorque.o \
 		cnreset.o \
		hksub.o \
		ae1.o \
		les8.o

PRINT	      = pr

PROGRAM       = cnt

SHELL	      = /bin/sh

SRCS	      = cnt.f \
 		cntgam.f \
 		cntpset.f \
 		cntptim.f \
 		cntflux.f \
 		cnterr.f \
 		cntresa.f \
 		cntresb.f \
 		cntresc.f \
 		cntresd.f \
 		cntpick.f \
 		cntorque.f \
 		cnreset.f \
		hksub.f \
		ae1.f \
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
