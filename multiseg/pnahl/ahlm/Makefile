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

OBJS	      = ahlm.o \
 		ahlmnewt.o \
 		ahlmgam.o \
 		ahlmpset.o \
 		ahlmptim.o \
 		ahlmflux.o \
 		ahlmerr.o \
 		ahlmresa.o \
 		ahlmresb.o \
 		ahlmresc.o \
 		ahlmresd.o \
 		ahlmresf.o \
 		ahlmpick.o \
 		ahlmreset.o \
		hksub.o \
		tsc.o \
		nhe3.o \
		nkcc.o \
		kcc.o \
		les8.o

PRINT	      = pr

PROGRAM       = ahlm

SHELL	      = /bin/sh

SRCS	      = ahlm.f \
 		ahlmnewt.f \
 		ahlmgam.f \
 		ahlmpset.f \
 		ahlmptim.f \
 		ahlmflux.f \
 		ahlmerr.f \
 		ahlmresa.f \
 		ahlmresb.f \
 		ahlmresc.f \
 		ahlmresd.f \
 		ahlmresf.f \
 		ahlmpick.f \
 		ahlmreset.f \
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
