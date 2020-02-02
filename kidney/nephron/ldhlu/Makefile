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

OBJS	      = ldhlu.o \
 		ldhlunewt.o \
 		ldhlugam.o \
 		ldhlupset.o \
 		ldhluptim.o \
 		ldhluflux.o \
 		ldhluerr.o \
 		ldhluresa.o \
 		ldhluresb.o \
 		ldhluresc.o \
 		ldhluresd.o \
 		ldhluresf.o \
 		ldhlupick.o \
 		ldhlureset.o \
		les8.o

PRINT	      = pr

PROGRAM       = ldhlu

SHELL	      = /bin/sh

SRCS	      = ldhlu.f \
 		ldhlunewt.f \
 		ldhlugam.f \
 		ldhlupset.f \
 		ldhluptim.f \
 		ldhluflux.f \
 		ldhluerr.f \
 		ldhluresa.f \
 		ldhluresb.f \
 		ldhluresc.f \
 		ldhluresd.f \
 		ldhluresf.f \
 		ldhlupick.f \
 		ldhlureset.f \
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
