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

OBJS	      = aresul.o \
		bresul.o \
		cresul.o \
		fresul.o \
		comp.o \
		newton.o \
		errvec.o \
		hksub.o \
		gamset.o \
		les8.o \
		nhe3.o \
		pset.o \
		ptim.o \
		reset.o \
		tsc.o 

PRINT	      = pr

PROGRAM       = comp

SHELL	      = /bin/sh

SRCS	      = aresul.f \
		bresul.f \
		cresul.f \
		fresul.f \
		comp.f \
		newton.f \
		errvec.f \
		hksub.f \
		gamset.f \
		les8.f \
		nhe3.f \
		pset.f \
		ptim.f \
		reset.f \
		tsc.f 

all:		$(PROGRAM)

$(PROGRAM):     $(OBJS) $(LIBS)
		@echo "Linking $(PROGRAM) ..."
		@$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)

$(OBJS):	$(SRCS)
		f90 +io77 $(FFLAGS) $(SRCS)      
		@echo "done"

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
