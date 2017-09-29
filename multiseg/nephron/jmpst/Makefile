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

OBJS	      = jmpst.o \
 		jmpstnewt.o \
 		jmpstgam.o \
 		jmpstpset.o \
 		jmpstptim.o \
 		jmpstflux.o \
 		jmpsterr.o \
 		jmpstresa.o \
 		jmpstresb.o \
 		jmpstresc.o \
 		jmpstresd.o \
 		jmpstresf.o \
 		jmpstpick.o \
 		jmpsreset.o \
 		jmpstqo2.o \
		nhe3.o \
		les8.o

PRINT	      = pr

PROGRAM       = jmpst

SHELL	      = /bin/sh

SRCS	      = jmpst.f \
 		jmpstnewt.f \
 		jmpstgam.f \
 		jmpstpset.f \
 		jmpstptim.f \
 		jmpstflux.f \
 		jmpsterr.f \
 		jmpstresa.f \
 		jmpstresb.f \
 		jmpstresc.f \
 		jmpstresd.f \
 		jmpstresf.f \
 		jmpstpick.f \
 		jmpsreset.f \
 		jmpstqo2.f \
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
