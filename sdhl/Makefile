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

OBJS	      = sdhl.o \
 		sdhlnewt.o \
 		sdhlgam.o \
 		sdhlpset.o \
 		sdhlptim.o \
 		sdhlflux.o \
 		sdhlerr.o \
 		sdhlresa.o \
 		sdhlresb.o \
 		sdhlresc.o \
 		sdhlresd.o \
 		sdhlresf.o \
 		sdhlpick.o \
 		sdhlreset.o \
		les8.o

PRINT	      = pr

PROGRAM       = sdhl

SHELL	      = /bin/sh

SRCS	      = sdhl.f \
 		sdhlnewt.f \
 		sdhlgam.f \
 		sdhlpset.f \
 		sdhlptim.f \
 		sdhlflux.f \
 		sdhlerr.f \
 		sdhlresa.f \
 		sdhlresb.f \
 		sdhlresc.f \
 		sdhlresd.f \
 		sdhlresf.f \
 		sdhlpick.f \
 		sdhlreset.f \
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
