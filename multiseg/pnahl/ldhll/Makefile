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

OBJS	      = ldhll.o \
 		ldhllnewt.o \
 		ldhllgam.o \
 		ldhllpset.o \
 		ldhllptim.o \
 		ldhllflux.o \
 		ldhllerr.o \
 		ldhllresa.o \
 		ldhllresb.o \
 		ldhllresc.o \
 		ldhllresd.o \
 		ldhllresf.o \
 		ldhllpick.o \
 		ldhllreset.o \
		les8.o

PRINT	      = pr

PROGRAM       = ldhll

SHELL	      = /bin/sh

SRCS	      = ldhll.f \
 		ldhllnewt.f \
 		ldhllgam.f \
 		ldhllpset.f \
 		ldhllptim.f \
 		ldhllflux.f \
 		ldhllerr.f \
 		ldhllresa.f \
 		ldhllresb.f \
 		ldhllresc.f \
 		ldhllresd.f \
 		ldhllresf.f \
 		ldhllpick.f \
 		ldhllreset.f \
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
