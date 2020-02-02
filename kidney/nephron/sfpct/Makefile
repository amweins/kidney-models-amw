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

OBJS	      = sfpct.o \
 		sfpctnewt.o \
 		sfpctgam.o \
 		sfpctpset.o \
 		sfpctptim.o \
 		sfpctflux.o \
 		sfpcterr.o \
 		sfpctresa.o \
 		sfpctresb.o \
 		sfpctresc.o \
 		sfpctresd.o \
 		sfpctresf.o \
 		sfpctpick.o \
 		sfpcreset.o \
 		sfpctqo2.o \
		sfpctorque.o \
		nhe3.o \
		les8.o

PRINT	      = pr

PROGRAM       = sfpct

SHELL	      = /bin/sh

SRCS	      = sfpct.f \
 		sfpctnewt.f \
 		sfpctgam.f \
 		sfpctpset.f \
 		sfpctptim.f \
 		sfpctflux.f \
 		sfpcterr.f \
 		sfpctresa.f \
 		sfpctresb.f \
 		sfpctresc.f \
 		sfpctresd.f \
 		sfpctresf.f \
 		sfpctpick.f \
 		sfpcreset.f \
 		sfpctqo2.f \
		sfpctorque.f \
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
