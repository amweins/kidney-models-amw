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

OBJS	      = jmpct.o \
 		jmpctnewt.o \
 		jmpctgam.o \
 		jmpctpset.o \
 		jmpctptim.o \
 		jmpctflux.o \
 		jmpcterr.o \
 		jmpctresa.o \
 		jmpctresb.o \
 		jmpctresc.o \
 		jmpctresd.o \
 		jmpctresf.o \
 		jmpctpick.o \
 		jmpcreset.o \
 		jmpctqo2.o \
		jmpctorque.o \
		nhe3.o \
		les8.o

PRINT	      = pr

PROGRAM       = jmpct

SHELL	      = /bin/sh

SRCS	      = jmpct.f \
 		jmpctnewt.f \
 		jmpctgam.f \
 		jmpctpset.f \
 		jmpctptim.f \
 		jmpctflux.f \
 		jmpcterr.f \
 		jmpctresa.f \
 		jmpctresb.f \
 		jmpctresc.f \
 		jmpctresd.f \
 		jmpctresf.f \
 		jmpctpick.f \
 		jmpcreset.f \
 		jmpctqo2.f \
		jmpctorque.f \
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
