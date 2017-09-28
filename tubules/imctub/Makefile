FFLAGS	      = -g

CFLAGS	      = 

DEST	      = .

EXTHDRS	      =

HDRS	      =

INSTALL	      = install

LD	      = f77 

LDFLAGS	      = 

LIBS	      =

MAKEFILE      = Makefile

OBJS	      = imct.o \
 		imctgam.o \
 		imctpset.o \
 		imctptim.o \
 		imctflux.o \
 		imcterr.o \
 		imctresa.o \
 		imctresb.o \
 		imctresc.o \
 		imctresd.o \
 		imctpick.o \
 		imcreset.o \
		hksub.o \
		les8.o

PRINT	      = pr

PROGRAM       = imct

SHELL	      = /bin/sh

SRCS	      = imct.f \
 		imctgam.f \
 		imctpset.f \
 		imctptim.f \
 		imctflux.f \
 		imcterr.f \
 		imctresa.f \
 		imctresb.f \
 		imctresc.f \
 		imctresd.f \
 		imctpick.f \
 		imcreset.f \
		hksub.f \
		les8.f

all:		$(PROGRAM)

$(PROGRAM):     $(OBJS) $(LIBS)
		@echo "Linking $(PROGRAM) ..."
		@$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)
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

