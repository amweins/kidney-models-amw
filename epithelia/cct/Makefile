F95 = gfortran

FFLAGS	      = -g -w

CFLAGS	      = 

DEST	      = .

EXTHDRS	      =

HDRS	      =

INSTALL	      = install

LD	      = gfortran

LDFLAGS	      = 

LIBS	      =

MAKEFILE      = Makefile

OBJS	      = ae1.o \
		aresul.o \
		bresul.o \
		cresul.o \
		fresul.o \
		comp.o \
		errvec.o \
		hksub.o \
		gamset.o \
		les8.o \
		pset.o \
		ptim.o \
		reset.o

PRINT	      = pr

PROGRAM       = comp

SHELL	      = /bin/sh

SRCS	      = ae1.f \
		aresul.f \
		bresul.f \
		cresul.f \
		fresul.f \
		comp.f \
		errvec.f \
		hksub.f \
		gamset.f \
		les8.f \
		pset.f \
		ptim.f \
		reset.f

all:		$(PROGRAM)

$(PROGRAM):     $(OBJS) $(LIBS)
		@echo "Linking $(PROGRAM) ..."
		@$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)
		@echo "done"

$(OBJS):	$(SRCS)
		gfortran -c $(FFLAGS) $(SRCS)

clean:;		@rm -f $(OBJS) core

clobber:;	@rm -f $(OBJS) $(PROGRAM) core tags

install:	$(PROGRAM)
		echo Installing $(PROGRAM) in $(DEST)
		if [ $(DEST) != . ]; then \
		(rm -f $(DEST)/$(PROGRAM); $(INSTALL) -f $(DEST) $(PROGRAM)); fi


