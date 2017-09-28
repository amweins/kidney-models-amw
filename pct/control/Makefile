F95 = f95

FFLAGS	      = -g

CFLAGS	      = 

DEST	      = .

EXTHDRS	      =

HDRS	      =

INSTALL	      = install

LD	      = f95

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
		gamset.o \
		les8.o \
		nhe3.o \
		pset.o \
		ptim.o \
		reset.o

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
		gamset.f \
		les8.f \
		nhe3.f \
		pset.f \
		ptim.f \
		reset.f

$(PROGRAM):     $(OBJS) $(LIBS)
		echo "Linking $(PROGRAM) ..."
		$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)
		echo "done"

$(OBJS):	$(SRCS)
		f95 -c $(FFLAGS) $(SRCS)

install:	$(PROGRAM)
		echo Installing $(PROGRAM) in $(DEST)
		if [ $(DEST) != . ]; then \
		(rm -f $(DEST)/$(PROGRAM); $(INSTALL) -f $(DEST) $(PROGRAM)); fi

clean:;		rm -f $(OBJS) core

