# Use the following two lines for GCC
#CC = gcc
#OPT = -O2

# Use the following two lines for SGI's CC compiler
CC = cc
OPT = -O2

# Use the following two lines for Sun's Workshop compiler
#CC = cc
#OPT = -FAST

DBG =
CFLAGS = $(DBG) $(OPT)
CCLIBS =  -lm

BBKNNOBJS = util.o bbknn.o
BBKNNSRCS = ${BBKNNOBJS:.o=.c}
BBKNNHDRS = defines.h util.h

TOOLSRCS = Nearby.c Conswat.c All.c Scale.c
TOOLEXES = ${TOOLSRCS:.c=.EXE}

TOOLSCRIPTS = allwats consolv2pdb extract pdbdehydrogen sumpred

BINDIR = ../bin
BINFILES = bbknn $(TOOLEXES) $(TOOLSCRIPTS)

all:	bbknn consolv

$(BBKNNOBJS):  $(BBKNNSRCS) $(BBKNNHDRS)

bbknn:	$(BBKNNOBJS)
	$(CC) $(CFLAGS) -o $@ $(BBKNNOBJS) $(CCLIBS) 

$(TOOLEXES):  $(TOOLSRCS)
	$(CC) $(CFLAGS) -DAHP_LIST_FILE=\"$(AHPLIST)\" \
	-DPDBLOC=\"$(PDBDEFLOC)\" -o $@ $(@:.EXE=.c) $(CCLIBS) 

consolv:  $(TOOLEXES)

install: $(BINFILES)
	cp -f $(BINFILES) $(BINDIR)

clean: 
	rm -f bbknn $(TOOLEXES) *.o
