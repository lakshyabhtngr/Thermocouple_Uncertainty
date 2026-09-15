INCLUDES	=	-I. -I./include 

LDLIBS		=	-lm $(VXI_LIBS)

OBJS		=	ht.o lstsqr.o

SRCS		=	ht.c lstsqr.c

SOURCES		=	$(SRCS)

CFLAGS          =       -g $(CDEBUGFLAGS)

PROGS		=	ht

all::  $(PROGS)

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDES) $<
