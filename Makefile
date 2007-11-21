CC=gcc
CFLAGS=-O2 -Wall -std=c99 -I. `pkg-config --cflags glib-2.0`
LIBS=`pkg-config --libs glib-2.0`

ifdef CYGWIN
CFLAGS=-static -s -O3 -Wall -std=c99 -I. `pkg-config --cflags glib-2.0`
endif

glookup: glookup.o 
	$(CC) $(CFLAGS) $^ $(LIBS) -o $@

glookup.o: glookup.c foreach.h procinfo.h ghashx.h
	$(CC) -c $(CFLAGS) $< -o $@

