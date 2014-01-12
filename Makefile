CC=gcc
CFLAGS=-O3 -Wall -std=gnu99 -I. `pkg-config --cflags glib-2.0`
LIBS=`pkg-config --libs glib-2.0`
VERSION=1.12

ifdef CYGWIN
CFLAGS=-static -s -O3 -Wall -std=c99 -I. `pkg-config --cflags glib-2.0`
endif

glookup: glookup.o
	$(CC) $(CFLAGS) $^ $(LIBS) -o $@

glookup.o: glookup.c foreach.h procinfo.h ghashx.h
	$(CC) -c $(CFLAGS) $< -o $@

glookup.txt: glookup.1
	groff -t -e -mandoc -Tascii $^ | col -bx > $@

release: Makefile README glookup.c glookup.1 glookup.txt foreach.h ghashx.h procinfo.h
	-rm -rf glookup-${VERSION}
	mkdir glookup-${VERSION}
	cp $^ glookup-${VERSION}
	-rm -rf release/glookup.tgz
	tar cvzf release/glookup.tgz glookup-${VERSION}
