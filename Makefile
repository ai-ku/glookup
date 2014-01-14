VERSION=2.0
CC=gcc
CFLAGS=-O3 -D_GNU_SOURCE -std=gnu99 -pedantic -Wall -save-temps
LIBS=

all: glookup glookup.txt model.txt

glookup: glookup.o dlib.o
	$(CC) $(CFLAGS) $^ $(LIBS) -o $@

glookup.o: glookup.c dlib.h
	$(CC) -c $(CFLAGS) $< -o $@

dlib.o: dlib.c dlib.h
	$(CC) -c $(CFLAGS) $< -o $@

glookup.txt: glookup.1
	groff -t -e -mandoc -Tascii $< | col -bx > $@

model.txt: model.pl
	perldoc $< > $@

release: Makefile README glookup.c glookup.1 dlib.c dlib.h
	-rm -rf glookup-${VERSION}
	mkdir glookup-${VERSION}
	cp $^ glookup-${VERSION}
	-rm -rf release/glookup.tgz
	tar cvzf release/glookup.tgz glookup-${VERSION}

clean:
	-rm *.o *.i *.s glookup README
