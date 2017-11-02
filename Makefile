include make.inc

LIB = libfsi

all: example_

example_: $(LIB)
	(cd example; $(MAKE))

libfsi:
	(cd src; $(MAKE))
clean:
	(cd src; $(MAKE) clean)
	(cd example; $(MAKE) clean)
