TOP := $(CURDIR)/..

include $(TOP)/Makefile.compiler
include $(TOP)/Makefile.include

LIBRARY = wrapper.a
FILES = wrapper.o

.PHONY: all clean

all: $(LIBRARY)

Makefile.depend:

clean:
	rm -f *.o $(LIBRARY)

$(LIBRARY): $(FILES)

include Makefile.depend
