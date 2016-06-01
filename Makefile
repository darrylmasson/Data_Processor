CC = g++
ROOT = $(shell root-config --cflags --libs)
OBJDIR = obj
SRCDIR = src
INCDIR = inc
CFLAGS = -g -Wall -Iinc -std=c++11 -O2
CPPFLAGS = $(CFLAGS) $(ROOT) $(DEFS)
DEFS = -DCONFIGDIR=\"$(CONFDIR)\"
#DEFS := $(DEFS) -DCCM_ONLY
INSTALL = -o $(INSTALLDIR)/NG_dp
TEST = -o test_exe
sources := $(wildcard src/*.cpp)
objects := $(sources:.cpp=.o)
VPATH = src:inc
include Paths.conf

# if there are ROOT errors in the linking phase, use this shell command:
# $ g++ -std=c++11 -g -Wall -Iinc -O2 -o -DCONFIGDIR=\"(configdir)\" NG_dp src/*.cpp `root-config --cflags --libs`


test :
	$(CC) $(CPPFLAGS) $(TEST) $(sources)

install : $(objects)
	$(CC) $(CPPFLAGS) $(INSTALL) $(objects)

$(L)%.o : %.cpp %.h %.d
	$(CC) $(CPPFLAGS) -c $< -o $@

$(L)%.d : %.cpp %.h
	$(CC) -MM $(CPPFLAGS) $< -o $@

.PHONY: clean

clean:
	-rm -f $(objects)
