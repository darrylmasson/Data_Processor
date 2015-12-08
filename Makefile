CC = g++
ROOT = $(shell root-config --cflags --libs)
OBJDIR = obj
SRCDIR = src
INCDIR = inc
CFLAGS = -g -Wall -Iinc -std=c++11
CPPFLAGS = $(CFLAGS) $(ROOT) $(DEFS)
DEFS = -DWORKING_DIR=\"$(WORKDIR)\" -DCONFIG_DIR=\"$(CONFDIR)\"
#DEFS := $(DEFS) -DCCM_ONLY
INSTALL = -o $(WORKDIR)/NG_dp
OPT = -O2
TEST = -o test_exe
sources := $(wildcard src/*.cpp)
objects := $(sources:.cpp=.o)
VPATH = src:inc
include Paths.conf
include $(sources:.cpp=.d)

test :
	$(CC) $(CPPFLAGS) $(TEST) $(sources)

install : $(objects)
	$(CC) $(CPPFLAGS) $(OPT) $(INSTALL) $(objects)

$(L)%.o : %.cpp %.h %.d
	$(CC) $(CPPFLAGS) $(OPT) -c $< -o $@

$(L)%.d : %.cpp %.h
	$(CC) -MM $(CPPFLAGS) $< -o $@

.PHONY: clean

clean:
	-rm -f $(objects)
