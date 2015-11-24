CC = g++
ROOT = $(shell root-config --cflags --libs)
OBJDIR = obj
SRCDIR = src
INCDIR = inc
CFLAGS = -g -Wall -Iinc -std=c++11
DEFS = -DWORKING_DIR=\"$(WORKDIR)\" -DCONFIG_DIR=\"$(CONFDIR)\"
#DEFS := $(DEFS) -DCCM_ONLY
INSTALL = -o $(WORKDIR)/NG_dp
OPT = -O2
TEST = -o test_exe
sources := $(wildcard src/*.cpp)
objects := $(sources:.cpp=.o)
VPATH = src
include Paths.conf

test :
	$(CC) $(CFLAGS) $(DEFS) $(TEST) $(sources) $(ROOT)

install : $(objects)
	$(CC) $(CFLAGS) $(OPT) $(DEFS) $(INSTALL) $(objects) $(ROOT)

%.o : %.cpp %.d
	$(CC) $(CFLAGS) $(OPT) $(DEFS) -c $< -o $@ $(ROOT)

%.d : %.cpp
	$(CC) -MM $(CFLAGS) $(ROOT) $< -o $@

.PHONY: clean

clean:
	rm -f $(OBJDIR)/*.o test_exe ../NG_DP
