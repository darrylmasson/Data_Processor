CC = g++
ROOT = $(shell root-config --cflags --libs)
OBJDIR = obj
SRCDIR = src
CFLAGS = -g -Wall -Iinc -O2 -std=c++11
DEFS = -DWORKING_DIR=\"$(shell pwd)\"
OUTDIR = $(shell pwd)/../
INSTALL = -o $(OUTDIR)NG_dp
TEST = -o test_exe
SRCS = Digitizer.cpp \
	Event.cpp \
	Event_ave.cpp \
	CCM.cpp \
	DFT.cpp \
	XSQ.cpp \
	LAP.cpp \
	Processor.cpp \
	NG_dp.cpp
OBJS = $(SRCS:.cpp=.o)
VPATH = src

test :
	$(CC) $(CFLAGS) $(DEFS) $(TEST) $(addprefix $(SRCDIR)/,$(SRCS)) $(ROOT)

install : $(L)$(OBJS)
	mv $(OUTDIR)NG_dp $(OUTDIR)NG_dp_prev
	$(CC) $(CFLAGS) $(DEFS) $(INSTALL) $(addprefix $(OBJDIR)/,$(OBJS)) $(ROOT)

$(L)%.o : %.cpp
	$(CC) $(CFLAGS) $(DEFS) -c $< -o $(OBJDIR)/$@ $(ROOT)

.PHONY: clean

clean:
	rm -f $(OBJDIR)/*.o test_exe ../NG_DP
