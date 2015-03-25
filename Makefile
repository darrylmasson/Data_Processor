CC = g++
ROOT = $(shell root-config --cflags --libs)
OBJDIR = obj
SRCDIR = src
CFLAGS = -g -Wall -Iinc -O
OUTDIR = /data/NeutronGenerator
INSTALL = -o $(OUTDIR)/NGrawdp
TEST = -o test_exe
SRCS = Digitizer.cpp \
			 Event.cpp \
			 Event_ave.cpp \
			 CCM.cpp \
			 DFT.cpp \
			 XSQ.cpp \
			 LAP.cpp \
			 Config.cpp \
			 Processor.cpp \
			 NG_dp.cpp
OBJS = $(SRCS:.cpp=.o)
VPATH = src

test :
	$(CC) $(CFLAGS) $(TEST) $(addprefix $(SRCDIR)/,$(SRCS)) $(ROOT)

install : $(L)$(OBJS)
	$(CC) $(CFLAGS) $(INSTALL) $(addprefix $(OBJDIR)/,$(OBJS)) $(ROOT)

$(L)%.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $(OBJDIR)/$@ $(ROOT)

.PHONY: clean

clean:
	rm -f $(OBJDIR)/*.o test_exe ../NGrawDP
