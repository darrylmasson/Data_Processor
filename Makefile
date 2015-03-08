CC = g++
ROOT = $(shell root-config --cflags --libs)
OBJDIR = obj
CFLAGS = -g -Wall -Iinc
OUTDIR = /data/NeutronGenerator
EXE = -o $(OUTDIR)/NG_dp
SRCS = Digitizer.cpp \
			 Event.cpp \
			 CCM.cpp \
			 DFT.cpp \
			 XSQ.cpp \
			 LAP.cpp \
			 Config.cpp \
			 Processor.cpp \
			 NG_dp.cpp
OBJS = $(SRCS:.cpp=.o)
VPATH = src

all : $(L)$(OBJS)
	$(CC) $(CFLAGS) $(EXE) $(addprefix $(OBJDIR)/,$(OBJS)) $(ROOT)

$(L)%.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $(OBJDIR)/$@ $(ROOT)

.PHONY: clean

clean:
	rm -f $(OBJDIR)/*.o
