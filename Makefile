
# Makefile for AFFCK 2.1

OST := $(shell uname -s)
ifeq ($(OST),Darwin)
	CXX = g++-5
else
	CXX = g++
endif
MKDIR = mkdir
RM = rm -fr
CFLAGS = -O3
SRC_PATH = ./src
BIN_PATH = ./bin
SRCS = $(wildcard $(SRC_PATH)/*.cpp)
OBJS = $(patsubst %.cpp, %.o, $(SRCS))
TARGET = $(BIN_PATH)/affck.exe

all : $(BIN_PATH) $(TARGET)

$(BIN_PATH) :
	$(MKDIR) $@

$(TARGET) : $(OBJS)
	$(CXX) $(CFLAGS) $^ -o $@
	@echo "compile completed."

.cpp.o:
	$(CXX) $(CFLAGS) -c -o $@ $<

clean:
	$(RM) $(OBJS) $(TARGET) $(BIN_PATH)
	@echo "ok."

.PHONY : all clean
