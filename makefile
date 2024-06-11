CC=g++
CFLAGS=-fPIC -Wall -Wextra -O2
LDFLAGS=-dynamiclib
RM=rm -f
CPP_DIR=cpp
INCLUDE_PATH=$(shell python3 -c "from sysconfig import get_paths; print(get_paths()['include'])")

CFLAGS += -I$(INCLUDE_PATH)
CFLAGS += -I$(CPP_DIR)/include
CFLAGS += -I$(CPP_DIR)/lib/schifra/include
CFLAGS += -I /Users/mguyot/Documents/fork_hedges/venv/lib/python3.12/site-packages/numpy/core/include
TARGET_LIB1=$(CPP_DIR)/NRpyDNAcode.so
TARGET_LIB2=$(CPP_DIR)/NRpyRS.so
SRCS1=$(CPP_DIR)/src/NRpyDNAcode.cpp
SRCS2=$(CPP_DIR)/src/NRpyRS.cpp
OBJS1=$(SRCS1:.cpp=.o)
OBJS2=$(SRCS2:.cpp=.o)

.PHONY: all
all: ${TARGET_LIB1} ${TARGET_LIB2}

$(TARGET_LIB1): $(OBJS1)
	$(CC) $(LDFLAGS) -o $@ $^

$(TARGET_LIB2): $(OBJS2)
	$(CC) $(LDFLAGS) -o $@ $^

$(CPP_DIR)/%.o: $(CPP_DIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(TARGET_LIB1) $(OBJS1) $(TARGET_LIB2) $(OBJS2)

.PHONY: test
test:
	python python/test_program.py
