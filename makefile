CC=g++
CFLAGS2 += -I$(CPP_DIR)/lib/schifra/include
CFLAGS2 += -I$(CPP_DIR)/include
#LDFLAGS=-dynamiclib
#LIB_PATH=$(shell python3 -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))")
LIB_PATH=$(shell python -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))")
LDFLAGS=-shared -L $(LIB_PATH) -lpython2.7
RM=rm -f
CPP_DIR=cpp
#INCLUDE_PATH=$(shell python3 -c "from sysconfig import get_paths; print(get_paths()['include'])")
INCLUDE_PATH=$(shell python -c "from sysconfig import get_paths; print(get_paths()['include'])")

# Get the numpy include path
NUMPY_INCLUDE_PATH=$(shell python -c "import numpy; print(numpy.get_include())")

SRC = cpp/lib/schifra/src/schifra_reed_solomon_codec_validation.cpp

CFLAGS += -fPIC -O2
CFLAGS += -I$(INCLUDE_PATH)
CFLAGS += -I$(CPP_DIR)/include
CFLAGS += -I$(CPP_DIR)/lib/schifra/include
CFLAGS += -I$(INCLUDE_PATH)/site-packages/numpy/core/include
CFLAGS += -I$(NUMPY_INCLUDE_PATH) # Add the numpy include path
#CFLAGS += -I ./venv/lib/python3.9/site-packages/numpy/core/include  
TARGET_LIB1=$(CPP_DIR)/NRpyDNAcode.so
TARGET_LIB2=$(CPP_DIR)/NRpyRS.so
SRCS1=$(CPP_DIR)/src/NRpyDNAcode.cpp
SRCS2=$(CPP_DIR)/src/NRpyRS.cpp
OBJS1=$(SRCS1:.cpp=.o)
OBJS2=$(SRCS2:.cpp=.o)

.PHONY: all
all: ${TARGET_LIB2} ${TARGET_LIB1}

validator: $(SRC)
	$(CC) $(CFLAGS2) -o $@ $^

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