.PHONY: all tests clean release remake

GTEST_DIR := $(HOME)/development_cpp/gtest/googletest

GTEST_HEADERS := $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h
GTEST_SRCS_ := $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

CC := gcc-9
CXX := g++-9
RM := rm -f
TARGET := qinfer

SRC := src
INC := -I include
OBJ := obj
BIN := bin
UT := unittest

SRCS = $(shell find $(SRC) -name '*.cpp')
UT_SRCS = $(shell find $(UT)/$(SRC) -name '*.cpp')
OBJS = $(subst $(SRC)/,$(OBJ)/,$(subst .cpp,.o,$(SRCS)))
UT_OBJS = $(subst $(SRC)/,$(OBJ)/,$(subst .cpp,.o,$(UT_SRCS)))

ifeq ($(RELEASE),true)
	CPPFLAGS = -O3 -DNDEBUG -Wall -Wextra -pthread -std=c++17
	LDFLAGS = -s -lstdc++fs
else
	CPPFLAGS = -g -O0 -Wall -Wextra -pthread -std=c++17
	LDFLAGS = -lstdc++fs
endif

all: $(BIN)/$(TARGET)

release:
	make RELEASE=true

tests: $(BIN)/unittest
	./$(BIN)/unittest

$(BIN)/$(TARGET): $(OBJS)
	mkdir -p $(BIN)
	$(CXX) $(LDFLAGS) -o $@ $(OBJS)

$(BIN)/unittest : $(filter-out $(OBJ)/main.o,$(OBJS)) $(UT_OBJS) $(BIN)/gtest_main.a
	mkdir -p $(BIN)
	$(CXX) -isystem $(GTEST_DIR)/include $(CPPFLAGS) -lpthread $^ -o $@

$(OBJ)/%.o: $(SRC)/%.cpp
	mkdir -p $(OBJ)
	$(CXX) $(CPPFLAGS) $(INC) -c $< -o $@

$(UT)/$(OBJ)/%.o: $(UT)/$(SRC)/%.cpp
	mkdir -p $(UT)/$(OBJ)
	$(CXX) $(CPPFLAGS) -isystem $(GTEST_DIR)/include $(INC) -c $< -o $@

$(UT)/$(OBJ)/gtest-all.o : $(GTEST_SRCS_)
	$(CXX) -isystem $(GTEST_DIR)/include -I$(GTEST_DIR) $(CPPFLAGS) -c \
            $(GTEST_DIR)/src/gtest-all.cc -o $@

$(UT)/$(OBJ)/gtest_main.o : $(GTEST_SRCS_)
	$(CXX) -isystem $(GTEST_DIR)/include -I$(GTEST_DIR) $(CPPFLAGS) -c \
            $(GTEST_DIR)/src/gtest_main.cc -o $@

$(UT)/$(OBJ)/gtest.a : gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

$(BIN)/gtest_main.a : $(UT)/$(OBJ)/gtest-all.o $(UT)/$(OBJ)/gtest_main.o
	mkdir -p $(BIN)
	$(AR) $(ARFLAGS) $@ $^

remake: clean all

clean:
	rm -rf $(OBJ)
	rm -rf $(BIN)
	rm -rf $(UT)/$(OBJ)
