PYTHON_VERSION ?= python3

PYTHON_CFLAGS  := $(shell $(PYTHON_VERSION)-config --includes)
PYTHON_LDFLAGS := $(shell $(PYTHON_VERSION)-config --embed --ldflags 2>/dev/null || $(PYTHON_VERSION)-config --ldflags)

CXX = g++
CXXFLAGS = -Wall -std=c++11 -Iinclude -DWITHOUT_NUMPY $(PYTHON_CFLAGS)
LDFLAGS  = $(PYTHON_LDFLAGS)

SRC_DIR = src
BUILD_DIR = build

SRC = $(wildcard $(SRC_DIR)/*.cpp)
OBJ = $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SRC))

tsp-study: $(OBJ)
	$(CXX) $(OBJ) -o $@ $(LDFLAGS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR) *.o tsp-study