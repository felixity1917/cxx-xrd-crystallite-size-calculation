# Compiler and flags
CXX := g++
CXXFLAGS := -std=c++20 -Wall -Wextra -O2 -Iinclude
LDFLAGS :=

# Directory structure
SRC_DIR := src
INC_DIR := include
OBJ_DIR := build
BIN_DIR := bin
DIST_DIR := dist

# Project name
PROJECT_NAME := cxx-xrd-crystallite-size-calculation

# Files
TARGET := $(BIN_DIR)/$(PROJECT_NAME)
SRCS := $(wildcard $(SRC_DIR)/*.cxx)
OBJS := $(patsubst $(SRC_DIR)/%.cxx, $(OBJ_DIR)/%.o, $(SRCS))

# Default target
all: $(TARGET)

# Link the final binary
$(TARGET): $(OBJS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Compile source files to object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cxx | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Ensure bin and build directories exist
$(BIN_DIR) $(OBJ_DIR) $(DIST_DIR):
	mkdir -p $@

# Package the binary and headers into a .tar.gz for distribution
dist: all | $(DIST_DIR)
	tar -czf $(DIST_DIR)/$(PROJECT_NAME).tar.gz -C $(BIN_DIR) $(PROJECT_NAME) -C ../include .

# Clean build and output files
clean:
	rm -rf $(OBJ_DIR)/*.o $(BIN_DIR)/$(PROJECT_NAME)

# Clean everything including dist
dist-clean: clean
	rm -rf $(DIST_DIR)/*.tar.gz
