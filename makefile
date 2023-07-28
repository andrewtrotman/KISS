# Compiler and flags
CC = g++
CFLAGS = -O3 -std=c++11
#CFLAGS = -g -o -std=c++11

# Source directory and files
SOURCE_DIR = .
SOURCES = main.cpp encode_kmer_2bit.cpp hash.cpp indexGenome.cpp serialiseKmersMap.cpp packGenomeBlob.cpp

# Header directory
HEADER_DIR = headers

# Object directory and files
OBJECT_DIR = objects
OBJECTS = $(SOURCES:%.cpp=$(OBJECT_DIR)/%.o)

# Executable name
EXECUTABLE = indexReference

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $@

$(OBJECTS): $(OBJECT_DIR)/%.o: $(SOURCE_DIR)/%.cpp
	$(CC) $(CFLAGS) -I$(HEADER_DIR) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

