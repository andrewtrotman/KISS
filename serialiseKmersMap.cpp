/*
	SERIALISEKMERSMAP.CPP
	---------------------
	indexReference

	Created by Shlomo Geva on 19/7/2023.
*/
#include <fstream>
#include <iostream>

#include "serialiseKmersMap.hpp"

void serializeMap(const std::vector<std::vector<uint32_t>>& kmersMap, const std::string& innerMapFilename, const std::string& outerMapFilename) {
    // Set the buffer size to 8192 bytes (for example)
    constexpr std::streamsize bufferSize = 8192;

    // Open the files with custom buffer sizes
    std::ofstream innerMapFile(innerMapFilename, std::ios::binary);
    innerMapFile.rdbuf()->pubsetbuf(nullptr, bufferSize);

    std::ofstream outerMapFile(outerMapFilename, std::ios::binary);
    outerMapFile.rdbuf()->pubsetbuf(nullptr, bufferSize);

    uint32_t offset = 0;
    for (const auto& innerVector : kmersMap) {
        // Write inner vector data
        innerMapFile.write(reinterpret_cast<const char*>(innerVector.data()), innerVector.size() * sizeof(uint32_t));

        // Write the offset for the outer map
        outerMapFile.write(reinterpret_cast<const char*>(&offset), sizeof(uint32_t));

        // Calculate the offset for the next inner vector
        offset += static_cast<uint32_t>(innerVector.size());
    }

    innerMapFile.close();
    outerMapFile.close();
}

void deserializeMap(const std::string& innerMapFilename, const std::string& outerMapFilename, std::vector<uint32_t>& innerMapBlob, std::vector<uint32_t>& outerMapBlob) {
    std::ifstream innerMapFile(innerMapFilename, std::ios::binary);
    std::ifstream outerMapFile(outerMapFilename, std::ios::binary);

    // Get the size of the files to determine the required memory
    innerMapFile.seekg(0, std::ios::end);
    size_t innerBlobSize = innerMapFile.tellg();
    innerMapFile.seekg(0, std::ios::beg);

    outerMapFile.seekg(0, std::ios::end);
    size_t outerBlobSize = outerMapFile.tellg();
    outerMapFile.seekg(0, std::ios::beg);

    // Read the blobs into memory
    innerMapBlob.resize(innerBlobSize / sizeof(uint32_t));
    innerMapFile.read(reinterpret_cast<char*>(innerMapBlob.data()), innerBlobSize);

    outerMapBlob.resize(outerBlobSize / sizeof(uint32_t));
    outerMapFile.read(reinterpret_cast<char*>(outerMapBlob.data()), outerBlobSize);

    innerMapFile.close();
    outerMapFile.close();
}

// Helper function to access the index
std::vector<uint32_t> getInnerVector(const std::vector<uint32_t>& innerMapBlob, const std::vector<uint32_t>& outerMapBlob, size_t index) {
    std::vector<uint32_t> innerVector;

    if (index < outerMapBlob.size()) {
        size_t startOffset = outerMapBlob[index];
        size_t endOffset = (index + 1 < outerMapBlob.size()) ? outerMapBlob[index + 1] : innerMapBlob.size();
        innerVector = std::vector<uint32_t>(innerMapBlob.begin() + startOffset, innerMapBlob.begin() + endOffset);
    }

    return innerVector;
}

bool writeTextBlobToFile(const char* text, std::size_t length, const std::string& filename) {
    std::ofstream outputFile(filename, std::ios::binary);
    if (!outputFile.is_open()) {
        return false; // Failed to open the file
    }

    outputFile.write(text, length);

    outputFile.close();

    return true; // File write successful
}

std::pair<char*, std::size_t> readTextBlobFromFile(const std::string& filename) {
    std::ifstream inputFile(filename, std::ios::binary | std::ios::ate);
    if (!inputFile.is_open()) {
        return {nullptr, 0}; // Failed to open the file
    }

    std::size_t length = inputFile.tellg(); // Get the file size
    inputFile.seekg(0); // Move the file pointer back to the beginning

    std::vector<char> buffer(length); // Create a buffer to hold the content

    if (!inputFile.read(buffer.data(), length)) {
        return {nullptr, 0}; // Failed to read the content
    }

    inputFile.close();

    // Allocate memory for the text blob and copy the content from the buffer
    char* textBlob = new char[length];
    std::memcpy(textBlob, buffer.data(), length);

    return {textBlob, length}; // Return the pointer to the blob and its size
}

