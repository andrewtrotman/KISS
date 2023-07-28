//
//  packGenomeBlob.cpp
//  indexReference
//
//  Created by Shlomo Geva on 22/7/2023.
//

#include "packGenomeBlob.hpp"

#include <cstdint>
#include <cstring>
#include <iostream>
#include <map>
#include <string>

std::size_t packGenome(char* genome, std::map<std::uint32_t, std::string>& referenceIDMap) {
    // Helper function to skip the ID line (if any) and return the index of the next line
    auto skipIDLine = [](char* ptr) {
        std::size_t index = 0;
        while (ptr[index] != '\n') {
            ++index;
        }
        return index + 1; // Move past '\n'
    };

    // Initialize pointers for the old and new genomes
    char* oldPtr = genome;
    char* newPtr = genome;

    // Initialize variables to keep track of the new genome length and the current ID start position
    std::uint32_t newLength = 0;
    std::uint32_t currentIDStart = 0;

    // Traverse the old genome, character by character
    for (std::size_t i = 0; oldPtr[i] != '\0'; ++i) {
        char c = oldPtr[i];

        // Check for the ID line (beginning with '>')
        if (c == '>') {
            // Save the ID line to the referenceIDMap and update the currentIDStart
            referenceIDMap[newLength] = std::string(oldPtr + i, skipIDLine(oldPtr + i));
            currentIDStart = newLength;
            i += skipIDLine(oldPtr + i) - 1; // Move to the end of the ID line
            continue;
        }

        // Skip newline characters and 'N' characters
        if (c == '\n' || c == 'N') {
            continue;
        }

        // Copy the DNA characters to the new position
        newPtr[newLength] = c;
        ++newLength;
    }

    // Null-terminate the new genome explicitly to make it a valid C-string
    newPtr[newLength] = '\0';

    // Return the new length of the modified text
    return newLength;
}

int mainTest() {
    // Sample usage
    char genome[] = ">ID1\nAGCT\n>NID2\nNNNN\nATGC\n";

    std::map<std::uint32_t, std::string> referenceIDMap;

    std::size_t newLength = packGenome(genome, referenceIDMap);

    // Access the modified genome using the original pointer 'genome'
    std::cout << "Modified Genome: " << genome << std::endl;
    std::cout << "New Genome Length: " << newLength << std::endl;

    // The original pointer 'genome' is still valid, and you can use it as needed.

    // Free memory (if applicable) using the original pointer 'genome'
    // ...

    return 0;
}

