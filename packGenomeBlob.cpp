/*
	PACKGENOMEBLOB.CPP
	------------------
	indexReference

	Created by Shlomo Geva on 22/7/2023.
*/
#include <stdint.h>
#include <string.h>
#include <sys/stat.h>

#include <map>
#include <chrono>
#include <string>
#include <iostream>

#include "packGenomeBlob.hpp"

/*
	PACKGENOME()
	------------
*/
size_t packGenome(char *genome, uint64_t genome_size, std::map<std::uint32_t, std::string> &referenceIDMap)
	{
	char *from;
	char *to;
	char *end = genome + genome_size;

	// Traverse the old genome, character by character
	for (from = to = genome; from < end; from++)
		{
		char c = *from;
		if (c == '\n' || c == 'N')
			{/* Nothing */}
		else if (c == '>')
			{
			// Find the start and end of the line
			char *start = from;
			do
				from++;
			while (*from != '\n');

			// Save the ID line to the referenceIDMap
			referenceIDMap[to - genome] = std::string(start, from - start + 1);
			}
		else
			*to++ = c;			// Copy the DNA characters to the new position
		}

	// Null-terminate the new genome explicitly to make it a valid C-string
	*to = '\0';

	// Return the new length of the modified text
	return to - genome;
	}

/*
	MAIN_UNITTEST()
	---------------
*/
int main_unittest()
	{
	// Sample usage
	char genome[] = ">ID1\nAGCT\n>NID2\nNNNN\nATGC\n";
	std::cout << "Original Genome: " << genome << std::endl;
	std::cout << "Old Genome Length: " << strlen(genome) << std::endl;
	uint64_t file_size = strlen(genome);

	std::map<std::uint32_t, std::string> referenceIDMap;
	size_t newLength = packGenome(genome, file_size, referenceIDMap);

	// Access the modified genome using the original pointer 'genome'
	std::cout << "Modified Genome  : " << genome << std::endl;
	std::cout << "New Genome Length: " << newLength << std::endl;

	for (const auto &thang : referenceIDMap)
		std::cout << thang.first << " " << thang.second;

	return 0;
	}
