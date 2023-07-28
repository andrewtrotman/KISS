/*
	INDEXGENOME.HPP
	---------------
	indexKISS

	Created by Shlomo Geva on 13/7/2023.
*/
#pragma once

#include <math.h>
#include <string.h>
#include <sys/stat.h>

#include <map>
#include <vector>

char *read_entire_file(const char *filename, uint64_t& fileSize);
char* index_kmers(const std::string& fastaFile, std::map<uint32_t, std::string>& referenceIDMap, std::vector<std::vector<uint32_t>>& kmersMap, uint32_t& MASK, uint64_t& genomeSize);
