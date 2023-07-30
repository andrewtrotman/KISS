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

#include "protected_vector.hpp"

char *read_entire_file(const char *filename, uint64_t& fileSize);
char *load_genome_file(const std::string &fastaFile, std::map<uint32_t, std::string> &referenceIDMap, uint64_t &genomeSize);
void index_kmers(char *genome, uint64_t genomeSize, std::vector<protected_vector<uint32_t>> &kmersMap, uint32_t MASK);

