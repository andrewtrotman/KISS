/*
	PACKGENOMEBLOB.HPP
	------------------
	indexReference

	Created by Shlomo Geva on 22/7/2023.
*/
#pragma once

#include <stdint.h>

#include <map>

size_t packGenome(char *genome, uint64_t genome_size, std::map<std::uint32_t, std::string> &referenceIDMap);
