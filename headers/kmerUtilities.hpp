/*
	KMERUTILITIES.HPP
	-----------------
	indexKISS

	Created by Shlomo Geva on 13/7/2023.
*/
#pragma once

#include <stdint.h>

#include <string>

uint64_t packKmer(const char *sequence);
std::string unpackKmer(uint64_t packed_sequence);
