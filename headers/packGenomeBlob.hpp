/*
	PACKGENOMEBLOB.HPP
	------------------
	indexReference

	Created by Shlomo Geva on 22/7/2023.
*/
#pragma once

#include <stdint.h>

#include <map>

std::size_t packGenome(char* genome, std::map<std::uint32_t, std::string>& referenceIDMap);
