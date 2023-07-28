/*
	HASH.HPP
	--------
	indexKISS

	Created by Shlomo Geva on 13/7/2023.
*/
#pragma once

#include <stdio.h>
#include <stdint.h>

uint32_t murmurHash3(uint64_t key);
uint32_t xorHash(uint64_t packedKmer);
