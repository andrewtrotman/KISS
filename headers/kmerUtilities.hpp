//
//  kmerUtilities.hpp
//  indexKISS
//
//  Created by Shlomo Geva on 13/7/2023.
//

#ifndef kmerUtilities_hpp
#define kmerUtilities_hpp

#include <stdio.h>

uint64_t packKmer(const char *sequence);
std::string unpackKmer(uint64_t packed_sequence);

#endif /* kmerUtilities_hpp */
