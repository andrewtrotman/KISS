//
//  indexGenome.hpp
//  indexKISS
//
//  Created by Shlomo Geva on 13/7/2023.
//

#ifndef indexGenome_hpp
#define indexGenome_hpp

#include <iostream>
#include <iomanip>
#include <sys/stat.h>
#include <map>
#include <vector>
#include <cmath>
#include <limits>
#include <cstring>
#include <algorithm>

char *read_entire_file(const char *filename, uint64_t& fileSize);
char* index_kmers(const std::string& fastaFile, int KMERSIZE,
                  std::map<uint32_t, std::string>& referenceIDMap,
                  std::vector<std::vector<uint32_t>>& kmersMap,
                  uint32_t& MASK,
                  uint64_t& genomeSize);

#endif /* indexGenome_hpp */

