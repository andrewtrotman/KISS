//
//  serialiseKmersMap.hpp
//  indexReference
//
//  Created by Shlomo Geva on 19/7/2023.
//

#ifndef serialiseKmersMap_hpp
#define serialiseKmersMap_hpp

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <cstring>

void serializeMap(const std::vector<std::vector<uint32_t>>& kmersMap, const std::string& innerMapFilename, const std::string& outerMapFilename);

void deserializeMap(const std::string& innerMapFilename, const std::string& outerMapFilename, std::vector<uint32_t>& innerMapBlob, std::vector<uint32_t>& outerMapBlob);

std::vector<uint32_t> getInnerVector(const std::vector<uint32_t>& innerMapBlob, const std::vector<uint32_t>& outerMapBlob, size_t index);

bool writeTextBlobToFile(const char* text, std::size_t length, const std::string& filename);

std::pair<char*, std::size_t> readTextBlobFromFile(const std::string& filename);

#endif /* serialiseKmersMap_hpp */

