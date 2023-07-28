//
//  packGenomeBlob.hpp
//  indexReference
//
//  Created by Shlomo Geva on 22/7/2023.
//

#ifndef packGenomeBlob_hpp
#define packGenomeBlob_hpp

#include <stdio.h>
#include <cstdint>
#include <iostream>
#include <map>

std::size_t packGenome(char* genome, std::map<std::uint32_t, std::string>& referenceIDMap);

#endif /* packGenomeBlob_hpp */

