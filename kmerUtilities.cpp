//
//  kmerUtilities.cpp
//  indexKISS
//
//  Created by Shlomo Geva on 13/7/2023.
//

#include <cstdint>
#include <string>

extern int KMERSIZE;
extern uint64_t kmer_encoding_table[];
/*
    PACK KMER
    -----------
*/
uint64_t packKmer(const char *sequence)
    {
    uint64_t packed = 0;

    for (int pos = 0; pos < KMERSIZE; pos++)
        packed = (packed << 2) | kmer_encoding_table[(size_t)sequence[pos]];

    return packed;
    }

/*
    UNPACK KMER
    -------------
*/
std::string unpackKmer(uint64_t packed_sequence)
    {
    std::string sequence;
    for (int32_t pos = KMERSIZE-1; pos >= 0; pos--)
        {
        switch ((packed_sequence >> (pos * 2)) & 3)
            {
            case 0:
                sequence += 'A';
                break;
            case 1:
                sequence += 'C';
                break;
            case 2:
                sequence += 'G';
                break;
            case 3:
                sequence += 'T';
                break;
            }
        }
    return sequence;
    }


