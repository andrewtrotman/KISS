//
//  indexGenome.cpp
//  indexKISS
//
//  Created by Shlomo Geva on 13/7/2023.
//

#include <iostream>
#include <iomanip>

#include "hash.hpp"
#include "indexGenome.hpp"
#include "packGenomeBlob.hpp"
#include "encode_kmer_2bit.h"

void displayProgress(uint64_t current, uint64_t total, int desiredUpdateInterval) {
    uint64_t percent = (current * 100) / total;
    static uint64_t lastDisplayedPercent = -desiredUpdateInterval; // Initialize to a value that will trigger the first update
    if (percent - lastDisplayedPercent >= desiredUpdateInterval) {
        std::cout << "Progress: " << std::setw(3) << percent << "%" << std::endl;
        std::cout.flush();
        lastDisplayedPercent = percent;
    }
    /*
     // example of use:
     int totalIterations = 1000;
     int percentUpdateInterval = 10; // Update progress every 10%
     
     for (int i = 0; i < totalIterations; ++i) {
     // Your loop processing here
     
     displayProgress(i + 1, totalIterations, percentUpdateInterval);
     }
     */
}

/*
     CONVERT TO UPPERCASE
 */
void convertToUppercase(char* text, std::size_t size) {
    for (std::size_t i = 0; i < size; ++i) {
        text[i] = std::toupper(static_cast<unsigned char>(text[i]));
    }
}

/*
    READ_ENTIRE_FILE
 */
char *read_entire_file(const char *filename, uint64_t& fileSize)
  {
  FILE *fp;
  struct stat details;
  char *contents = NULL;
  fileSize = 0;
  if ((fp = fopen(filename, "rb")) != NULL)
  {
      if (fstat(fileno(fp), &details) == 0)
      {
          if (details.st_size != 0 || details.st_size > UINT32_MAX)
          {
              contents = (char *)malloc(details.st_size);
              if (fread(contents, details.st_size, 1, fp) != 1)
              {
                  free(contents);
                  contents = NULL;
              }
              else
              {
                  fileSize = details.st_size;
              }
          }
      }
      fclose(fp);
   }
    convertToUppercase(contents, fileSize);
   return contents;
}

/*
 INDEX_KMERS
 */
char* index_kmers(const std::string& fastaFile, int KMERSIZE,
                  std::map<uint32_t, std::string>& referenceIDMap,
                  std::vector<std::vector<uint32_t>>& kmersMap,
                  uint32_t& MASK,
                  uint64_t& genomeSize)
{
//std::vector<std::vector<uint32_t>> index_kmers(const std::string& fastaFile, int KMERSIZE, char *genome) {
    // Create a kmer index for the genome.
    // Index information is limited to just a unit32, holding the position of the kmer in the genome fasta file.
    // Size is upper-limited by the size of the genome in bytes (assuming one byte per base).
    // For instance, 2^20 map size will be created for a file of size 1MB
    // - there can be at most 1M unique kmers in such a reference.
    // However, if the filesize is 3.5GB (e.g. human genome) then the index will have room for 2^32-1 kmers.
    // Each kmer is in the current implementation is restricted to the range of 16mer to 32mer.
    /*
    // A kmer position in the map is found by packing the kmer onto uint64, and xoring it with its packed uint64 reverse complement.  This gives a canonical representation to kmer.
    // Shorter than 32mers will map to fewer bits, using the least signficant bits.
    // Hashing the uint64 onto a 32bit integer follows, using murmurHash3().
    */
    // A kmer entry consists of the positions in the genome where the kmer is found.
    // The position is also limited by the filesize to 2^32, so it fits on a uint32.
    //
    // Note that the kmer size is limited to 32mers because we pack a base on 2 bits of a uint64.
    // This restriction can be lifted by using kmer random projection onto 64 bits,
    //     and then the kmer length can be arbitrary.

    // Read reference file into memory
    std::cout << std::endl << "Loading References: " << fastaFile << std::endl;
    std::string line;
    std::string sequence; // Variable to collect the DNA sequence
    uint64_t fileSize;
    char *genome = read_entire_file(fastaFile.c_str(),fileSize);
    if (genome==NULL) {
        std::cerr << "Failed to read " << fastaFile << " - Either missing, or larger than 4GB" << std::endl;
        exit(1);
    }
    std::cout << "Reference file size on disk " << fileSize << std::endl;
//    genomeSize = packGenome(genome, fileSize, referenceIDMap);
    genomeSize = packGenome(genome, referenceIDMap);
    std::cout << "        Reference blob size " << genomeSize << std::endl;
    fileSize = genomeSize;
   // Calculate the number of elements to reserve in kmersIndex based on genome size
    int numBitsToKeep = ::ceil(::log2(fileSize));
    if (numBitsToKeep == 32) {
        MASK = UINT32_MAX; // Set all bits to 1
    } else {
        MASK = (1 << numBitsToKeep) - 1;
    }
    std::cout << "Keeping " << numBitsToKeep << " bits in kmerHash" << std::endl;
    kmersMap.resize(pow(2,numBitsToKeep));
    // indexing proceeds by while keeping track of character position in the genome file, skipping header lines
    // in the fasta file, and indexing kmers of the size specified by KMERSIZE.
    // A kmer is first hashed into a 32 bits integer - a kmerIndex into the kmersMap vector.
    // Then the position of the kmer is entered into the set in that position.
    // example:  kmerMap[kmerIndex].insert(kmerPos);
    // All positions are relative to the begning of the genome file in memory (char*).
    // In this manner, the index can be trivially serialised as a binary dump to file.
    // When the genome is reloaded into memory in a different location, the index is still valid.

    // just for counting the number of unique kmers
//    std::set<uint64_t> uniqueKmers; // The map will store unique integers as keys and their count as values.
//    std::set<uint32_t> uniqueKmerHashes; // The map will store unique integers as keys and their count as values.

    uint32_t kmerHash;
    uint64_t kmerCount=0;
    uint64_t pkmer;
    std::string kMer;
    std::string genomeStr(genome);
    for (uint32_t pos=0; pos < genomeStr.size(); pos++)
    {
        // Generate and store packed kmers and their occurrences
        kMer = (genomeStr.substr(pos, KMERSIZE));
        // Add the occurrence to the kMerMap;
        pkmer = encode_kmer_2bit::pack_32mer(kMer.c_str());
        pkmer = pkmer ^ encode_kmer_2bit::reverse_complement_32mer(pkmer);// canonical kmer representation
        if (MASK<UINT32_MAX)
            kmerHash = murmurHash3(pkmer) & MASK;
        else
            kmerHash = murmurHash3(pkmer);
        kmersMap[kmerHash].push_back(pos);
        ++kmerCount;
        displayProgress((u_int64_t) pos, genomeStr.size(), 10);
    }

    uint64_t kmersInMap=0;
    int count=0;
    for (int i=0; i<kmersMap.size(); i++) {
        if (kmersMap[i].empty()) {
            kmersInMap++;
            continue;
        }
        count++;
   }
    std::cout  << "Map size " << kmersMap.size() << ", kmersCount " << kmerCount
               << ", kmers in Map " << kmersInMap << std::endl;
    return genome;
}


