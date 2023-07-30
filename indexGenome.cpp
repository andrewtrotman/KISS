/*
	INDEXGENOME.CPP
	---------------
	indexKISS

	Created by Shlomo Geva on 13/7/2023.
*/

#include <chrono>
#include <thread>
#include <iomanip>
#include <iostream>

#include "hash.hpp"
#include "indexGenome.hpp"
#include "packGenomeBlob.hpp"
#include "encode_kmer_2bit.h"

/*
	DISPLAYPROGRESS()
	-----------------
	// example of use:
	int totalIterations = 1000;
	int percentUpdateInterval = 10; // Update progress every 10%

	for (int i = 0; i < totalIterations; ++i)
		{
		displayProgress(i + 1, totalIterations, percentUpdateInterval);
		}
*/
void displayProgress(std::chrono::time_point<std::chrono::steady_clock> &start, uint64_t &lastDisplayedPercent, uint64_t current, uint64_t total, int desiredUpdateInterval)
	{
	uint64_t percent = (current * 100) / total;
	if (percent - lastDisplayedPercent >= desiredUpdateInterval)
		{
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();
		printf("Progress: %3llu%% (%lld milliseconds)\n", percent, duration);
		lastDisplayedPercent = percent;
		}
	}

/*
    READ_ENTIRE_FILE()
    ------------------
*/
char *read_entire_file(const char *filename, uint64_t &fileSize)
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
				contents = (char *)malloc(details.st_size + 1);
				if (fread(contents, details.st_size, 1, fp) != 1)
					{
					free(contents);
					contents = NULL;
					}
				else
					{
					fileSize = details.st_size;
					contents[details.st_size] = '\0';
					}
				}
			}
		fclose(fp);
		}
	return contents;
	}

/*
	LOAD_GENOME_FILE()
	------------------
*/
char *load_genome_file(const std::string &fastaFile, std::map<uint32_t, std::string> &referenceIDMap, uint64_t &genomeSize)
	{
	/*
		Load the genome file
	*/
	std::cout << std::endl << "Loading References: " << fastaFile << std::endl;
	uint64_t fileSize;
	char *genome = read_entire_file(fastaFile.c_str(), fileSize);
	if (genome == NULL)
		{
		std::cerr << "Failed to read " << fastaFile << " - Either missing, or larger than 4GB" << std::endl;
		exit(1);
		}
	std::cout << "Reference file size on disk " << fileSize << std::endl;

	/*
		Remove all non-base data from the file.
	*/
	genomeSize = packGenome(genome, fileSize, referenceIDMap);
	std::cout << "        Reference blob size " << genomeSize << std::endl;

	return genome;
	}

/*
	INDEX_KMERS_THREAD()
	--------------------
*/
void index_kmers_thread(char *genome, uint64_t genomeSize, std::vector<protected_vector<uint32_t>> &kmersMap, uint32_t MASK)
	{
//printf("%llu bytes from %p\n", genomeSize, genome);

    auto start = std::chrono::steady_clock::now();
	uint64_t lastDisplayedPercent = -10; // Initialize to a value that will trigger the first update
	/*
		Index by sliding a windows over the genome.  As the reverse complement is also needed, its done by
		keeping two "running windows" and shifting them then adding to the end.  That is, (pkmer << 2 | new_base)
		where new_base is the 2-bit encoding of the new base to add to the window and the encoding is 2 bits per base.
	*/
	uint64_t pkmer = encode_kmer_2bit::pack_32mer(genome);
	uint64_t remkp = encode_kmer_2bit::reverse_complement_32mer(pkmer);
	pkmer >>= 2;
	remkp <<= 2;
	char *encode_pos = genome + 31;
	for (uint32_t pos = 0; pos < genomeSize - 32; pos++)
		{
		uint64_t new_base = encode_kmer_2bit::pack_1mer(*encode_pos++);
		pkmer = (pkmer << 2) | new_base;
		remkp = (remkp >> 2) | (~new_base << 62);
		uint64_t cononical = pkmer ^ remkp;
		uint32_t kmerHash = murmurHash3(cononical) & MASK;
		kmersMap[kmerHash].push_back(pos);
		displayProgress(start, lastDisplayedPercent, pos, genomeSize, 10);
		}
	}

/*
	INDEX_KMERS()
	-------------
*/
void index_kmers(char *genome, uint64_t genomeSize, std::vector<protected_vector<uint32_t>> &kmersMap, uint32_t MASK)
	{
	size_t thread_count = std::thread::hardware_concurrency();
//	size_t thread_count = 2;
	uint64_t chunk_size = genomeSize / thread_count;
	uint64_t start = 0;

	/*
		Allocate the thread pool
	*/
	std::vector<std::thread> threads;

	/*
		Launch each thread
	*/
	std::cout << "Launching " << thread_count << " threads each with " << chunk_size << " pieces\n";
	for (size_t i = 0; i < thread_count - 1; i++)
		{
		threads.push_back(std::thread(index_kmers_thread, genome + start, chunk_size, std::ref(kmersMap), MASK));
		start += chunk_size;
		}
	index_kmers_thread(genome + start, genomeSize - start, kmersMap, MASK);

	/*
		Wait for each thread to terminate
	*/
	for (auto &thread : threads)
		thread.join();
	}
