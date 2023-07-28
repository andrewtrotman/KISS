/*
	MAIN.CPP
	--------
	indexReference

	Created by Shlomo Geva on 16/7/2023.
*/

#include <map>
#include <chrono>
#include <sstream>
#include <fstream>
#include <iostream>

#include "indexGenome.hpp"
#include "serialiseKmersMap.hpp"


using std::cout;
using std::endl;
using std::cerr;
std::vector<std::string> sampleSequences;
std::vector<std::string> fileNames;

// some global default values (overide with cmd line arguments)
unsigned int KMERSIZE; // seed kmer size
std::string REFERENCE; // file name for reference file to match against
std::string OUTPUT_DIR; // directory name for index
uint64_t kmer_encoding_table[256]; // for packing a kmer into uint64

/*
    WRITE MAP TO FILE
 */
void writeMapToFile(const std::string& filename, const std::map<uint32_t, std::string>& referenceIDMap) {
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error opening the file: " << filename << std::endl;
        return;
    }

    for (const auto& entry : referenceIDMap) {
        outFile << entry.first << " " << entry.second; // << std::endl;
    }

    outFile.close();
}


/*
    EXTRACT REF NAME
 */
std::string extractRefName(const std::string& filePath) {
    // Extract the bare file name from the file path
    std::string fileName = filePath.substr(filePath.find_last_of('/') + 1);

    // Remove the file extension
    std::string baseName = fileName.substr(0, fileName.find_last_of('.'));

    return baseName;
}
/*
    GET BASE NAME
 */
std::string getBaseName(const std::string& filePath) {
    std::stringstream ss(filePath);
    std::string baseName;
    std::getline(ss, baseName, '.');
    return baseName;
}

/*
    GET REFERENCE
*/
void getReference(std::string inputFile,
                  std::map<uint32_t, std::string> referenceIDMap,
                  std::vector<std::vector<uint32_t>>& kmersMap,
                  char*& genome,
                  uint32_t& MASK) {
    // load the map
    auto start = std::chrono::steady_clock::now();
    uint64_t genomeSize;
    genome = index_kmers(REFERENCE,referenceIDMap,kmersMap,MASK,genomeSize);
    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    int minutes = (int) duration.count() / (1000 * 60);
    float seconds = ((duration.count() - minutes * 1000 * (float) 60))/1000;
    std::cout << "Building time: " << minutes << " min " << seconds << " sec" << std::endl;

    // Serialize the map
    start = std::chrono::steady_clock::now();
    std::string outerMapFilename = getBaseName(inputFile)+"_"+std::to_string(KMERSIZE)+"_OuterBlob.idx";
    std::string innerMapFilename = getBaseName(inputFile)+"_"+std::to_string(KMERSIZE)+"_InnerBlob.idx";
    std::string genomeFilename = getBaseName(inputFile)+"_genome.idx";
    std::string refIDFilename = getBaseName(inputFile)+"_refID.idx";

    std::cout << "Serialising genome to " << genomeFilename << " and " << innerMapFilename << std::endl;
    writeTextBlobToFile(genome, genomeSize, genomeFilename);
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    minutes = (int) duration.count() / (1000 * 60);
    seconds = ((duration.count() - minutes * 1000 * (float) 60))/1000;
    std::cout << "Serialising genome time: " << minutes << " min " << seconds << " sec" << std::endl;

    std::cout << "Serialising map to " << outerMapFilename << " and " << innerMapFilename << std::endl;
    serializeMap(kmersMap, innerMapFilename, outerMapFilename);
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    minutes = (int) duration.count() / (1000 * 60);
    seconds = ((duration.count() - minutes * 1000 * (float) 60))/1000;
    std::cout << "Serialising Maps time: " << minutes << " min " << seconds << " sec" << std::endl;
    
    std::cout << "Serialising ReferenceIDMap" << std::endl;
    writeMapToFile(refIDFilename, referenceIDMap);
        
    // DeSerialize the genome
    start = std::chrono::steady_clock::now();
    // Read the text blob from the file and directly assign to genome and textLength
    std::tie(genome, genomeSize) = readTextBlobFromFile(genomeFilename);
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    minutes = (int) duration.count() / (1000 * 60);
    seconds = ((duration.count() - minutes * 1000 * (float) 60))/1000;
    std::cout << "DeSerialising genome time: " << minutes << " min " << seconds << " sec" << std::endl;

    // DeSerialize the map
    start = std::chrono::steady_clock::now();
    std::vector<uint32_t> innerMapBlob;
    std::vector<uint32_t> outerMapBlob;
    deserializeMap(innerMapFilename, outerMapFilename, innerMapBlob, outerMapBlob);
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    minutes = (int) duration.count() / (1000 * 60);
    seconds = ((duration.count() - minutes * 1000 * (float) 60))/1000;
    std::cout << "DeSerialising Maps time: " << minutes << " min " << seconds << " sec" << std::endl;

/*  SANITY TEST CODE, ignore
    // test the index
    // test the index
    // Use the helper function to access the map elements
    size_t index;
    for (index=1000000; index<10000000; index++) {
        if (!(kmersMap[index].empty())) {
            break;
        }
    }
    // Use the helper function to access the map elements
    std::vector<uint32_t> innerVector = getInnerVector(innerMapBlob, outerMapBlob, index);
    // Output the inner vector retrieved from the blob
    for (const auto& value : innerVector) {
        std::cout << value << " ";
    }
    std::cout << std::endl;
    for (const auto& value : kmersMap[index]) {
        std::cout << value << " ";
    }
    std::cout << std::endl;
 */
}

/*
    INITIALISE
 */

void intialise(int argc, char* argv[]) {
    std::string arg;
    if (argc>1) {
        arg = argv[1];
    }
    if ((argc==1) || arg=="-help") {
        cout << "Usage: ./indexReference -argName1 -argValue1 -argName2 -argValue2 ... " << endl;
        cout << "ArgMames are as follows:" << endl;
        cout << "-kmerSize is the number of bases in the seed kmer (recommended 16 to 32)" << endl;
        cout << "-refernce is the name of the reference file" << endl;
        cout << "-out_dir is the name of directory where the index is stored" << endl;
        cout << "example ./indexReference -kmer_size 30 -reference CutibacteriumGenome.fasta -out_dir temp" <<endl;
        exit(0);
    }
    // some default values (overide with cmd line arguments, see below)
    KMERSIZE = 30;
    REFERENCE = "";/* "/Users/geva/Crispr/AMR106.fasta"; */
    OUTPUT_DIR = "/Users/geva/temp";

    for (int i = 1; i < argc; ++i) {
        arg = argv[i];
        if (i + 1 >= argc) {
            std::cout << "Error: Missing value for " << arg << " option." << std::endl;
            continue;
        }

        std::string value = argv[i + 1];
        ++i;  // Skip the next argument

        // Process the command line option
        if (arg == "-kmerSize") {
            KMERSIZE = std::stoi(value);
        } else if (arg == "-reference") {
            REFERENCE = value;
        } else if (arg == "-out_dir") {
            OUTPUT_DIR = value;
        } else {
            cerr << "Error: Unknown option: " << arg << std::endl;
        }
    }
    // set up gobal encodingTable for packKmer
    kmer_encoding_table[(size_t)'A'] = 0;
    kmer_encoding_table[(size_t)'C'] = 1;
    kmer_encoding_table[(size_t)'G'] = 2;
    kmer_encoding_table[(size_t)'T'] = 3;
    
    cout << "indexReference run parameters" << endl;
    cout << "kmerSize: " << KMERSIZE << endl;
    cout << "reference: " << REFERENCE << endl;
    cout << "output directory: " << OUTPUT_DIR << endl;
}

/*
     MAIN
 */
int main(int argc, char *argv[]) {

    std::map<uint32_t, std::string> referenceIDMap;
    std::vector<std::vector<uint32_t>> kmersMap;
    char* genome=NULL;
    uint32_t MASK;

    // Start the timer
    auto startAll = std::chrono::steady_clock::now();

    // set up KISS parameters
    intialise(argc, argv);
    getReference(REFERENCE,referenceIDMap,kmersMap,genome,MASK); // load the reference collection index
    
    // report overall program duration
    auto endAll = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endAll - startAll);
    int minutes = (int) duration.count() / (1000 * 60);
    float seconds = ((duration.count() - minutes * 1000 * (float) 60))/1000;
    std::cout << endl << "Total time: " << minutes << " min " << seconds << " sec" << std::endl;
    return 0;
}


