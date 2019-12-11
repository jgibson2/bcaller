//
// Created by john on 11/13/19.
//
#include <iostream>
#include <string>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Easy.hpp>

#define BASAL_MUTATION_RATE 1.1E-8

int main(int argc, char** argv) {
    if(argc != 5) {
        std::cout << "Args: 1. Input HDF file 2. Output HDF File 3. Contig name 4. Window size" << std::endl;
        return 1;
    }

    H5Easy::File inFile(argv[1], H5Easy::File::ReadOnly);
    H5Easy::File outFile(argv[2], H5Easy::File::ReadWrite | H5Easy::File::Create);
    std::string contig(argv[3]);
    auto windowSize = std::stoul(argv[4]);
    size_t currentLoc = 0;

    //this should be sorted already
    auto locAcAnData = H5Easy::load<std::vector<std::vector<size_t>>>(inFile, "/" + contig);
    std::vector<unsigned int> numAltInWindow;
    std::vector<unsigned int> totalInWindow;
    std::vector<std::vector<double>> windowPositionAndVariantProbs;
    double basalMutationRate = BASAL_MUTATION_RATE;

    for(auto& vec : locAcAnData) {
        auto& loc = vec[0];
        auto& ac = vec[1];
        auto& an = vec[2];
        if( currentLoc % (windowSize * 1000) == 0) {
            std::cout << "Current location: " << loc << std::endl << std::flush;
        }
        while(loc - currentLoc >= windowSize){
            auto acSum = std::accumulate(numAltInWindow.cbegin(), numAltInWindow.cend(), 0, std::plus<>());
            auto anMean = totalInWindow.empty() ? 0.0 : static_cast<double>(std::accumulate(totalInWindow.cbegin(), totalInWindow.cend(), 0, std::plus<>())) / totalInWindow.size();
            windowPositionAndVariantProbs.push_back(
                    std::vector<double>{
                        static_cast<double>(currentLoc),
                        std::max(acSum == 0 ? 0.0 :static_cast<double>(acSum) / (windowSize * anMean), basalMutationRate)
                    });

            numAltInWindow.clear();
            totalInWindow.clear();
            currentLoc += windowSize;
        }
        numAltInWindow.push_back(ac);
        totalInWindow.push_back(an);
    }
    auto acSum = std::accumulate(numAltInWindow.cbegin(), numAltInWindow.cend(), 0, std::plus<>());
    auto anMean = static_cast<double>(std::accumulate(totalInWindow.cbegin(), totalInWindow.cend(), 0, std::plus<>())) / totalInWindow.size();
    windowPositionAndVariantProbs.push_back(
            std::vector<double>{
                    static_cast<double>(currentLoc),
                    std::max(acSum == 0 ? 0.0 :static_cast<double>(acSum) / (windowSize * anMean), basalMutationRate)
            });



    H5Easy::dump(outFile, "/" + contig, windowPositionAndVariantProbs, H5Easy::DumpMode::Overwrite);
    if(outFile.hasAttribute("/" + contig + "/window_size")) {
        outFile.getAttribute("/" + contig + "/window_size").write(windowSize);
    } else {
        outFile.createAttribute("/" + contig + "/window_size", windowSize);
    }
    for(auto& a : std::vector<char>{'A', 'T', 'C', 'G'}) {
        for(auto& b : std::vector<char>{'A', 'T', 'C', 'G'}) {
            if(a != b) {
                auto path = std::string("/freqs/") + a + '/' + b;
                auto count = H5Easy::load<uint64_t>(inFile, path);
                auto loadedCounts = 0;
                if(outFile.exist(path)) {
                    loadedCounts = H5Easy::load<size_t>(outFile, path);
                }
                H5Easy::dump(outFile, path, loadedCounts + count, H5Easy::DumpMode::Overwrite);
            }
        }
    }

    return 0;
}