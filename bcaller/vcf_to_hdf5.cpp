#include <iostream>
#include <Variant.h>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Easy.hpp>

struct hash_pair {
    template <class T1, class T2>
    size_t operator()(const pair<T1, T2>& p) const
    {
        auto hash1 = hash<T1>{}(p.first);
        auto hash2 = hash<T2>{}(p.second);
        return hash1 ^ hash2;
    }
};

template <typename T>
void writeData(HighFive::File& outFile, const std::string& curSequence, const std::vector<std::vector<T>>& data) {
    if(!outFile.exist("/" + curSequence)){
        std::vector<size_t> dims(2);
        dims[0] = data.size();
        dims[1] = data[0].size();
        const std::string dsetName = "/" + curSequence;
        auto dspace = HighFive::DataSpace(dims);
        outFile.createDataSet<T>(dsetName, dspace);
    }
    auto dset = outFile.getDataSet("/" + curSequence);
    dset.write(data);
}

int main(int argc, char** argv) {

    if(argc != 3) {
        std::cout << "Args: 1. VCF File 2. Output HDF file" << std::endl;
        return 1;
    }

    vcflib::VariantCallFile variantFile;
    std::vector<std::vector<size_t>> locAndData;
    locAndData.reserve(10000);
    std::unordered_map<std::pair<std::string, std::string>, size_t, hash_pair> mutFrequencies;
    HighFive::File outFile(argv[2], H5Easy::File::ReadWrite | H5Easy::File::Create);
    std::string curSequence;

    std::string vcfFilename = argv[1];
    variantFile.open(vcfFilename);

    if (!variantFile.is_open()) {
        return 1;
    }

    vcflib::Variant var(variantFile);
    for(size_t i = 0; variantFile.getNextVariant(var); ++i) {
        if(i % 100 == 0){
            std::cout << "Read " << i << " records" << std::endl << std::flush;
        }
        if(curSequence != var.sequenceName) {
            if(!curSequence.empty()) {

                writeData(outFile, curSequence, locAndData);
                locAndData.clear();
            }
            curSequence = var.sequenceName;
        }
        auto pos = static_cast<size_t>(var.position);
        auto ac = static_cast<size_t>(var.getInfoValueFloat("AC", 0));
        auto an = static_cast<size_t>(var.getInfoValueFloat("AN", 0));
        locAndData.push_back(std::vector<size_t>{pos, ac, an});
        for(auto& alt : var.alt) {
            if(var.ref.length() == 1 && alt.length() == 1) {
                auto key = std::make_pair(var.ref, alt);
                //size_t defaults to zero
                auto &val = mutFrequencies[key];
                mutFrequencies[key] = val + 1;
            }
        }

    }
    if(!curSequence.empty()) {
        writeData(outFile, curSequence, locAndData);
    }

    for(auto& pair : mutFrequencies) {
        auto& refAndAlt = pair.first;
        auto& count = pair.second;
        auto path = "/freqs/" + refAndAlt.first + "/" + refAndAlt.second;
        auto loadedCounts = 0;
        if(outFile.exist(path)) {
            loadedCounts = H5Easy::load<size_t>(outFile, path);
        }
        H5Easy::dump(outFile, path, loadedCounts + count, H5Easy::DumpMode::Overwrite);
    }


    return 0;

}
