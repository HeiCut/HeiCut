/******************************************************************************
 * mt_kahypar_io.h
 * *
 * Contains the I/O functions for the Mt-KaHyPar library.
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

#ifndef SMHM_MTKAHYPAR_IO_H
#define SMHM_MTKAHYPAR_IO_H

#include <cassert>
// Own headers
#include "lib/utils/output.h"
// Mt-KaHyPar headers
#include "mt-kahypar-library/libmtkahypar.h"
#include "mt-kahypar/partition/context_enum_classes.h"
#include "mt-kahypar/io/hypergraph_factory.h"

class MtKaHyParIO
{

private:
    // Method taken from mt-kahypar/lib/libmtkahypar.cpp (no modifications made)
    static mt_kahypar::PresetType to_preset_type(mt_kahypar_preset_type_t preset)
    {
        switch (preset)
        {
        case mt_kahypar_preset_type_t::DETERMINISTIC:
            return mt_kahypar::PresetType::deterministic;
        case mt_kahypar_preset_type_t::LARGE_K:
            return mt_kahypar::PresetType::large_k;
        case mt_kahypar_preset_type_t::DEFAULT:
            return mt_kahypar::PresetType::default_preset;
        case mt_kahypar_preset_type_t::QUALITY:
            return mt_kahypar::PresetType::quality;
        case mt_kahypar_preset_type_t::HIGHEST_QUALITY:
            return mt_kahypar::PresetType::highest_quality;
        }
        return mt_kahypar::PresetType::UNDEFINED;
    }

public:
    template <typename Config>
    static mt_kahypar_hypergraph_t read_hypergraph_from_file(Config &config)
    {
        /// Get the number of available threads
        const size_t numAvailableThreads = std::thread::hardware_concurrency();

        // Set the number of threads to use
        config.numThreads  = (config.numThreads > 0 && config.numThreads <= numAvailableThreads) ? config.numThreads : numAvailableThreads;

        // Initialize thread pool
        mt_kahypar_initialize_thread_pool(config.numThreads, true /* activate interleaved NUMA allocation policy */);

        std::cout << "num_threads \t\t\t\t" << config.numThreads << std::endl;

        // NB: We do not remove single pin hyperedges during IO, since they are removed in the pruner
        bool remove_single_pin_hes = false;

        std::cout << "file_format \t\t\t\t" << config.hypergraphFileFormat << std::endl;
        std::cout << "preset_type \t\t\t\t" << config.presetType << std::endl;

        mt_kahypar_hypergraph_t hypergraphWrapper = mt_kahypar_read_hypergraph_from_file(config.hypergraphFileName, config.presetType, config.hypergraphFileFormat, remove_single_pin_hes);

        // Make sure that the hypergraph is compatible with the preset type
        assert(mt_kahypar_check_compatibility(hypergraphWrapper, config.presetType));

        return hypergraphWrapper;
    };

    // Method taken from mt-kahypar/lib/libmtkahypar.cpp
    // Changes:
    //  - Additional possibility to specify whether to remove single pin hyperedges
    //  - The instance type is always hypergraph, even if the input format is METIS
    static mt_kahypar_hypergraph_t mt_kahypar_read_hypergraph_from_file(const char *file_name,
                                                                        const mt_kahypar_preset_type_t preset,
                                                                        const mt_kahypar_file_format_type_t file_format,
                                                                        const bool remove_single_pin_hes)
    {
        const mt_kahypar::PresetType config = to_preset_type(preset);
        const mt_kahypar::InstanceType instance = mt_kahypar::InstanceType::hypergraph;
        const mt_kahypar::FileFormat format = file_format == mt_kahypar_file_format_type_t::HMETIS ? mt_kahypar::FileFormat::hMetis : mt_kahypar::FileFormat::Metis;
        const bool stable_construction = preset == mt_kahypar_preset_type_t::DETERMINISTIC ? true : false;
        try
        {
            return mt_kahypar::io::readInputFile(file_name, config, instance, format, stable_construction, remove_single_pin_hes);
        }
        catch (std::exception &ex)
        {
            LOG << ex.what();
        }
        return mt_kahypar_hypergraph_t{nullptr, NULLPTR_HYPERGRAPH};
    }
};

#endif // end of SMHM_MTKAHYPAR_IO_H