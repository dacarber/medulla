/**
 * @file trees.h
 * @brief Header file for the trees namespace.
 * @details This file contains the header for the trees namespace. The trees
 * namespace contains functions that read and interface with the TTrees
 * produced by the CAFAna analysis framework. Different "copying" actions can
 * be performed on the TTrees, such as a simple copy or adding systematics to
 * the output file based on the selected signal candidates and the configured
 * systematics.
 * @author mueller@fnal.gov
 */
#ifndef TREES_H
#define TREES_H
#include <iostream>
#include <cstdint>
#include <functional>
#include <string>
#include <tuple>
#include <utility>

#include "detsys.h"
#include "configuration.h"

#include "TFile.h"

/**
 * @namespace sys::trees
 * @brief Namespace for functions that read and interface with the TTrees
 * produced by the CAFAna analysis framework.
 * @details This namespace contains functions that read and interface with the
 * TTrees produced by the CAFAna analysis framework. Different "copying" 
 * actions can be performed on the TTrees, such as a simple copy or the
 * addition of systematics to the output file based on the selected signal
 * candidates and the configured systematics.
 */
namespace sys::trees
{
    /**
     * @brief Type definitions for systematic indexing (variable name and
     * index).
     */
    typedef std::pair<std::string, int64_t> syst_t;
    typedef std::tuple<uint64_t, uint64_t, uint64_t, uint64_t, double> index_t;

    /**
     * @brief Hash functor for the index_t tuple key.
     * @details Combines the run, subrun, event, neutrino-id, and energy fields
     * using a 64-bit boost-style mixing constant. Used so that the selected
     * signal candidate lookup can be performed in O(1) via an unordered_map /
     * unordered_set instead of an O(log N) red-black tree.
     */
    struct index_hash
    {
        std::size_t operator()(const index_t & k) const noexcept
        {
            auto mix = [](std::size_t h, std::size_t v) noexcept
            {
                return h ^ (v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
            };
            std::size_t h = std::hash<uint64_t>{}(std::get<0>(k));
            h = mix(h, std::hash<uint64_t>{}(std::get<1>(k)));
            h = mix(h, std::hash<uint64_t>{}(std::get<2>(k)));
            h = mix(h, std::hash<uint64_t>{}(std::get<3>(k)));
            h = mix(h, std::hash<double>{}(std::get<4>(k)));
            return h;
        }
    };

    /**
     * @brief Hash functor for the syst_t pair key (variable name, index).
     * @details Used so that the per-systematic result histogram lookup can be
     * performed in O(1) instead of paying an O(log N) string comparison on
     * every universe.
     */
    struct syst_hash
    {
        std::size_t operator()(const syst_t & k) const noexcept
        {
            std::size_t h = std::hash<std::string>{}(k.first);
            return h ^ (std::hash<int64_t>{}(k.second) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
        }
    };

    /**
     * @brief Copy the input TTree to the output TTree.
     * @details This function copies the input TTree to the output TTree. The
     * function loops over the input TTree and copies the values of the branches
     * to the output TTree. The output TTree is created with the same branches
     * as the input TTree.
     * @param table The table that contains the configuration for the tree.
     * @param output The output TFile.
     * @param input The input TFile.
     * @return void
     */
    void copy_tree(cfg::ConfigurationTable & table, TFile * output, TFile * input);

    /**
     * @brief Add reweightable systematics to the output TTree.
     * @details This function adds reweightable systematics to the output
     * TTree. The function loops over the input TTree to build a map for the
     * selected signal candidates to their index in the input TTree. The
     * function then loops over the neutrinos in the CAF input files and
     * populates the output TTree with the selected signal candidates and the
     * universe weights for matched neutrinos.
     * @param table The table that contains the configuration for the tree.
     * @param output The output TFile.
     * @param input The input TFile.
     * @return void
     */
    void copy_with_weight_systematics(cfg::ConfigurationTable & config, cfg::ConfigurationTable & table, TFile * output, TFile * input, sys::detsys::DetsysCalculator & calc);
}
#endif