// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "active_region_generator.hpp"

#include <iterator>

#include "utils/repeat_finder.hpp"
#include "utils/mappable_algorithms.hpp"
#include "utils/append.hpp"

namespace octopus { namespace coretools {

ActiveRegionGenerator::ActiveRegionGenerator(const ReferenceGenome& reference, Options options)
: reference_ {reference}
, options_ {options}
, assembler_name_ {"LocalReassembler"}
, cigar_scanner_name_ {"CigarScanner"}
, using_assembler_ {false}
, assembler_active_region_generator_ {}
, max_read_length_ {}
{}

void ActiveRegionGenerator::add_generator(const std::string& name)
{
    if (is_assembler(name)) {
        using_assembler_ = true;
        if (!options_.assemble_all) {
            assembler_active_region_generator_ = AssemblerActiveRegionGenerator {reference_};
        }
    }
}

void ActiveRegionGenerator::add_read(const SampleName& sample, const AlignedRead& read)
{
    max_read_length_ = std::max(max_read_length_, sequence_size(read));
    if (assembler_active_region_generator_) assembler_active_region_generator_->add(sample, read);
}

auto merge(std::vector<GenomicRegion> lhs, std::vector<GenomicRegion> rhs)
{
    auto itr = utils::append(std::move(rhs), lhs);
    std::inplace_merge(std::begin(lhs), itr, std::end(lhs));
    return extract_covered_regions(lhs);
}

template <typename Range, typename T>
auto append(const Range& range, std::vector<T>& result)
{
    return result.insert(std::cend(result), std::cbegin(range), std::cend(range));
};

auto find_minisatellites(const std::vector<TandemRepeat>& repeats, const GenomicRegion& region,
                         const std::size_t max_read_length)
{
    InexactRepeatDefinition minisatellite_def {};
    minisatellite_def.min_exact_repeat_seed_length = 2 * max_read_length / 3;
    minisatellite_def.min_exact_repeat_seed_periods = 3;
    minisatellite_def.max_seed_join_distance = max_read_length / 3;
    minisatellite_def.min_joined_repeat_length = 2 * max_read_length;
    return join(find_repeat_regions(repeats, region, minisatellite_def), max_read_length / 2);
}

auto find_compound_microsatellites(const std::vector<TandemRepeat>& repeats, const GenomicRegion& region,
                                   const std::size_t max_read_length)
{
    InexactRepeatDefinition compound_microsatellite_def {};
    compound_microsatellite_def.max_exact_repeat_seed_period = 6;
    compound_microsatellite_def.min_exact_repeat_seed_length = 4;
    compound_microsatellite_def.min_exact_repeat_seed_periods = 4;
    compound_microsatellite_def.min_exact_seeds = 2;
    compound_microsatellite_def.max_seed_join_distance = 1;
    compound_microsatellite_def.min_joined_repeat_length = max_read_length / 4;
    return join(find_repeat_regions(repeats, region, compound_microsatellite_def), max_read_length / 2);
}

auto find_adjacent_homopolymer_regions(const std::vector<TandemRepeat>& homopolymers)
{
    std::vector<GenomicRegion> blocks {};
    blocks.reserve(homopolymers.size());
    for (auto lhs_itr = std::cbegin(homopolymers), rhs_itr = std::next(lhs_itr); rhs_itr != std::cend(homopolymers); ++lhs_itr, ++rhs_itr) {
        if (are_adjacent(*lhs_itr, *rhs_itr)) blocks.push_back(encompassing_region(*lhs_itr, *rhs_itr));
    }
    return extract_covered_regions(blocks);
}

auto find_adjacent_homopolymer_regions(const ReferenceGenome& reference, const GenomicRegion& region,
                                       const unsigned min_repeat_len)
{
    auto homopolymers = find_exact_tandem_repeats(reference, region, 1);
    homopolymers.erase(std::remove_if(std::begin(homopolymers), std::end(homopolymers),
                                      [=] (const auto& homopolymer) { return region_size(homopolymer) >= min_repeat_len; }),
                                      std::end(homopolymers));
    return find_adjacent_homopolymer_regions(homopolymers);
}

auto merged_covered_regions(std::vector<GenomicRegion> lhs, std::vector<GenomicRegion> rhs)
{
    auto itr = utils::append(std::move(rhs), lhs);
    std::inplace_merge(std::begin(lhs), itr, std::end(lhs));
    return extract_covered_regions(lhs);
}

std::vector<GenomicRegion> ActiveRegionGenerator::generate(const GenomicRegion& region, const std::string& generator) const
{
    if (is_assembler(generator) && assembler_active_region_generator_) {
        auto default_assembler_regions = assembler_active_region_generator_->generate(region);
        auto homopolymer_regions = find_adjacent_homopolymer_regions(reference_, region, 4);
        return merged_covered_regions(std::move(default_assembler_regions), std::move(homopolymer_regions));
    } else {
        return {region};
    }
}

void ActiveRegionGenerator::clear() noexcept
{
    if (assembler_active_region_generator_) assembler_active_region_generator_->clear();
}

// private methods

bool ActiveRegionGenerator::is_cigar_scanner(const std::string& generator) const noexcept
{
    return generator == cigar_scanner_name_;
}

bool ActiveRegionGenerator::is_assembler(const std::string& generator) const noexcept
{
    return generator == assembler_name_;
}

bool ActiveRegionGenerator::using_assembler() const noexcept
{
    return using_assembler_;
}

} // namespace coretools
} // namespace octopus
