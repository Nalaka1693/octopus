//
//  basic_caller.cpp
//  Octopus
//
//  Created by Daniel Cooke on 15/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "basic_caller.hpp"

#include <unordered_map>
#include <numeric>
#include <algorithm>

#include "common.hpp"
#include "genomic_region.hpp"
#include "read_manager.hpp"
#include "allele.hpp"
#include "variant.hpp"
#include "haplotype.hpp"
#include "genotype.hpp"
#include "haplotype_tree.hpp"
#include "search_regions.hpp"
#include "vcf_record.hpp"

#include "mappable_algorithms.hpp"
#include "variant_utils.hpp"
#include "genotype_model.hpp"
#include "population_genotype_model.hpp"

#include <iostream> // TEST

BasicVariantCaller::BasicVariantCaller(ReferenceGenome& reference, ReadManager& read_manager,
                                       ReadFilter read_filter, ReadTransform read_transform,
                                       CandidateVariantGenerator& candidate_generator)
:
VariantCaller {reference, read_manager, read_filter, read_transform, candidate_generator}
{}

GenomicRegion BasicVariantCaller::get_init_region(const GenomicRegion& region)
{
    return region;
}

GenomicRegion BasicVariantCaller::get_next_region(const GenomicRegion& current_region)
{
    return GenomicRegion {"TEST", 0, 0};
}

std::unordered_map<Haplotype, double>
get_haplotype_posteriors(const std::vector<Haplotype>& haplotypes,
                         const Octopus::GenotypeModel::SampleGenotypeProbabilities& genotype_posteriors)
{
    std::unordered_map<Haplotype, double> result {};
    result.reserve(haplotypes.size());
    
    for (const auto& haplotype : haplotypes) {
        for (const auto& genotype_posterior: genotype_posteriors) {
            if (genotype_posterior.first.contains(haplotype)) {
                result[haplotype] += genotype_posterior.second;
            }
        }
    }
    
    return result;
}

double marginalise(const Allele& allele, const std::unordered_map<Haplotype, double>& haplotype_posteriors)
{
    double result {};
    
    for (const auto& haplotype_posterior : haplotype_posteriors) {
        if (haplotype_posterior.first.contains(allele)) {
            result += haplotype_posterior.second;
        }
    }
    
    return result;
}

std::unordered_map<Allele, double>
get_allele_posteriors(const std::unordered_map<Haplotype, double>& haplotype_posteriors,
                      const std::vector<Variant>& variants)
{
    auto alleles = decompose(variants);
    
    std::unordered_map<Allele, double> result {};
    result.reserve(alleles.size());
    
    for (const auto& allele : alleles) {
        result.emplace(allele, marginalise(allele, haplotype_posteriors));
    }
    
    return result;
}

double marginalise(const Allele& allele, const Octopus::GenotypeModel::SampleGenotypeProbabilities& genotype_posteriors)
{
    double result {};
    
    for (const auto& genotype_posterior : genotype_posteriors) {
        if (std::any_of(std::cbegin(genotype_posterior.first), std::cend(genotype_posterior.first),
                        [&allele] (const auto& haplotype) {
                            return haplotype.contains(allele);
                        })) {
                            result += genotype_posterior.second;
                        }
    }
    
    return result;
}

std::unordered_map<Allele, double>
compute_sample_allele_posteriors(const Octopus::GenotypeModel::SampleGenotypeProbabilities& genotype_posteriors,
                                 const std::vector<Allele>& alleles)
{
    std::unordered_map<Allele, double> result {};
    result.reserve(alleles.size());
    
    for (const auto& allele : alleles) {
        result.emplace(allele, marginalise(allele, genotype_posteriors));
    }
    
    return result;
}

auto call_genotype(const Octopus::GenotypeModel::SampleGenotypeProbabilities& genotype_posteriors)
{
    return *std::max_element(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors),
                             [] (const auto& lhs, const auto& rhs) {
                                 return lhs.second < rhs.second;
                             });
}

using GenotypeCalls = std::unordered_map<Octopus::SampleIdType, std::pair<Genotype<Allele>, double>>;

std::vector<GenotypeCalls>
call_genotypes(const Octopus::GenotypeModel::GenotypeProbabilities& genotype_posteriors,
               const std::vector<GenomicRegion>& segments)
{
    std::vector<GenotypeCalls> result(segments.size());
    
    for (const auto sample_genotype_posteriors : genotype_posteriors) {
        const auto& sample_genotype_call = call_genotype(sample_genotype_posteriors.second);
        for (size_t i {}; i < segments.size(); ++i) {
            result[i].emplace(sample_genotype_posteriors.first,
                              std::make_pair(splice<Allele>(sample_genotype_call.first, segments[i]),
                                             sample_genotype_call.second));
        }
    }
    
    return result;
}

using AllelePosteriors = std::unordered_map<Octopus::SampleIdType, std::unordered_map<Allele, double>>;

AllelePosteriors compute_allele_posteriors(const Octopus::GenotypeModel::GenotypeProbabilities& genotype_posteriors,
                                           const std::vector<Haplotype>& haplotypes,
                                           const std::vector<Allele>& alleles)
{
    AllelePosteriors result {};
    result.reserve(genotype_posteriors.size());
    
    for (const auto sample_genotype_posteriors : genotype_posteriors) {
        auto haplotype_posteriors = get_haplotype_posteriors(haplotypes, sample_genotype_posteriors.second);
        result.emplace(sample_genotype_posteriors.first,
                       compute_sample_allele_posteriors(sample_genotype_posteriors.second, alleles));
    }
    
    return result;
}

std::vector<VcfRecord::SequenceType> to_vcf_genotype(const Genotype<Allele>& genotype)
{
    std::vector<VcfRecord::SequenceType> result {};
    result.reserve(genotype.ploidy());
    for (const auto& allele : genotype) result.push_back(allele.get_sequence());
    return result;
}

unsigned to_phred_quality(double p)
{
    return -10 * static_cast<unsigned>(std::log10(1.0 - p));
}

VcfRecord call_segment(const std::vector<Allele>& segment,
                       const GenotypeCalls& genotype_calls,
                       const AllelePosteriors& allele_posteriors)
{
    const auto& ref_allele = segment.front();
    
    auto result = VcfRecord::Builder();
    
    result.set_chromosome(get_contig_name(ref_allele));
    result.set_position(get_begin(ref_allele));
    result.set_ref_allele(ref_allele.get_sequence());
    result.set_alt_allele(segment.back().get_sequence());
    result.set_quality(60);
    result.add_info("NS", std::to_string(genotype_calls.size()));
    result.set_format({"GT", "GQ"});
    
    for (const auto& sample_call : genotype_calls) {
        const auto& sample = sample_call.first;
        result.add_genotype(sample, to_vcf_genotype(sample_call.second.first), false);
        result.add_genotype_field(sample, "GQ", std::to_string(to_phred_quality(sample_call.second.second)));
    }
    
    return result.build_once();
}

std::vector<VcfRecord> BasicVariantCaller::call_variants(const GenomicRegion& region,
                                                         const std::vector<Variant>& candidates,
                                                         const ReadMap& reads)
{
    std::vector<VcfRecord> result {};
    
    Octopus::HaplotypeTree tree {reference_};
    extend_tree(candidates, tree);
    
    auto haplotypes = tree.get_haplotypes(region);
    
    std::cout << "there are " << haplotypes.size() << " haplotypes" << std::endl;
    
    auto genotype_model = std::make_unique<Octopus::PopulationGenotypeModel>(1, 2);
    
    auto genotype_posteriors = genotype_model->evaluate(haplotypes, reads);
    
    auto alleles = decompose(candidates);
    
    auto allele_posteriors = compute_allele_posteriors(genotype_posteriors, haplotypes, alleles);
    
    auto segments = segment(alleles);
    
    std::vector<GenomicRegion> regions {};
    regions.reserve(segments.size());
    for (const auto& segment : segments) regions.push_back(segment.front().get_region());
    
    auto genotype_calls = call_genotypes(genotype_posteriors, regions);
    auto it = std::cbegin(genotype_calls);
    
    for (const auto& segment : segments) {
        result.emplace_back(call_segment(segment, *it, allele_posteriors));
        ++it;
    }
    
    return result;
}
