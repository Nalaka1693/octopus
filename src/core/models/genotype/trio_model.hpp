// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef trio_model_hpp
#define trio_model_hpp

#include <vector>
#include <functional>
#include <cstddef>

#include <boost/optional.hpp>

#include "basics/trio.hpp"
#include "core/types/haplotype.hpp"
#include "population_prior_model.hpp"
#include "core/models/mutation/denovo_model.hpp"
#include "core/models/haplotype_likelihood_cache.hpp"
#include "core/types/genotype.hpp"
#include "logging/logging.hpp"

namespace octopus { namespace model {

class TrioModel
{
public:
    using GenotypeVector = std::vector<Genotype<Haplotype>>;
    
    struct Latents
    {
        struct JointProbability
        {
            using GenotypeReference = std::reference_wrapper<const Genotype<Haplotype>>;
            GenotypeReference maternal, paternal, child;
            double probability;
        };
        using JointProbabilityVector = std::vector<JointProbability>;
        JointProbabilityVector joint_genotype_probabilities;
    };
    
    struct InferredLatents
    {
        Latents posteriors;
        double log_evidence;
    };
    
    struct Options
    {
        std::size_t max_joint_genotypes;
        double max_individual_mass_loss = 1e-80, max_joint_mass_loss = 1e-200;
    };
    
    TrioModel() = delete;
    
    TrioModel(const Trio& trio,
              const PopulationPriorModel& prior_model,
              const DeNovoModel& mutation_model,
              Options options,
              boost::optional<logging::DebugLogger> debug_log = boost::none);
    
    TrioModel(const TrioModel&)            = delete;
    TrioModel& operator=(const TrioModel&) = delete;
    TrioModel(TrioModel&&)                 = delete;
    TrioModel& operator=(TrioModel&&)      = delete;
    
    ~TrioModel() = default;
    
    const PopulationPriorModel& prior_model() const noexcept;
    
    InferredLatents evaluate(const GenotypeVector& maternal_genotypes,
                             const GenotypeVector& paternal_genotypes,
                             const GenotypeVector& child_genotypes,
                             const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
    
    // Use if all samples have same ploidy
    InferredLatents evaluate(const GenotypeVector& genotypes,
                             const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
    
    InferredLatents evaluate(const GenotypeVector& genotypes,
                             std::vector<GenotypeIndex>& genotype_indices,
                             const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
    
private:
    const Trio& trio_;
    const PopulationPriorModel& prior_model_;
    const DeNovoModel& mutation_model_;
    Options options_;
    mutable boost::optional<logging::DebugLogger> debug_log_;
    
    InferredLatents evaluate_allosome(const GenotypeVector& parent_genotypes,
                                      const GenotypeVector& child_genotypes,
                                      const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
};
    
} // namespace model
} // namespace octopus

#endif
