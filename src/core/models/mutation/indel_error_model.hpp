// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef indel_error_model_hpp
#define indel_error_model_hpp

#include <vector>
#include <array>
#include <cstdint>

namespace octopus {

class Haplotype;

class IndelErrorModel
{
public:
    using PenaltyType = std::int8_t;
    using PenaltyVector = std::vector<PenaltyType>;
    
    IndelErrorModel() = default;
    
    IndelErrorModel(const IndelErrorModel&)            = default;
    IndelErrorModel& operator=(const IndelErrorModel&) = default;
    IndelErrorModel(IndelErrorModel&&)                 = default;
    IndelErrorModel& operator=(IndelErrorModel&&)      = default;
    
    virtual ~IndelErrorModel() = default;
    
    void set_penalities(const Haplotype& haplotype, PenaltyVector& gap_open, PenaltyVector& gap_extend) const;
    
private:
    virtual void do_set_penalities(const Haplotype& haplotype, PenaltyVector& gap_open, PenaltyVector& gap_extend) const = 0;
};

} // namespace octopus

#endif
