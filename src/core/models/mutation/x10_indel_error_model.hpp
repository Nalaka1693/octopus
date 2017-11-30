// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef x10_indel_error_model_hpp
#define x10_indel_error_model_hpp

#include "indel_error_model.hpp"

namespace octopus {

class Haplotype;

class X10IndelErrorModel : public IndelErrorModel
{
public:
    using IndelErrorModel::PenaltyType;
    using IndelErrorModel::PenaltyVector;
    
    X10IndelErrorModel() = default;
    
    X10IndelErrorModel(const X10IndelErrorModel&)            = default;
    X10IndelErrorModel& operator=(const X10IndelErrorModel&) = default;
    X10IndelErrorModel(X10IndelErrorModel&&)                 = default;
    X10IndelErrorModel& operator=(X10IndelErrorModel&&)      = default;

private:
    static constexpr std::array<PenaltyType, 50> homopolymerErrors_ =
    {{
     60,59,48,43,38,32,28,23,20,18,16,15,14,13,12,11,10,10,9,
     9,8,8,8,7,7,7,6,6,6,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1
     }};
    static constexpr std::array<PenaltyType, 50> diNucleotideTandemRepeatErrors_ =
    {{
     60,58,47,42,37,31,27,22,19,18,16,15,14,13,12,11,10,10,10,
     9,9,8,8,7,7,7,6,6,6,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1
     }};
    static constexpr std::array<PenaltyType, 50> triNucleotideTandemRepeatErrors_ =
    {{
     60,57,46,41,36,30,28,23,20,19,17,16,15,14,13,12,11,11,10,
     9,9,8,8,7,7,7,6,6,6,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1
     }};
    static constexpr std::array<PenaltyType, 50> polyNucleotideTandemRepeatErrors_ =
    {{
     60,60,51,45,45,45,45,45,23,20,19,17,16,15,14,13,12,11,11,10,
     9,9,8,8,7,7,7,6,6,6,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1
     }};
    static constexpr std::array<PenaltyType, 15> extendPenalties_ =
    {{
     10,10,9,8,7,5,4,4,4,3,3,3,2,2,2
     }};
    
    virtual void do_set_penalities(const Haplotype& haplotype, PenaltyVector& gap_open, PenaltyVector& gap_extend) const;
};
    
} // namespace octopus

#endif
