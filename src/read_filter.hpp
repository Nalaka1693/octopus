//
//  read_filter.hpp
//  Octopus
//
//  Created by Daniel Cooke on 06/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__read_filter__
#define __Octopus__read_filter__

#include <vector>
#include <functional>
#include <algorithm>
#include <iterator>
#include <type_traits>
#include <utility>
#include <unordered_map>

#include "context_iterators.hpp"
#include "mappable_flat_multi_set.hpp"
#include "mappable_map.hpp"
#include "type_tricks.hpp"

#include <iostream> // DEBUG

namespace Octopus
{
/*
 ReadFilter stores a collection of filter functions, which can either be non-context-based (can be 
 applied to a single read), or context-based (conditional on what has already been filtered).
 
 This is a template class as the type of iterator used for context-based filteration needs to be 
 known at compile time. Essentially the class needs to know what container it is going to be operating on.
 */
template <typename BidirIt>
class ReadFilter
{
public:
    using Iterator = BidirIt;
    
    using ContextFreeFilter  = std::function<bool(const AlignedRead&)>;
    using ContextBasedFilter = std::function<bool(const AlignedRead&, Iterator, Iterator)>;
    
    ReadFilter()  = default;
    ~ReadFilter() = default;
    
    ReadFilter(const ReadFilter&)            = default;
    ReadFilter& operator=(const ReadFilter&) = default;
    ReadFilter(ReadFilter&&)                 = default;
    ReadFilter& operator=(ReadFilter&&)      = default;
    
    void register_filter(ContextFreeFilter filter);
    void register_filter(ContextBasedFilter filter);
    
    unsigned num_filters() const noexcept;
    
    // removes all failing reads, which may not be preserved
    Iterator partition(Iterator first, Iterator last) const;
    
    // copy partitions
    template <typename InputIt, typename OutputIt1, typename OutputIt2>
    std::pair<OutputIt1, OutputIt2>
    partition_copy(InputIt first, InputIt last, OutputIt1 good_reads, OutputIt2 bad_reads) const;
    
private:
    std::vector<ContextFreeFilter> noncontext_filters_;
    std::vector<ContextBasedFilter> context_filters_;
    
    bool filter_read(const AlignedRead& read, BidirIt first_good, BidirIt previous_good) const;
};

template <typename BidirIt>
void ReadFilter<BidirIt>::register_filter(ContextFreeFilter filter)
{
    noncontext_filters_.emplace_back(std::move(filter));
}

template <typename BidirIt>
void ReadFilter<BidirIt>::register_filter(ContextBasedFilter filter)
{
    context_filters_.emplace_back(std::move(filter));
}

template <typename BidirIt>
unsigned ReadFilter<BidirIt>::num_filters() const noexcept
{
    return static_cast<unsigned>(noncontext_filters_.size() + context_filters_.size());
}

namespace detail
{
    template <typename T>
    using IsContextIterator = std::is_void<typename T::value_type>;
    
    template <typename T>
    auto get_first(T first, T last, std::true_type)
    {
        return std::begin(last);
    }
    
    template <typename T>
    auto get_first(T first, T last, std::false_type)
    {
        return first;
    }
    
    template <typename T>
    auto get_first(T first, T last)
    {
        return get_first(first, last, IsContextIterator<T> {});
    }
    
    template <typename T>
    auto get_last(T first, T last, std::true_type)
    {
        using std::begin; using std::end;
        return (begin(last) != end(last)) ? std::prev(end(last)) : begin(last);
    }
    
    template <typename T>
    auto get_last(T first, T last, std::false_type)
    {
        return (first != last) ? std::prev(last) : last;
    }
    
    template <typename T>
    auto get_last(T first, T last)
    {
        return get_last(first, last, IsContextIterator<T> {});
    }
} // namespace detail

template <typename BidirIt>
BidirIt ReadFilter<BidirIt>::partition(BidirIt first, BidirIt last) const
{
    using std::all_of; using std::cbegin; using std::cend;
    
    auto it = std::remove_if(first, last,
                        [this] (const auto& read) {
                            return std::any_of(cbegin(noncontext_filters_), cend(noncontext_filters_),
                                                [&read] (const auto& filter) { return !filter(read); });
                        });
    
    auto prev = first;
    
    if (!noncontext_filters_.empty()) {
        it = std::unique(first, it, IsDuplicate {}); // TODO
    }
    
    return it;
//    
//    return std::remove_if(first, it,
//                     [this, first, &prev] (const auto& read) {
//                         auto all_good = all_of(cbegin(context_filters_), cend(context_filters_),
//                                                [&read, first, prev] (const auto& filter) {
//                                                    return filter(read, first, prev);
//                                                });
//                         
//                         if (all_good) ++prev;
//                         
//                         return all_good;
//                     });
}

template <typename BidirIt>
template <typename InputIt, typename OutputIt1, typename OutputIt2>
std::pair<OutputIt1, OutputIt2>
ReadFilter<BidirIt>::partition_copy(InputIt first, InputIt last,
                                    OutputIt1 good_reads, OutputIt2 bad_reads) const
{
    auto good_reads_last = good_reads;
    
    return std::partition_copy(first, last, good_reads, bad_reads,
        [this, good_reads, &good_reads_last] (const AlignedRead& read) {
            if (filter_read(read, detail::get_first(good_reads, good_reads_last),
                            detail::get_last(good_reads, good_reads_last))) {
            ++good_reads_last;
            return true;
        }
        return false;
    });
}

// private member methods

template <typename Iterator>
bool ReadFilter<Iterator>::filter_read(const AlignedRead& read,
                                       Iterator first_good, Iterator prev_good) const
{
    using std::cbegin; using std::cend; using std::all_of;
    
    return all_of(cbegin(noncontext_filters_), cend(noncontext_filters_),
                       [&read] (const auto& filter) { return filter(read); })
        && all_of(cbegin(context_filters_), cend(context_filters_),
                  [&read, first_good, prev_good] (const auto& filter) {
                      return filter(read, first_good, prev_good);
                  });
}

// non-member methods

template <typename ReadFilter>
std::pair<MappableFlatMultiSet<AlignedRead>, MappableFlatMultiSet<AlignedRead>>
partition_copy(const MappableFlatMultiSet<AlignedRead>& reads, const ReadFilter& filter)
{
    MappableFlatMultiSet<AlignedRead> good_reads {}, bad_reads {};
    
    good_reads.reserve(reads.size());
    bad_reads.reserve(reads.size() / 10);
    
    filter.partition_copy(std::cbegin(reads), std::cend(reads),
                          context_inserter(good_reads), context_inserter(bad_reads));
    
    good_reads.shrink_to_fit();
    bad_reads.shrink_to_fit();
    
    return std::make_pair(std::move(good_reads), std::move(bad_reads));
}

template <typename ReadFilter>
std::pair<MappableFlatMultiSet<AlignedRead>, MappableFlatMultiSet<AlignedRead>>
partition_copy(MappableFlatMultiSet<AlignedRead>&& reads, const ReadFilter& filter)
{
    MappableFlatMultiSet<AlignedRead> good_reads {}, bad_reads {};
    
    good_reads.reserve(reads.size());
    bad_reads.reserve(reads.size() / 10);
    
    filter.partition_copy(std::make_move_iterator(std::begin(reads)),
                          std::make_move_iterator(std::end(reads)),
                          context_inserter(good_reads),
                          context_inserter(bad_reads));
    
    reads.clear();
    
    good_reads.shrink_to_fit();
    bad_reads.shrink_to_fit();
    
    return std::make_pair(std::move(good_reads), std::move(bad_reads));
}

template <typename ReadFilter>
auto partition(MappableFlatMultiSet<AlignedRead>& reads, const ReadFilter& filter)
{
    return filter.partition(std::begin(reads), std::end(reads));
}

template <typename KeyType, typename ReadFilter>
std::pair<MappableMap<KeyType, AlignedRead>, MappableMap<KeyType, AlignedRead>>
partition_copy(const MappableMap<KeyType, AlignedRead>& reads, const ReadFilter& filter)
{
    MappableMap<KeyType, AlignedRead> good_reads {reads.size()}, bad_reads {reads.size()};
    
    for (const auto& sample_reads : reads) {
        auto sample_filtered_reads = filter_reads(sample_reads.second, filter);
        good_reads.emplace(sample_reads.first, std::move(sample_filtered_reads.first));
        bad_reads.emplace(sample_reads.first, std::move(sample_filtered_reads.second));
    }
    
    return std::make_pair(std::move(good_reads), std::move(bad_reads));
}

template <typename KeyType, typename ReadFilter>
std::pair<MappableMap<KeyType, AlignedRead>, MappableMap<KeyType, AlignedRead>>
partition_copy(MappableMap<KeyType, AlignedRead>&& reads, const ReadFilter& filter)
{
    MappableMap<KeyType, AlignedRead> good_reads {reads.size()}, bad_reads {reads.size()};
    
    for (auto&& sample_reads : reads) {
        auto sample_filtered_reads = filter_reads(std::move(sample_reads.second), filter);
        good_reads.emplace(sample_reads.first, std::move(sample_filtered_reads.first));
        bad_reads.emplace(sample_reads.first, std::move(sample_filtered_reads.second));
    }
    
    return std::make_pair(std::move(good_reads), std::move(bad_reads));
}

template <typename KeyType>
using PartitionPointMap = std::unordered_map<KeyType, typename MappableMap<KeyType, AlignedRead>::mapped_type::iterator>;

template <typename KeyType, typename ReadFilter>
PartitionPointMap<KeyType>
partition(MappableMap<KeyType, AlignedRead>& reads, const ReadFilter& filter)
{
    PartitionPointMap<KeyType> result {reads.size()};
    
    for (auto& p : reads) {
        result.emplace(p.first, partition(p.second, filter));
    }
    
    return result;
}

template <typename KeyType>
std::size_t erase_filtered_reads(MappableMap<KeyType, AlignedRead>& reads,
                                 const PartitionPointMap<KeyType>& partition_points)
{
    std::size_t result {0};
    
    for (auto& p : reads) {
        const auto it = partition_points.at(p.first);
        result += std::distance(it, std::end(p.second));
        p.second.erase(it, std::end(p.second));
    }
    
    return result;
}

} // namespace Octopus

#endif /* defined(__Octopus__read_filter__) */
