// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "octopus.hpp"

#include <vector>
#include <deque>
#include <queue>
#include <map>
#include <set>
#include <algorithm>
#include <numeric>
#include <memory>
#include <functional>
#include <cstddef>
#include <typeinfo>
#include <thread>
#include <future>
#include <condition_variable>
#include <mutex>
#include <atomic>
#include <chrono>
#include <sstream>
#include <iostream>
#include <cassert>

#include <boost/optional.hpp>

#include "config/common.hpp"
#include "basics/genomic_region.hpp"
#include "concepts/mappable.hpp"
#include "containers/mappable_flat_multi_set.hpp"
#include "containers/mappable_map.hpp"
#include "io/reference/reference_genome.hpp"
#include "io/read/read_manager.hpp"
#include "readpipe/read_pipe_fwd.hpp"
#include "readpipe/buffered_read_pipe.hpp"
#include "utils/mappable_algorithms.hpp"
#include "utils/read_stats.hpp"
#include "utils/append.hpp"
#include "config/octopus_vcf.hpp"
#include "core/callers/caller_factory.hpp"
#include "core/callers/caller.hpp"
#include "utils/maths.hpp"
#include "logging/progress_meter.hpp"
#include "logging/logging.hpp"
#include "logging/error_handler.hpp"
#include "core/tools/vcf_header_factory.hpp"
#include "io/variant/vcf.hpp"
#include "utils/timing.hpp"
#include "exceptions/program_error.hpp"
#include "csr/filters/variant_call_filter.hpp"
#include "csr/filters/variant_call_filter_factory.hpp"
#include "readpipe/buffered_read_pipe.hpp"

#include "timers.hpp" // BENCHMARK

namespace octopus {

using logging::get_debug_log;

namespace {

template <typename S>
void print_input_regions(S&& stream, const InputRegionMap& regions)
{
    stream << "All input regions:" << '\n';
    for (const auto& p : regions) {
        stream << "Contig " << p.first << '\n';
        for (const auto& region : p.second) {
            stream << region << ' ';
        }
        stream << '\n';
    }
}

void print_input_regions(const InputRegionMap& regions)
{
    print_input_regions(std::cout, regions);
}

bool apply_csr(const GenomeCallingComponents& components) noexcept
{
    return static_cast<bool>(components.filtered_output());
}

using CallTypeSet = std::set<std::type_index>;

VcfHeader make_vcf_header(const std::vector<SampleName>& samples,
                          const std::vector<GenomicRegion::ContigName>& contigs,
                          const ReferenceGenome& reference,
                          const CallTypeSet& call_types,
                          const std::string& command)
{
    auto builder = vcf::make_header_template().set_samples(samples);
    for (const auto& contig : contigs) {
        builder.add_contig(contig, {{"length", std::to_string(reference.contig_size(contig))}});
    }
    builder.add_basic_field("reference", reference.name());
    builder.add_structured_field("octopus", {{"command", '"' + command + '"'}});
    VcfHeaderFactory factory {};
    for (const auto& type : call_types) {
        factory.register_call_type(type);
    }
    factory.annotate(builder);
    return builder.build_once();
}

VcfHeader make_vcf_header(const std::vector<SampleName>& samples,
                          const GenomicRegion::ContigName& contig,
                          const ReferenceGenome& reference,
                          const CallTypeSet& call_types,
                          const std::string& command)
{
    return make_vcf_header(samples, std::vector<GenomicRegion::ContigName> {contig},
                           reference, call_types, command);
}

bool has_reads(const GenomicRegion& region, ContigCallingComponents& components)
{
    return components.read_manager.get().has_reads(components.samples.get(), region);
}

auto get_call_types(const GenomeCallingComponents& components, const std::vector<ContigName>& contigs)
{
    CallTypeSet result {};
    
    for (const auto& contig : components.contigs()) {
        const auto tmp_caller = components.caller_factory().make(contig);
        auto caller_call_types = tmp_caller->call_types();
        result.insert(std::begin(caller_call_types), std::end(caller_call_types));
    }
    
    return result;
}

void write_caller_output_header(GenomeCallingComponents& components, const std::string& command)
{
    const auto call_types = get_call_types(components, components.contigs());
    if (components.sites_only() && !apply_csr(components)) {
        components.output() << make_vcf_header({}, components.contigs(), components.reference(),
                                               call_types, command);
    } else {
        components.output() << make_vcf_header(components.samples(), components.contigs(),
                                               components.reference(), call_types, command);
    }
}

std::string get_caller_name(const GenomeCallingComponents& components)
{
    assert(!components.contigs().empty());
    const auto& test_contig = components.contigs().front();
    const auto test_caller = components.caller_factory().make(test_contig);
    return test_caller->name(); // There can only be one caller type per run
}

void log_startup_info(const GenomeCallingComponents& components)
{
    logging::InfoLogger log {};
    std::ostringstream ss {};
    if (!components.samples().empty()) {
        const auto num_samples = components.samples().size();
        if (num_samples == 1) {
            ss << "Detected 1 sample: ";
        } else {
            ss << "Detected " << num_samples << " samples: ";
        }
    }
    std::transform(std::cbegin(components.samples()), std::cend(components.samples()),
                   std::ostream_iterator<std::string> {ss, " "},
                   [] (const auto& sample) -> std::string {
                       return "\"" + sample + "\"";
                   });
    auto str = ss.str();
    str.pop_back(); // the extra whitespace
    log << str;
    stream(log) << "Invoked calling model: " << get_caller_name(components);
    const auto search_size = utils::format_with_commas(sum_region_sizes(components.search_regions()));
    const auto num_threads = components.num_threads();
    if (num_threads) {
        if (*num_threads == 1) {
            stream(log) << "Processing " << search_size << "bp with a single thread";
        } else {
            stream(log) << "Processing " << search_size << "bp with " << *num_threads << " threads";
        }
    } else {
        stream(log) << "Processing " << search_size << "bp with automatic thread management";
    }
    auto sl = stream(log);
    const bool is_filtered_run {components.filtered_output()};
    if (is_filtered_run) {
        sl << "Writing filtered calls to ";
    } else {
        sl << "Writing unfiltered calls to ";
    }
    auto output_path = components.output().path();
    if (is_filtered_run) {
        output_path = components.filtered_output()->path();
    }
    if (output_path) {
        sl << *output_path;
    } else {
        sl << "stdout";
    }
}

void write_calls(std::vector<VcfRecord>&& calls, VcfWriter& out)
{
    static auto debug_log = get_debug_log();
    if (debug_log) stream(*debug_log) << "Writing " << calls.size() << " calls to output";
    write(calls, out);
    calls.clear();
    calls.shrink_to_fit();
}

auto find_max_window(const ContigCallingComponents& components,
                     const GenomicRegion& remaining_call_region)
{
    const auto& rm = components.read_manager.get();
    if (!rm.has_reads(components.samples.get(), remaining_call_region)) {
        return remaining_call_region;
    }
    auto result = rm.find_covered_subregion(components.samples, remaining_call_region,
                                            components.read_buffer_size);
    if (ends_before(result, remaining_call_region)) {
        auto rest = right_overhang_region(remaining_call_region, result);
        if (!rm.has_reads(components.samples.get(), rest)) {
            result = remaining_call_region;
        }
    }
    return result;
}

auto propose_call_subregion(const ContigCallingComponents& components,
                            const GenomicRegion& remaining_call_region,
                            boost::optional<GenomicRegion::Size> min_size = boost::none)
{
    if (is_empty(remaining_call_region)) {
        return remaining_call_region;
    }
    const auto max_window = find_max_window(components, remaining_call_region);
    if (ends_before(remaining_call_region, max_window)) {
        return remaining_call_region;
    }
    if (min_size && size(max_window) < *min_size) {
        if (size(remaining_call_region) < *min_size) {
            return remaining_call_region;
        }
        return expand_rhs(head_region(max_window), *min_size);
    }
    return max_window;
}

auto propose_call_subregion(const ContigCallingComponents& components,
                            const GenomicRegion& current_subregion,
                            const GenomicRegion& input_region,
                            boost::optional<GenomicRegion::Size> min_size = boost::none)
{
    assert(contains(input_region, current_subregion));
    return propose_call_subregion(components, right_overhang_region(input_region, current_subregion), min_size);
}

void buffer_connecting_calls(std::vector<VcfRecord>& calls,
                             const GenomicRegion& next_calling_region,
                             std::vector<VcfRecord>& buffer)
{
    const auto it = std::find_if(std::begin(calls), std::end(calls),
                                 [&next_calling_region] (const auto& call) {
                                     return mapped_end(call) > next_calling_region.begin();
                                 });
    buffer.insert(std::end(buffer),
                  std::make_move_iterator(it),
                  std::make_move_iterator(std::end(calls)));
    calls.erase(it, std::end(calls));
}

void buffer_connecting_calls(const GenomicRegion& buffered_region,
                             std::vector<VcfRecord>& calls,
                             std::vector<VcfRecord>& buffer)
{
    const auto it = std::find_if_not(std::begin(calls), std::end(calls),
                                     [&buffered_region] (const auto& call) {
                                         return mapped_begin(call) < buffered_region.end();
                                     });
    buffer.insert(std::end(buffer),
                  std::make_move_iterator(std::begin(calls)),
                  std::make_move_iterator(it));
    calls.erase(std::begin(calls), it);
}

bool is_consistent(const std::deque<VcfRecord>& merged_calls)
{
    return true; // TODO
}

void resolve_connecting_calls(std::vector<VcfRecord>& old_connecting_calls,
                              std::vector<VcfRecord>& calls,
                              const ContigCallingComponents& components)
{
    using std::begin; using std::end; using std::make_move_iterator;
    
    if (!old_connecting_calls.empty()) {
        std::vector<VcfRecord> new_connecting_calls {};
        const auto old_connecting_calls_region = encompassing_region(old_connecting_calls);
        buffer_connecting_calls(old_connecting_calls_region, calls, new_connecting_calls);
        std::deque<VcfRecord> merged_calls {};
        std::set_union(begin(old_connecting_calls), end(old_connecting_calls),
                       begin(new_connecting_calls), end(new_connecting_calls),
                       std::back_inserter(merged_calls));
        old_connecting_calls.clear();
        old_connecting_calls.shrink_to_fit();
        new_connecting_calls.clear();
        new_connecting_calls.shrink_to_fit();
        if (is_consistent(merged_calls)) {
            calls.insert(begin(calls),
                         make_move_iterator(begin(merged_calls)),
                         make_move_iterator(end(merged_calls)));
        } else {
            const auto unresolved_region = encompassing_region(merged_calls);
            merged_calls.clear();
            merged_calls.shrink_to_fit();
            auto new_calls = components.caller.call(unresolved_region);
            // TODO: we need to make sure the new calls don't contain any calls
            // outside the unresolved_region, and also possibly adjust phase regions
            // in calls past unresolved_region.
            calls.insert(begin(calls),
                         make_move_iterator(begin(new_calls)),
                         make_move_iterator(end(new_calls)));
        }
    }
}

void run_octopus_on_contig(ContigCallingComponents&& components)
{
    // TODO: refactor to use connection resolution developed for multithreaded version
    static auto debug_log = get_debug_log();
    
    assert(!components.regions.empty());
    
    std::vector<VcfRecord> calls;
    std::vector<VcfRecord> connecting_calls {};
    auto input_region = components.regions.front();
    auto subregion    = propose_call_subregion(components, input_region);
    auto first_input_region      = std::cbegin(components.regions);
    const auto last_input_region = std::cend(components.regions);
    
    while (first_input_region != last_input_region && !is_empty(subregion)) {
        if (debug_log) stream(*debug_log) << "Processing subregion " << subregion;
        
        try {
            calls = components.caller.call(subregion);
        } catch(...) {
            // TODO: which exceptions can we recover from?
            throw;
        }
        resolve_connecting_calls(connecting_calls, calls, components);
        
        auto next_subregion = propose_call_subregion(components, subregion, input_region);
        
        if (is_empty(next_subregion)) {
            ++first_input_region;
            if (first_input_region != last_input_region) {
                input_region = *first_input_region;
                next_subregion = propose_call_subregion(components, input_region);
            }
        }
        assert(connecting_calls.empty());
        
        buffer_connecting_calls(calls, next_subregion, connecting_calls);
        try {
            write_calls(std::move(calls), components.output);
        } catch(...) {
            // TODO: which exceptions can we recover from?
            throw;
        }
        subregion = std::move(next_subregion);
    }
}

void run_octopus_single_threaded(GenomeCallingComponents& components)
{
    #ifdef BENCHMARK
    init_timers();
    #endif
    components.progress_meter().start();
    for (const auto& contig : components.contigs()) {
        run_octopus_on_contig(ContigCallingComponents {contig, components});
    }
    components.progress_meter().stop();
    #ifdef BENCHMARK
        print_all_timers();
    #endif
}

VcfWriter create_unique_temp_output_file(const GenomicRegion& region,
                                         const GenomeCallingComponents& components)
{
    auto path = *components.temp_directory();
    const auto& contig = region.contig_name();
    const auto begin   = std::to_string(region.begin());
    const auto end     = std::to_string(region.end());
    boost::filesystem::path file_name {contig + "_" + begin + "-" + end + "_temp.bcf"};
    path /= file_name;
    const auto call_types = get_call_types(components, {region.contig_name()});
    auto header = make_vcf_header(components.samples(), contig, components.reference(), call_types,
                                  "octopus-internal");
    return VcfWriter {std::move(path), std::move(header)};
}

VcfWriter create_unique_temp_output_file(const GenomicRegion::ContigName& contig,
                                         const GenomeCallingComponents& components)
{
    return create_unique_temp_output_file(components.reference().contig_region(contig),
                                          components);
}

using TempVcfWriterMap = std::unordered_map<ContigName, VcfWriter>;

TempVcfWriterMap make_temp_vcf_writers(const GenomeCallingComponents& components)
{
    if (!components.temp_directory()) {
        throw std::runtime_error {"Could not make temp writers"};
    }
    TempVcfWriterMap result {};
    result.reserve(components.contigs().size());
    for (const auto& contig : components.contigs()) {
        result.emplace(contig, create_unique_temp_output_file(contig, components));
    }
    return result;
}

struct Task : public Mappable<Task>
{
    GenomicRegion region;
    ExecutionPolicy policy;
    
    Task() = delete;
    
    Task(GenomicRegion region, ExecutionPolicy policy = ExecutionPolicy::seq)
    : region {std::move(region)}
    , policy {policy}
    {};
    
    const GenomicRegion& mapped_region() const noexcept { return region; }
};

std::ostream& operator<<(std::ostream& os, const Task& task)
{
    os << task.region;
    return os;
}

struct ContigOrder
{
    using ContigName = GenomicRegion::ContigName;
    
    template <typename Container>
    ContigOrder(const Container& contigs)
    : contigs_ {std::cbegin(contigs)
    , std::cend(contigs)}
    {}
    
    bool operator()(const ContigName& lhs, const ContigName& rhs) const
    {
        const auto it1 = std::find(std::cbegin(contigs_), std::cend(contigs_), lhs);
        const auto it2 = std::find(std::cbegin(contigs_), std::cend(contigs_), rhs);
        return it1 < it2;
    }
    
private:
    std::vector<ContigName> contigs_;
};

using TaskQueue = std::queue<Task>;
using TaskMap   = std::map<ContigName, TaskQueue, ContigOrder>;

auto count_tasks(const TaskMap& tasks) noexcept
{
    return std::accumulate(std::cbegin(tasks), std::cend(tasks), std::size_t {0},
                           [] (auto curr, const auto& p) noexcept { return curr + p.second.size(); });
}

struct TaskMakerSyncPacket
{
    TaskMakerSyncPacket() : batch_size_hint {1}, waiting {true}, num_tasks {0}, finished {}, all_done {false} {}
    std::condition_variable cv;
    std::mutex mutex;
    bool ready = true;
    std::atomic_uint batch_size_hint; // Only read by task maker
    std::atomic_bool waiting; // Only read by task maker
    std::atomic_uint num_tasks;
    std::unordered_map<ContigName, bool> finished;
    std::atomic_bool all_done;
};

void make_region_tasks(const GenomicRegion& region, const ContigCallingComponents& components, const ExecutionPolicy policy,
                       TaskQueue& result, TaskMakerSyncPacket& sync, const bool last_region_in_contig, const bool last_contig)
{
    static constexpr GenomicRegion::Size minTaskSize {5'000};
    std::unique_lock<std::mutex> lock {sync.mutex, std::defer_lock};
    auto subregion = propose_call_subregion(components, region, minTaskSize);
    if (ends_equal(subregion, region)) {
        lock.lock();
        sync.cv.wait(lock, [&] () { return sync.ready; });
        result.emplace(std::move(subregion), policy);
        ++sync.num_tasks;
        if (last_region_in_contig) {
            sync.finished.at(region.contig_name()) = true;
            if (last_contig) sync.all_done = true;
        }
        lock.unlock();
        sync.cv.notify_one();
    } else {
        std::deque<GenomicRegion> batch {};
        batch.push_back(subregion);
        bool done {false};
        while (true) {
            while (batch.size() < std::max(sync.batch_size_hint.load(), 1u) || !sync.waiting) {
                subregion = propose_call_subregion(components, subregion, region, minTaskSize);
                batch.push_back(subregion);
                assert(!ends_before(region, subregion));
                if (ends_equal(subregion, region)) {
                    done = true;
                    break;
                }
            }
            assert(!batch.empty());
            assert(!lock.owns_lock());
            lock.lock();
            sync.cv.wait(lock, [&] () { return sync.ready; });
            for (auto&& r : batch) result.emplace(std::move(r), policy);
            sync.num_tasks += batch.size();
            if (done) {
                if (last_region_in_contig) {
                    sync.finished.at(region.contig_name()) = true;
                    if (last_contig) sync.all_done = true;
                } else {
                    assert(!last_contig);
                }
                lock.unlock();
                sync.cv.notify_one();
                break;
            } else {
                batch.clear();
                lock.unlock();
                sync.cv.notify_one();
            }
        }
    }
}

void make_contig_tasks(const ContigCallingComponents& components, const ExecutionPolicy policy,
                       TaskQueue& result, TaskMakerSyncPacket& sync, const bool last_contig)
{
    if (components.regions.empty()) return;
    std::for_each(std::cbegin(components.regions), std::prev(std::cend(components.regions)), [&] (const auto& region) {
        make_region_tasks(region, components, policy, result, sync, false, last_contig);
    });
    make_region_tasks(components.regions.back(), components, policy, result, sync, true, last_contig);
}

ExecutionPolicy make_execution_policy(const GenomeCallingComponents& components)
{
    if (components.num_threads()) {
        return ExecutionPolicy::seq;
    }
    return ExecutionPolicy::par;
}

auto make_contig_components(const ContigName& contig, GenomeCallingComponents& components, const unsigned num_threads)
{
    ContigCallingComponents result {contig, components};
    result.read_buffer_size /= num_threads;
    return result;
}

void make_tasks_helper(TaskMap& tasks, std::vector<ContigName> contigs, GenomeCallingComponents& components,
                       const unsigned num_threads, ExecutionPolicy execution_policy, TaskMakerSyncPacket& sync)
{
    try {
        static auto debug_log = get_debug_log();
        assert(!contigs.empty());
        std::for_each(std::cbegin(contigs), std::prev(std::cend(contigs)), [&] (const auto& contig) {
            auto contig_components = make_contig_components(contig, components, num_threads);
            make_contig_tasks(contig_components, execution_policy, tasks[contig], sync, false);
            if (debug_log) *debug_log << "Finished making tasks for contig " << contig;
        });
        auto contig_components = make_contig_components(contigs.back(), components, num_threads);
        make_contig_tasks(contig_components, execution_policy, tasks[contigs.back()], sync, true);
        if (debug_log) *debug_log << "Finished making tasks";
    } catch (const Error& e) {
        log_error(e);
        logging::FatalLogger fatal_log {};
        fatal_log << "Encountered error in task maker thread. Calling terminate";
        std::terminate();
    } catch (const std::exception& e) {
        log_error(e);
        logging::FatalLogger fatal_log {};
        fatal_log << "Encountered error in task maker thread. Calling terminate";
        std::terminate();
    } catch (...) {
        logging::FatalLogger fatal_log {};
        fatal_log << "Encountered error in task maker thread. Calling terminate";
        std::terminate();
    }
}

std::thread make_tasks(TaskMap& tasks, GenomeCallingComponents& components, const unsigned num_threads,
                       TaskMakerSyncPacket& sync)
{
    auto contigs = components.contigs();
    if (contigs.empty()) {
        sync.all_done = true;
        return std::thread {};
    }
    sync.finished.reserve(contigs.size());
    for (const auto& contig : contigs) {
        sync.finished.emplace(contig, false);
    }
    return std::thread {make_tasks_helper, std::ref(tasks), std::move(contigs), std::ref(components),
                        num_threads, make_execution_policy(components), std::ref(sync)};
}

void log_num_cores(const unsigned num_cores)
{
    auto debug_log = logging::get_debug_log();
    if (debug_log) stream(*debug_log) << "Detected " << num_cores << " system cores";
}

void warn_undetected_cores()
{
    logging::WarningLogger log {};
    log << "Unable to detect the number of system cores,"
    " it may be better to run with a user number if the number of cores is known";
}

unsigned calculate_num_task_threads(const GenomeCallingComponents& components)
{
    if (components.num_threads()) {
        return *components.num_threads();
    }
    // TODO: come up with a better calculation
    const auto num_cores = std::thread::hardware_concurrency();
    if (num_cores > 0) {
        log_num_cores(num_cores);
        return num_cores;
    } else {
        warn_undetected_cores();
        return std::min(components.read_manager().num_files(), 8u);
    }
}

Task pop(TaskMap& tasks, TaskMakerSyncPacket& sync)
{
    assert(!tasks.empty());
    std::unique_lock<std::mutex> lock {sync.mutex};
    sync.ready = false;
    assert(sync.num_tasks > 0);
    const auto contig_task_itr = std::begin(tasks);
    assert(!contig_task_itr->second.empty());
    const auto result = std::move(contig_task_itr->second.front());
    contig_task_itr->second.pop();
    if (sync.finished.at(contig_task_itr->first) && contig_task_itr->second.empty()) {
        static auto debug_log = get_debug_log();
        if (debug_log) stream(*debug_log) << "Finished calling contig " << contig_task_itr->first;
        tasks.erase(contig_task_itr);
    }
    --sync.num_tasks;
    sync.ready = true;
    lock.unlock();
    sync.cv.notify_one();
    return result;
}

struct CompletedTask : public Task
{
    CompletedTask(Task task) : Task {std::move(task)}, calls {}, runtime {} {}
    std::vector<VcfRecord> calls;
    utils::TimeInterval runtime;
};

std::string duration(const CompletedTask& task)
{
    std::ostringstream ss {};
    ss << task.runtime;
    return ss.str();
}

std::ostream& operator<<(std::ostream& os, const CompletedTask& task)
{
    os << task.region;
    return os;
}

struct CallerSyncPacket
{
    CallerSyncPacket() : num_finished {0} {}
    std::condition_variable cv;
    std::mutex mutex;
    std::atomic_uint num_finished;
};

template<typename R>
bool is_ready(const std::future<R>& f)
{
    return f.valid() && f.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
}

auto run(Task task, ContigCallingComponents components, CallerSyncPacket& sync)
{
    static auto debug_log = get_debug_log();
    if (debug_log) stream(*debug_log) << "Spawning task " << task;
    return std::async(std::launch::async, [task = std::move(task), components = std::move(components), &sync] () {
        try {
            CompletedTask result {task};
            result.runtime.start = std::chrono::system_clock::now();
            result.calls = components.caller.call(task.region);
            result.runtime.end = std::chrono::system_clock::now();
            std::unique_lock<std::mutex> lock {sync.mutex};
            ++sync.num_finished;
            sync.cv.notify_all();
            lock.unlock();
            return result;
        } catch (...) {
            logging::ErrorLogger error_log {};
            stream(error_log) << "Encountered a problem whilst calling " << task;
            using namespace std::chrono_literals;
            std::this_thread::sleep_for(2s); // Try to make sure the error is logged before raising
            std::unique_lock<std::mutex> lock {sync.mutex};
            ++sync.num_finished;
            lock.unlock();
            sync.cv.notify_all();
            throw;
        }
    });
}

using CompletedTaskMap = std::map<ContigName, std::map<ContigRegion, CompletedTask>>;
using HoldbackTask = boost::optional<std::reference_wrapper<const CompletedTask>>;

auto get_writable_completed_tasks(CompletedTask&& task, CompletedTaskMap::mapped_type& buffered_tasks,
                                  TaskQueue& running_tasks, HoldbackTask& holdback)
{
    std::deque<CompletedTask> result {std::move(task)};
    while (!running_tasks.empty()) {
        const auto itr = buffered_tasks.find(contig_region(running_tasks.front()));
        if (itr != std::end(buffered_tasks)) {
            result.push_back(std::move(itr->second));
            buffered_tasks.erase(itr);
            running_tasks.pop();
        } else {
            break;
        }
    }
    if (holdback) {
        const auto itr = buffered_tasks.find(contig_region(holdback->get()));
        assert(itr->second < result.front());
        result.push_front(std::move(itr->second));
        buffered_tasks.erase(itr);
        holdback = boost::none;
    }
    return result;
}

using ContigCallingComponentFactory    = std::function<ContigCallingComponents()>;
using ContigCallingComponentFactoryMap = std::map<ContigName, ContigCallingComponentFactory>;

auto make_contig_calling_component_factory_map(GenomeCallingComponents& components)
{
    ContigCallingComponentFactoryMap result {};
    for (const auto& contig : components.contigs()) {
        result.emplace(contig, [&components, contig] () -> ContigCallingComponents
                       { return ContigCallingComponents {contig, components}; });
    }
    return result;
}

auto find_first_lhs_connecting(const std::vector<VcfRecord>& lhs_calls, const GenomicRegion& rhs_region)
{
    const auto rhs_begin = mapped_begin(rhs_region);
    return std::find_if(std::cbegin(lhs_calls), std::cend(lhs_calls),
                        [&rhs_begin] (const auto& call) { return mapped_end(call) > rhs_begin; });
}

auto find_last_rhs_connecting(const GenomicRegion& lhs_region, const std::vector<VcfRecord>& rhs_calls)
{
    const auto lhs_end = mapped_end(lhs_region);
    return std::find_if_not(std::cbegin(rhs_calls), std::cend(rhs_calls),
                            [&lhs_end] (const auto& call) { return mapped_begin(call) < lhs_end; });
}

void resolve_connecting_calls(CompletedTask& lhs, CompletedTask& rhs,
                              const ContigCallingComponentFactory& calling_components)
{
    static auto debug_log = get_debug_log();
    using std::begin; using std::end; using std::cbegin; using std::cend; using std::make_move_iterator;
    if (lhs.calls.empty() || rhs.calls.empty()) return;
    const auto first_lhs_connecting = find_first_lhs_connecting(lhs.calls, encompassing_region(rhs.calls));
    const auto last_rhs_connecting  = find_last_rhs_connecting(encompassing_region(lhs.calls), rhs.calls);
    if (first_lhs_connecting == cend(lhs.calls) && last_rhs_connecting == cbegin(rhs.calls)) {
        return;
    }
    if (debug_log) {
        stream(*debug_log) << "Resolving connecting calls between tasks " << lhs << " & " << rhs;
    }
    std::deque<VcfRecord> merged_calls {};
    std::set_union(first_lhs_connecting, cend(lhs.calls),
                   cbegin(rhs.calls), last_rhs_connecting,
                   std::back_inserter(merged_calls));
    lhs.calls.erase(first_lhs_connecting, cend(lhs.calls));
    rhs.calls.erase(cbegin(rhs.calls), last_rhs_connecting);
    if (is_consistent(merged_calls)) {
        rhs.calls.insert(begin(rhs.calls),
                         make_move_iterator(begin(merged_calls)),
                         make_move_iterator(end(merged_calls)));
    } else {
        const auto unresolved_region = encompassing_region(mapped_region(merged_calls.front()),
                                                           mapped_region(merged_calls.back()));
        const auto components = calling_components();
        auto num_unresolved_region_reads = components.read_manager.get().count_reads(components.samples, unresolved_region);
        if (num_unresolved_region_reads <= components.read_buffer_size) {
            merged_calls.clear();
            merged_calls.shrink_to_fit();
            if (debug_log) {
                stream(*debug_log) << "Calls are inconsistent in connecting region " << unresolved_region
                                   << ". Recalling the region";
            }
            logging::WarningLogger warn_log {};
            stream(warn_log) << "Recalling " << unresolved_region
                             << " due to call inconsistency between thread tasks. This may increase expected runtime";
            auto resolved_calls = components.caller.call(unresolved_region);
            if (!resolved_calls.empty()) {
                if (!contains(unresolved_region, encompassing_region(resolved_calls))) {
                    // TODO
                }
                // TODO: we may need to adjust phase regions in calls past the unresolved_region
                rhs.calls.insert(begin(rhs.calls),
                                 make_move_iterator(begin(resolved_calls)),
                                 make_move_iterator(end(resolved_calls)));
            }
        } else {
            // TODO: we could try to manually resolve the calls. Very difficult.
            logging::WarningLogger log {};
            stream(log) << "Skipping region " << unresolved_region
                        << " as there are too many reads to analyse the whole region, and partitions give inconsistent calls";
        }
    }
}

void resolve_connecting_calls(std::deque<CompletedTask>& adjacent_tasks,
                              const ContigCallingComponentFactory& calling_components)
{
    if (adjacent_tasks.size() < 2) return;
    assert(std::is_sorted(std::cbegin(adjacent_tasks), std::cend(adjacent_tasks)));
    auto lhs = std::begin(adjacent_tasks);
    std::for_each(std::next(lhs), std::end(adjacent_tasks),
                  [&] (auto& rhs) { resolve_connecting_calls(*lhs++, rhs, calling_components); });
}

struct TaskWriterSync
{
    std::condition_variable cv;
    std::mutex mutex;
    std::deque<CompletedTask> tasks = {};
    bool done = false;
};

void write(std::deque<CompletedTask>& tasks, TempVcfWriterMap& writers)
{
    static auto debug_log = get_debug_log();
    for (auto&& task : tasks) {
        if (debug_log) {
            stream(*debug_log) << "Writing completed task " << task << " that finished in " << duration(task);
        }
        auto& writer = writers.at(contig_name(task));
        write_calls(std::move(task.calls), writer);
    }
    tasks.clear();
}

void write_temp_vcf_helper(TempVcfWriterMap& writers, TaskWriterSync& sync)
{
    try {
        std::unique_lock<std::mutex> lock {sync.mutex, std::defer_lock};
        std::deque<CompletedTask> buffer {};
        while (!sync.done) {
            lock.lock();
            sync.cv.wait(lock, [&] () { return !sync.tasks.empty(); });
            assert(buffer.empty());
            std::swap(sync.tasks, buffer);
            lock.unlock();
            sync.cv.notify_one();
            write(buffer, writers);
        }
    } catch (const Error& e) {
        log_error(e);
        logging::FatalLogger fatal_log {};
        fatal_log << "Encountered error in task writer thread. Calling terminate";
        std::terminate();
    
    } catch (const std::exception& e) {
        log_error(e);
        logging::FatalLogger fatal_log {};
        fatal_log << "Encountered error in task writer thread. Calling terminate";
        std::terminate();
    } catch (...) {
        logging::FatalLogger fatal_log {};
        fatal_log << "Encountered error in task writer thread. Calling terminate";
        std::terminate();
    }
}

std::thread make_task_writer_thread(TempVcfWriterMap& temp_writers, TaskWriterSync& writer_sync)
{
    return std::thread {write_temp_vcf_helper, std::ref(temp_writers), std::ref(writer_sync)};
}

void write(std::deque<CompletedTask>&& tasks, VcfWriter& temp_vcf)
{
    static auto debug_log = get_debug_log();
    for (auto&& task : tasks) {
        if (debug_log) {
            stream(*debug_log) << "Writing completed task " << task << " that finished in " << duration(task);
        }
        write_calls(std::move(task.calls), temp_vcf);
    }
}

void write(std::deque<CompletedTask>&& tasks, TaskWriterSync& sync)
{
    std::unique_lock<std::mutex> lock {sync.mutex};
    utils::append(std::move(tasks), sync.tasks);
    lock.unlock();
    sync.cv.notify_one();
}

// A CompletedTask can only be written if all proceeding tasks have completed (either written or buffered)
void write_or_buffer(CompletedTask&& task, CompletedTaskMap::mapped_type& buffered_tasks,
                     TaskQueue& running_tasks, HoldbackTask& holdback,
                     TaskWriterSync& sync, const ContigCallingComponentFactory& calling_components)
{
    static auto debug_log = get_debug_log();
    if (is_same_region(task, running_tasks.front())) {
        running_tasks.pop();
        auto writable_tasks = get_writable_completed_tasks(std::move(task), buffered_tasks, running_tasks, holdback);
        assert(holdback == boost::none);
        resolve_connecting_calls(writable_tasks, calling_components);
        // Keep the last task buffered to enable connection resolution when the next task finishes
        const auto p = buffered_tasks.emplace(contig_region(writable_tasks.back()), std::move(writable_tasks.back()));
        assert(p.second);
        holdback = p.first->second;
        if (debug_log) stream(*debug_log) << "Holding back completed task " << *holdback;
        writable_tasks.pop_back();
        write(std::move(writable_tasks), sync);
    } else {
        if (debug_log) stream(*debug_log) << "Buffering completed task " << task;
        buffered_tasks.emplace(contig_region(task), std::move(task));
    }
}

void wait_until_finished(TaskWriterSync& sync)
{
    std::unique_lock<std::mutex> lock {sync.mutex};
    sync.cv.wait(lock, [&] () { return sync.tasks.empty(); });
    sync.done = true;
    lock.unlock();
    sync.cv.notify_one();
}

using FutureCompletedTasks = std::vector<std::future<CompletedTask>>;
using RemainingTaskMap = std::map<ContigName, std::deque<CompletedTask>>;

void extract_remaining_future_tasks(FutureCompletedTasks& futures, std::deque<CompletedTask>& result)
{
    const auto itr = std::remove_if(std::begin(futures), std::end(futures), [] (const auto& f) { return !f.valid(); });
    std::transform(std::begin(futures), itr, std::back_inserter(result), [] (auto& fut) { return fut.get(); });
    futures.clear();
    futures.shrink_to_fit();
}

void extract_buffered_tasks(CompletedTaskMap& buffered_tasks, std::deque<CompletedTask>& result)
{
    for (auto& p : buffered_tasks) {
        std::transform(std::make_move_iterator(std::begin(p.second)), std::make_move_iterator(std::end(p.second)),
                       std::back_inserter(result), [] (auto&& p) { return std::move(p.second); });
        p.second.clear();
    }
    buffered_tasks.clear();
}

void sort_ignoring_contig_name(std::deque<CompletedTask>& tasks)
{
    std::sort(std::begin(tasks), std::end(tasks),
              [] (const auto& lhs, const auto& rhs) { return contig_region(lhs) < contig_region(rhs); });
}

RemainingTaskMap make_map(std::deque<CompletedTask>& tasks)
{
    sort_ignoring_contig_name(tasks);
    RemainingTaskMap result {};
    for (auto&& task : tasks) {
        result[contig_name(task.region)].push_back(std::move(task));
    }
    return result;
}

RemainingTaskMap extract_remaining_tasks(FutureCompletedTasks& futures, CompletedTaskMap& buffered_tasks)
{
    std::deque<CompletedTask> tasks {};
    extract_remaining_future_tasks(futures, tasks);
    extract_buffered_tasks(buffered_tasks, tasks);
    return make_map(tasks);
}

void resolve_connecting_calls(RemainingTaskMap& remaining_tasks,
                              const ContigCallingComponentFactoryMap& calling_components)
{
    for (auto& p : remaining_tasks) {
        resolve_connecting_calls(p.second, calling_components.at(p.first));
    }
}

void write(RemainingTaskMap&& remaining_tasks, TempVcfWriterMap& temp_vcfs)
{
    for (auto& p : remaining_tasks) {
        write(std::move(p.second), temp_vcfs.at(p.first));
    }
}

void write_remaining_tasks(FutureCompletedTasks& futures, CompletedTaskMap& buffered_tasks, TempVcfWriterMap& temp_vcfs,
                           const ContigCallingComponentFactoryMap& calling_components)
{
    auto remaining_tasks = extract_remaining_tasks(futures, buffered_tasks);
    resolve_connecting_calls(remaining_tasks, calling_components);
    write(std::move(remaining_tasks), temp_vcfs);
}

auto extract_writers(TempVcfWriterMap&& vcfs)
{
    std::vector<VcfWriter> result {};
    result.reserve(vcfs.size());
    for (auto&& p : vcfs) {
        result.push_back(std::move(p.second));
    }
    vcfs.clear();
    return result;
}

auto extract_as_readers(TempVcfWriterMap&& vcfs)
{
    return writers_to_readers(extract_writers(std::move(vcfs)));
}

void merge(TempVcfWriterMap&& temp_vcf_writers, GenomeCallingComponents& components)
{
    auto temp_readers = extract_as_readers(std::move(temp_vcf_writers));
    merge(temp_readers, components.output(), components.contigs());
}

void run_octopus_multi_threaded(GenomeCallingComponents& components)
{
    using namespace std::chrono_literals;
    static auto debug_log = get_debug_log();
    
    const auto num_task_threads = calculate_num_task_threads(components);
    
    TaskMap pending_tasks {components.contigs()};
    TaskMakerSyncPacket task_maker_sync {};
    task_maker_sync.batch_size_hint = 2 * num_task_threads;
    std::unique_lock<std::mutex> pending_task_lock {task_maker_sync.mutex, std::defer_lock};
    auto maker_thread = make_tasks(pending_tasks, components, num_task_threads, task_maker_sync);
    if (maker_thread.joinable()) maker_thread.detach();
    
    FutureCompletedTasks futures(num_task_threads);
    TaskMap running_tasks {ContigOrder {components.contigs()}};
    CompletedTaskMap buffered_tasks {};
    std::map<ContigName, HoldbackTask> holdbacks {};
    // Populate all the maps first so we can make unchecked accesses
    for (const auto& contig : components.contigs()) {
        running_tasks.emplace(contig, TaskMap::mapped_type {});
        buffered_tasks.emplace(contig, CompletedTaskMap::mapped_type {});
        holdbacks.emplace(contig, boost::none);
    }
    
    CallerSyncPacket caller_sync {};
    const auto calling_components = make_contig_calling_component_factory_map(components);
    unsigned num_idle_futures {0};
    
    auto temp_writers = make_temp_vcf_writers(components);
    TaskWriterSync task_writer_sync {};
    auto task_writer_thread = make_task_writer_thread(temp_writers, task_writer_sync);
    task_writer_thread.detach();
    
    // Wait for the first task to be made
    const auto tasks_available = [&] () noexcept { return task_maker_sync.num_tasks > 0; };
    while(task_maker_sync.num_tasks == 0) {
        pending_task_lock.lock();
        task_maker_sync.cv.wait(pending_task_lock, tasks_available);
        pending_task_lock.unlock();
    }
    task_maker_sync.batch_size_hint = num_task_threads / 2;
    
    components.progress_meter().start();
    
    while (!task_maker_sync.all_done || task_maker_sync.num_tasks > 0) {
        pending_task_lock.lock();
        assert(count_tasks(pending_tasks) == task_maker_sync.num_tasks);
        if (!task_maker_sync.all_done && task_maker_sync.num_tasks == 0) {
            task_maker_sync.batch_size_hint = std::max(num_idle_futures, num_task_threads / 2);
            if (num_idle_futures < futures.size()) {
                // If there are running futures then it's good periodically check to see if
                // any have finished and process them while we wait for the task maker.
                while (task_maker_sync.num_tasks == 0 && caller_sync.num_finished == 0) {
                    auto now = std::chrono::system_clock::now();
                    task_maker_sync.cv.wait_until(pending_task_lock, now + 5s, tasks_available);
                }
            } else {
                task_maker_sync.cv.wait(pending_task_lock, tasks_available);
            }
        }
        pending_task_lock.unlock();
        num_idle_futures = 0;
        for (auto& future : futures) {
            if (is_ready(future)) {
                auto completed_task = future.get();
                const auto& contig = contig_name(completed_task.region);
                write_or_buffer(std::move(completed_task), buffered_tasks.at(contig),
                                running_tasks.at(contig), holdbacks.at(contig),
                                task_writer_sync, calling_components.at(contig));
                --caller_sync.num_finished;
            }
            if (!future.valid()) {
                pending_task_lock.lock();
                if (task_maker_sync.num_tasks > 0) {
                    pending_task_lock.unlock(); // As pop will need to lock the mutex too == deadlock
                    auto task = pop(pending_tasks, task_maker_sync);
                    future = run(task, calling_components.at(contig_name(task))(), caller_sync);
                    running_tasks.at(contig_name(task)).push(std::move(task));
                } else {
                    pending_task_lock.unlock();
                    ++num_idle_futures;
                }
            }
        }
        // If there are no idle futures then all threads are busy and we must wait for one to finish,
        // otherwise we must have run out of tasks, so we should wait for new ones.
        if (num_idle_futures == 0 && caller_sync.num_finished == 0) {
            task_maker_sync.waiting = false;
            std::unique_lock<std::mutex> lock {caller_sync.mutex};
            caller_sync.cv.wait(lock, [&] () { return caller_sync.num_finished > 0; });
            task_maker_sync.waiting = true;
        } else {
            if (debug_log) stream(*debug_log) << "There are " << num_idle_futures << " idle futures";
        }
    }
    assert(task_maker_sync.num_tasks == 0);
    assert(pending_tasks.empty());
    running_tasks.clear();
    holdbacks.clear(); // holdbacks are just references to buffered tasks
    wait_until_finished(task_writer_sync);
    write_remaining_tasks(futures, buffered_tasks, temp_writers, calling_components);
    components.progress_meter().stop();
    merge(std::move(temp_writers), components);
}

} // namespace

bool is_multithreaded(const GenomeCallingComponents& components)
{
    return !components.num_threads() || *components.num_threads() > 1;
}

void run_calling(GenomeCallingComponents& components)
{
    if (is_multithreaded(components)) {
        if (DEBUG_MODE) {
            logging::WarningLogger warn_log {};
            warn_log << "Running in parallel mode can make debug log difficult to interpret";
        }
        run_octopus_multi_threaded(components);
    } else {
        run_octopus_single_threaded(components);
    }
}

void destroy(VcfWriter& writer)
{
    VcfWriter tmp {};
    swap(writer, tmp);
}

void log_filtering_info(const GenomeCallingComponents& components)
{
    logging::InfoLogger log {};
    log << "Starting Call Set Refinement (CSR) filtering";
}

std::vector<GenomicRegion> extract_call_regions(VcfReader& vcf)
{
    std::deque<GenomicRegion> regions {};
    auto p = vcf.iterate(VcfReader::UnpackPolicy::sites);
    std::transform(std::move(p.first), std::move(p.second), std::back_inserter(regions),
                   [] (const VcfRecord& record) { return mapped_region(record); });
    return {std::make_move_iterator(std::begin(regions)), std::make_move_iterator(std::end(regions))};
}

std::vector<GenomicRegion> extract_call_regions(boost::filesystem::path vcf_path)
{
    VcfReader tmp {std::move(vcf_path)};
    return extract_call_regions(tmp);
}

template <typename Map>
std::size_t sum_mapped_container_size(const Map& map)
{
    return std::accumulate(std::cbegin(map), std::cend(map), std::size_t {0},
                           [] (auto curr, const auto& p) noexcept { return curr + p.second.size(); });
}

std::vector<GenomicRegion> flatten(const InputRegionMap& regions)
{
    std::vector<GenomicRegion> result {};
    result.reserve(sum_mapped_container_size(regions));
    for (const auto& p : regions) {
        std::copy(std::cbegin(p.second), std::cend(p.second), std::back_inserter(result));
    }
    return result;
}

bool use_unfiltered_call_region_hints_for_filtering(const GenomeCallingComponents& components)
{
    // TODO: no need to do this if reads won't be used, or if likely hints won't help because the calls
    // are very dense.
    return true;
}

void run_filtering(GenomeCallingComponents& components)
{
    if (apply_csr(components)) {
        log_filtering_info(components);
        ProgressMeter progress {components.search_regions()};
        const auto& filter_factory = components.call_filter_factory();
        const auto& filter_read_pipe = components.filter_read_pipe();
        auto unfiltered_output_path = components.output().path();
        assert(unfiltered_output_path); // cannot be stdout
        BufferedReadPipe buffered_rp {filter_read_pipe, components.read_buffer_size()};
        if (use_unfiltered_call_region_hints_for_filtering(components)) {
            buffered_rp.hint(extract_call_regions(*unfiltered_output_path));
        } else {
            buffered_rp.hint(flatten(components.search_regions()));
        }
        VariantCallFilter::OutputOptions output_config {};
        if (components.sites_only()) {
            output_config.emit_sites_only = true;
        }
        const auto filter = filter_factory.make(components.reference(), std::move(buffered_rp), output_config, progress);
        assert(filter);
        const VcfReader in {std::move(*unfiltered_output_path)};
        VcfWriter& out {*components.filtered_output()};
        progress.start();
        filter->filter(in, out);
        progress.stop();
    }
}

void convert_to_legacy(const boost::filesystem::path& src, const boost::filesystem::path& dest)
{
    const VcfReader in {src};
    VcfWriter out {dest};
    convert_to_legacy(in, out);
}

VcfWriter& get_final_output(GenomeCallingComponents& components)
{
    if (apply_csr(components)) {
        return *components.filtered_output();
    } else {
        return components.output();
    }
}

void run_legacy_generation(GenomeCallingComponents& components)
{
    if (components.legacy()) {
        VcfWriter& final_output {get_final_output(components)};
        const auto output_path = final_output.path();
        if (output_path) {
            destroy(final_output);
            convert_to_legacy(*output_path, *components.legacy());
        }
    }
}

void log_run_start(const GenomeCallingComponents& components)
{
    static auto debug_log = get_debug_log();
    log_startup_info(components);
    if (debug_log) print_input_regions(stream(*debug_log), components.search_regions());
}

class CallingBug : public ProgramError
{
    std::string do_where() const override { return "run_octopus"; }
    std::string do_why() const override
    {
        if (what_) {
            return "Encountered an exception during calling '" + *what_ + "'. This means there is a bug"
            " and your results are untrustworthy.";
        } else {
            return "Encountered an unknown error during calling. This means there is a bug"
            " and your results are untrustworthy.";
        }
    }
    
    boost::optional<std::string> what_;
public:
    CallingBug() = default;
    CallingBug(const std::exception& e) : what_ {e.what()} {}
};

void run_octopus(GenomeCallingComponents& components, std::string command)
{
    static auto debug_log = get_debug_log();
    logging::InfoLogger info_log {};
    using utils::TimeInterval;
    
    log_run_start(components);
    write_caller_output_header(components, command);
    const auto start = std::chrono::system_clock::now();
    try {
        run_calling(components);
    } catch (const ProgramError& e) {
        try {
            if (debug_log) *debug_log << "Encountered an error whilst calling, attempting to cleanup";
            cleanup(components);
        } catch (...) {}
        throw;
    } catch (const std::exception& e) {
        try {
            if (debug_log) *debug_log << "Encountered an error whilst calling, attempting to cleanup";
            cleanup(components);
        } catch (...) {}
        throw CallingBug {e};
    } catch (...) {
        try {
            if (debug_log) *debug_log << "Encountered an error whilst calling, attempting to cleanup";
            cleanup(components);
        } catch (...) {}
        throw CallingBug {};
    }
    components.output().close();
    try {
        run_filtering(components);
    } catch (...) {
        try {
            if (debug_log) *debug_log << "Encountered an error whilst filtering, attempting to cleanup";
            cleanup(components);
        } catch (...) {}
        throw CallingBug {};
    }
    try {
        run_legacy_generation(components);
    } catch (...) {
        logging::WarningLogger warn_log {};
        warn_log << "Failed to make legacy vcf";
    }
    const auto end = std::chrono::system_clock::now();
    const auto search_size = sum_region_sizes(components.search_regions());
    stream(info_log) << "Finished calling "
                     << utils::format_with_commas(search_size) << "bp, total runtime "
                     << TimeInterval {start, end};
    cleanup(components);
}

} // namespace octopus
