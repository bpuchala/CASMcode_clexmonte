#ifndef CASM_clexmonte_run_functions
#define CASM_clexmonte_run_functions

#include "casm/casm_io/Log.hh"
#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/misc/to_json.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/monte/MethodLog.hh"
#include "casm/monte/RunManager.hh"
#include "casm/monte/SamplingFixture.hh"
#include "casm/monte/events/OccLocation.hh"
#include "casm/monte/results/Results.hh"
#include "casm/monte/results/ResultsAnalysisFunction.hh"
#include "casm/monte/results/io/ResultsIO.hh"
#include "casm/monte/state/StateGenerator.hh"
#include "casm/monte/state/io/json/ValueMap_json_io.hh"

namespace CASM {
namespace clexmonte {

/// \brief Perform a series of runs, according to a state_generator
template <typename CalculationType>
void run_series(CalculationType &calculation,
                state_generator_type &state_generator,
                std::vector<monte::SamplingFixtureParams<config_type>> const
                    &sampling_fixture_params);

// --- Implementation ---

/// \brief Perform a series of runs, according to a state_generator
///
/// \param calculation A calculation instance, such as canonical::Canonical,
///     semi_grand_canonical::SemiGrandCanonical, or kinetic::Kinetic.
/// \param sampling_functions Functions of `state_type`, for sampling
/// \param analysis_functions Functions of `results_type`, for post-run
///     results analysis, such as `heat_capacity`
/// \param completion_check_params Specify when a single run should be
///     considered complete
/// \param state_generator A monte::StateGenerator, which produces a
///     a series of initial states
/// \param results_io A monte::ResultsIO, which outputs results
/// \param method_log A monte::MethodLog, where run status can be
///     output. Default does not write.
///
/// Requires:
/// - std::shared_ptr<system_type> CalculationType::system: Shared ptr
///   with system info
/// - CalculationType::run(...): Method to run a single calculation, see
///   canonical::Canonical<EngineType>::run for an example
/// - bool CalculationType::update_species: For occupant tracking,
///   should be true for KMC, false otherwise
template <typename CalculationType>
void run_series(CalculationType &calculation,
                state_generator_type &state_generator,
                std::vector<monte::SamplingFixtureParams<config_type>> const
                    &sampling_fixture_params) {
  auto &log = CASM::log();
  log.begin("Monte Carlo calculation series");

  monte::RunManager<config_type> run_manager(sampling_fixture_params);

  // Final states are made available to the state generator which can use them
  // to determine the next state and enable restarts
  log.indent() << "Checking for finished runs..." << std::endl;
  run_manager.read_final_states();
  log.indent() << "Found " << run_manager.final_states.size() << std::endl
               << std::endl;

  // For all states generated, prepare input and run canonical Monte Carlo
  while (!state_generator.is_complete(run_manager.final_states)) {
    // Get initial state for the next calculation
    log.indent() << "Generating next initial state..." << std::endl;
    state_type initial_state =
        state_generator.next_state(run_manager.final_states);
    log.indent() << to_json(initial_state.conditions) << std::endl;
    log.indent() << "Done" << std::endl;

    // Construct and initialize occupant tracking
    monte::Conversions const &convert =
        get_index_conversions(*calculation.system, initial_state);
    monte::OccCandidateList const &occ_candidate_list =
        get_occ_candidate_list(*calculation.system, initial_state);
    monte::OccLocation occ_location(convert, occ_candidate_list,
                                    calculation.update_species);
    occ_location.initialize(get_occupation(initial_state));

    // Run Monte Carlo at a single condition
    log.indent() << "Performing Run " << run_manager.final_states.size() + 1
                 << "..." << std::endl;
    calculation.run(initial_state, occ_location, run_manager);
    log.indent() << "Run " << run_manager.final_states.size() << " Done"
                 << std::endl;

    log.indent() << std::endl;
  }
  log.indent() << "Canonical Monte Carlo Done" << std::endl;
}

}  // namespace clexmonte
}  // namespace CASM

#endif
