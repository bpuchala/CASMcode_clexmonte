#include "casm/clexmonte/canonical/run.hh"

#include "casm/casm_io/Log.hh"
#include "casm/clexmonte/canonical/sampling_functions.hh"
#include "casm/clexmonte/system/OccSystem.hh"
#include "casm/clexmonte/system/enforce_composition.hh"
#include "casm/clexmonte/system/sampling_functions.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"
#include "casm/monte/Conversions.hh"
#include "casm/monte/checks/CompletionCheck.hh"
#include "casm/monte/events/OccCandidate.hh"
#include "casm/monte/methods/canonical.hh"
#include "casm/monte/results/Results.hh"
#include "casm/monte/results/io/ResultsIO.hh"
#include "casm/monte/state/StateGenerator.hh"
#include "casm/monte/state/StateSampler.hh"

namespace CASM {
namespace clexmonte {
namespace canonical {

/// \brief Run canonical Monte Carlo calculations
void run(std::shared_ptr<system_type> const &system_data,
         state_generator_type &state_generator,
         monte::StateSampler<config_type> &state_sampler,
         monte::CompletionCheck &completion_check,
         monte::ResultsIO<config_type> &results_io,
         MTRand &random_number_generator) {
  auto &log = CASM::log();
  log.begin("Cluster expansion canonical Monte Carlo");

  // Final states are made available to the state generator which can use them
  // to determine the next state
  std::vector<state_type> final_states;

  // Enable restarts: Check for a partically completed path
  log.indent() << "Checking for finished runs..." << std::endl;
  final_states = results_io.read_final_states();
  log.indent() << "Found " << final_states.size() << std::endl << std::endl;

  // For all states generated, prepare input and run canonical Monte Carlo
  while (!state_generator.is_complete(final_states)) {
    log.indent() << "Generating next initial state..." << std::endl;
    // Get initial state for the next calculation
    state_type initial_state = state_generator.next_state(final_states);
    log.indent() << "Done" << std::endl;

    // Make supercell-specific potential energy clex calculator
    // (equal to formation energy calculator now)
    clexulator::ClusterExpansion &potential_energy_clex_calculator =
        get_formation_energy_clex(*system_data, initial_state);

    // Prepare supercell-specific index conversions
    monte::Conversions convert{
        *get_shared_prim(*system_data),
        get_transformation_matrix_to_super(initial_state.configuration)};

    // Prepare list of allowed swaps -- currently using all allowed
    monte::OccCandidateList occ_candidate_list{convert};
    std::vector<monte::OccSwap> canonical_swaps =
        make_canonical_swaps(convert, occ_candidate_list);
    std::vector<monte::OccSwap> grand_canonical_swaps =
        make_grand_canonical_swaps(convert, occ_candidate_list);

    log.indent() << "Enforcing composition..." << std::endl;
    enforce_composition(get_occupation(initial_state.configuration),
                        initial_state.conditions.at("comp_n"),
                        get_composition_calculator(*system_data), convert,
                        grand_canonical_swaps, random_number_generator);
    log.indent() << "Done" << std::endl;

    // Run Monte Carlo at a single condition
    log.indent() << "Beginning run " << final_states.size() + 1 << std::endl;
    results_type result =
        monte::canonical(initial_state, potential_energy_clex_calculator,
                         convert, canonical_swaps, random_number_generator,
                         state_sampler, completion_check);
    log.indent() << "Run complete" << std::endl;

    // Store final state for state generation input
    final_states.push_back(*result.final_state);

    // Write results for this condition
    Index run_index = final_states.size();
    results_io.write(result, run_index);
  }
  log.indent() << "Canonical Monte Carlo Done" << std::endl;
}

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM
