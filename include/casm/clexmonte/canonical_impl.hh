#ifndef CASM_clexmonte_canonical_impl
#define CASM_clexmonte_canonical_impl

#include "casm/clexmonte/canonical.hh"
#include "casm/clexmonte/run/analysis_functions.hh"
#include "casm/clexmonte/run/functions.hh"
#include "casm/clexmonte/state/Conditions.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/state/enforce_composition.hh"
#include "casm/clexmonte/state/modifying_functions.hh"
#include "casm/clexmonte/state/sampling_functions.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/monte/Conversions.hh"
#include "casm/monte/events/OccEventProposal.hh"
#include "casm/monte/events/OccLocation.hh"
#include "casm/monte/methods/occupation_metropolis.hh"
#include "casm/monte/results/Results.hh"
#include "casm/monte/state/State.hh"
#include "casm/monte/state/StateSampler.hh"

namespace CASM {
namespace clexmonte {
namespace canonical {

template <typename EngineType>
Canonical<EngineType>::Canonical(
    std::shared_ptr<system_type> _system,
    std::shared_ptr<EngineType> _random_number_engine)
    : system(_system),
      random_number_generator(_random_number_engine),
      state(nullptr),
      transformation_matrix_to_supercell(Eigen::Matrix3l::Zero(3, 3)),
      occ_location(nullptr),
      state_sampler(nullptr) {
  if (!is_clex_data(*system, "formation_energy")) {
    throw std::runtime_error(
        "Error constructing Canonical: no 'formation_energy' clex.");
  }
}

/// \brief Perform a single run, evolving current state
///
/// Notes:
/// - state and occ_location are evolved and end in modified states
template <typename EngineType>
void Canonical<EngineType>::run(state_type &state,
                                monte::OccLocation &occ_location,
                                monte::RunManager<config_type> &run_manager) {
  if (!state.conditions.scalar_values.count("temperature")) {
    throw std::runtime_error(
        "Error in Canonical::run: state `temperature` not set.");
  }
  if (!state.conditions.vector_values.count("mol_composition")) {
    throw std::runtime_error(
        "Error in Canonical::run: state `mol_composition` conditions not set.");
  }

  this->state = &state;
  this->transformation_matrix_to_supercell =
      get_transformation_matrix_to_supercell(state);
  this->occ_location = &occ_location;
  this->conditions = make_conditions(*this->system, state);

  CanonicalPotential potential(this->system);
  potential.set(this->state, this->conditions);

  /// \brief Construct swaps
  monte::Conversions const &convert = get_index_conversions(*system, state);
  monte::OccCandidateList const &occ_candidate_list =
      get_occ_candidate_list(*system, state);

  std::vector<monte::OccSwap> canonical_swaps =
      make_canonical_swaps(convert, occ_candidate_list);

  std::vector<monte::OccSwap> grand_canonical_swaps =
      make_grand_canonical_swaps(convert, occ_candidate_list);

  // Enforce composition
  clexmonte::enforce_composition(
      get_occupation(state),
      state.conditions.vector_values.at("mol_composition"),
      get_composition_calculator(*system), grand_canonical_swaps, occ_location,
      random_number_generator);

  // Run Monte Carlo at a single condition
  typedef monte::RandomNumberGenerator<EngineType> generator_type;
  monte::occupation_metropolis(state, occ_location, potential, canonical_swaps,
                               monte::propose_canonical_event<generator_type>,
                               random_number_generator, run_manager);
}

/// \brief Perform a series of runs, according to a state generator
template <typename EngineType>
void Canonical<EngineType>::run_series(
    state_generator_type &state_generator,
    std::vector<monte::SamplingFixtureParams<config_type>> const
        &sampling_fixture_params) {
  clexmonte::run_series(*this, state_generator, sampling_fixture_params);
}

/// \brief Construct functions that may be used to sample various quantities of
///     the Monte Carlo calculation as it runs
///
/// \param calculation Shared pointer to Canonical calculation, which
///     can be used by sampling functions to access system and calculation data
///     such as the prim, the cluster expansion, and the composition axes.
///
template <typename EngineType>
monte::StateSamplingFunctionMap<Configuration>
Canonical<EngineType>::standard_sampling_functions(
    std::shared_ptr<Canonical<EngineType>> const &calculation) {
  std::vector<monte::StateSamplingFunction<Configuration>> functions = {
      make_temperature_f(calculation),
      make_mol_composition_f(calculation),
      make_param_composition_f(calculation),
      make_formation_energy_corr_f(calculation),
      make_formation_energy_f(calculation),
      make_potential_energy_f(calculation)};

  monte::StateSamplingFunctionMap<Configuration> function_map;
  for (auto const &f : functions) {
    function_map.emplace(f.name, f);
  }
  return function_map;
}

/// \brief Construct functions that may be used to analyze Monte Carlo
///     calculation results
template <typename EngineType>
monte::ResultsAnalysisFunctionMap<Configuration>
Canonical<EngineType>::standard_analysis_functions(
    std::shared_ptr<Canonical<EngineType>> const &calculation) {
  std::vector<monte::ResultsAnalysisFunction<Configuration>> functions = {
      make_heat_capacity_f(calculation), make_mol_susc_f(calculation),
      make_param_susc_f(calculation), make_mol_thermochem_susc_f(calculation),
      make_param_thermochem_susc_f(calculation)};

  monte::ResultsAnalysisFunctionMap<Configuration> function_map;
  for (auto const &f : functions) {
    function_map.emplace(f.name, f);
  }
  return function_map;
}

/// \brief Construct functions that may be used to modify states
template <typename EngineType>
monte::StateModifyingFunctionMap<config_type>
Canonical<EngineType>::standard_modifying_functions(
    std::shared_ptr<Canonical<EngineType>> const &calculation) {
  std::vector<monte::StateModifyingFunction<config_type>> functions = {
      make_set_mol_composition_f(calculation)};

  monte::StateModifyingFunctionMap<config_type> function_map;
  for (auto const &f : functions) {
    function_map.emplace(f.name, f);
  }
  return function_map;
}

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM

#endif
