#ifndef CASM_clexmonte_run_analysis_functions
#define CASM_clexmonte_run_analysis_functions

#include "casm/clexmonte/run/covariance_functions.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/composition/CompositionConverter.hh"
#include "casm/monte/results/ResultsAnalysisFunction.hh"

// debugging
#include "casm/casm_io/container/stream_io.hh"

namespace CASM {
namespace clexmonte {

// ---
// These methods are used to construct results analysis functions. They are
// templated so that they can be reused. The definition documentation should
// state interface requirements for the methods to be applicable and usable in
// a particular context.
//
// Example requirements are:
// - that a conditions `monte::ValueMap` contains scalar "temperature"
// - that the method `ClexData &get_clex(SystemType &,
//   StateType const &, std::string const &key)`
//   exists for template type `SystemType` (i.e. when
//   SystemType=clexmonte::System).
// ---

/// \brief Make heat capacity analysis function ("heat_capacity")
monte::ResultsAnalysisFunction<Configuration> make_heat_capacity_f();

/// \brief Make mol_composition susceptibility analysis function
/// ("mol_susc(A,B)")
template <typename SystemType>
monte::ResultsAnalysisFunction<Configuration> make_mol_susc_f(
    std::shared_ptr<SystemType> const &system);

/// \brief Make param_composition susceptibility analysis function
/// ("param_susc(a,b)")
template <typename SystemType>
monte::ResultsAnalysisFunction<Configuration> make_param_susc_f(
    std::shared_ptr<SystemType> const &system);

/// \brief Make mol_composition thermo-chemical susceptibility
///     analysis function ("mol_thermochem_susc(S,A)")
template <typename SystemType>
monte::ResultsAnalysisFunction<Configuration> make_mol_thermochem_susc_f(
    std::shared_ptr<SystemType> const &system);

/// \brief Make param_composition thermo-chemical susceptibility
///     analysis function ("param_thermochem_susc(S,a)")
template <typename SystemType>
monte::ResultsAnalysisFunction<Configuration> make_param_thermochem_susc_f(
    std::shared_ptr<SystemType> const &system);

// --- Inline definitions ---

/// \brief Calculates `(kB * temperature * temperature) / n_unitcells`
inline double heat_capacity_normalization_constant_f(
    monte::Results<Configuration> const &results) {
  // validate initial state
  if (!results.initial_state.has_value()) {
    std::stringstream msg;
    msg << "Results analysis error: heat_capacity requires saving initial "
           "state";
    throw std::runtime_error(msg.str());
  }
  auto const &state = *results.initial_state;

  // validate temperature
  if (!state.conditions.scalar_values.count("temperature")) {
    std::stringstream msg;
    msg << "Results analysis error: heat_capacity requires temperature "
           "condition";
    throw std::runtime_error(msg.str());
  }
  double temperature = state.conditions.scalar_values.at("temperature");

  // calculate
  double n_unitcells =
      get_transformation_matrix_to_supercell(state).determinant();
  return (CASM::KB * temperature * temperature) / n_unitcells;
}

/// \brief Make heat capacity analysis function ("heat_capacity")
///
/// Notes:
/// - Requires sampling "potential_energy" (as per unit cell energy)
/// - Requires scalar condition "temperature"
/// - Requires result "initial_state"
inline monte::ResultsAnalysisFunction<Configuration> make_heat_capacity_f() {
  return make_variance_f(
      "heat_capacity",
      "Heat capacity (per unit cell) = "
      "var(potential_energy_per_unitcell)*n_unitcells/(kB*T*T)",
      "potential_energy", {"0"}, {}, heat_capacity_normalization_constant_f);
}

/// \brief Calculates `(kB * temperature) / n_unitcells`
inline std::function<double(monte::Results<Configuration> const &results)>
make_susc_normalization_constant_f(std::string name) {
  return [=](monte::Results<Configuration> const &results) -> double {
    // validate initial state
    if (!results.initial_state.has_value()) {
      std::stringstream msg;
      msg << "Results analysis error: " << name
          << " requires saving initial "
             "state";
      throw std::runtime_error(msg.str());
    }
    auto const &state = *results.initial_state;

    // validate temperature
    if (!state.conditions.scalar_values.count("temperature")) {
      std::stringstream msg;
      msg << "Results analysis error: " << name
          << " requires temperature "
             "condition";
      throw std::runtime_error(msg.str());
    }
    double temperature = state.conditions.scalar_values.at("temperature");

    // calculate
    double n_unitcells =
        get_transformation_matrix_to_supercell(state).determinant();
    return (CASM::KB * temperature) / n_unitcells;
  };
}

/// \brief Make mol_composition susceptibility analysis function
/// ("mol_susc(A,B)")
///
/// Notes:
/// - Requires sampling "mol_composition"
/// - Requires scalar condition "temperature"
/// - Requires result "initial_state"
template <typename SystemType>
monte::ResultsAnalysisFunction<Configuration> make_mol_susc_f(
    std::shared_ptr<SystemType> const &system) {
  auto const &component_names = get_composition_converter(*system).components();
  return make_covariance_f(
      "mol_susc",
      "Chemical susceptibility (per unit cell) = "
      "cov(mol_composition_i, mol_composition_j)*n_unitcells/(kB*T)",
      "mol_composition", "mol_composition", component_names, component_names,
      make_susc_normalization_constant_f("mol_susc"));
}

/// \brief Make param_composition susceptibility analysis function
/// ("param_susc(a,b)")
///
/// Notes:
/// - Requires sampling "param_composition"
/// - Requires scalar condition "temperature"
/// - Requires result "initial_state"
template <typename SystemType>
monte::ResultsAnalysisFunction<Configuration> make_param_susc_f(
    std::shared_ptr<SystemType> const &system) {
  // name param_composition components "a", "b", ... for each independent
  // composition axis
  composition::CompositionConverter const &composition_converter =
      get_composition_converter(*system);
  std::vector<std::string> component_names;
  for (Index i = 0; i < composition_converter.independent_compositions(); ++i) {
    component_names.push_back(composition_converter.comp_var(i));
  }
  return make_covariance_f(
      "param_susc",
      "Chemical susceptibility (per unit cell) = "
      "cov(param_composition_i, param_composition_j)*n_unitcells/(kB*T)",
      "param_composition", "param_composition", component_names,
      component_names, make_susc_normalization_constant_f("param_susc"));
}

/// \brief Make mol_composition thermo-chemical susceptibility
///     analysis function ("mol_thermochem_susc(S,A)")
///
/// Notes:
/// - Requires sampling "potential_energy" (as per unit cell energy)
/// - Requires sampling "mol_composition"
/// - Requires scalar condition "temperature"
/// - Requires result "initial_state"
template <typename SystemType>
monte::ResultsAnalysisFunction<Configuration> make_mol_thermochem_susc_f(
    std::shared_ptr<SystemType> const &system) {
  std::vector<std::string> first_component_names = {"S"};

  auto const &second_component_names =
      get_composition_converter(*system).components();

  return make_covariance_f(
      "mol_thermochem_susc",
      "Thermo-chemical susceptibility (per unit cell) = "
      "cov(potential_energy, mol_composition)*n_unitcells/(kB*T)",
      "potential_energy", "mol_composition", first_component_names,
      second_component_names,
      make_susc_normalization_constant_f("mol_thermochem_susc"));
}

/// \brief Make param_composition thermo-chemical susceptibility
///     analysis function ("param_thermochem_susc(S,a)")
///
/// Notes:
/// - Requires sampling "potential_energy" (as per unit cell energy)
/// - Requires sampling "param_composition"
/// - Requires scalar condition "temperature"
/// - Requires result "initial_state"
template <typename SystemType>
monte::ResultsAnalysisFunction<Configuration> make_param_thermochem_susc_f(
    std::shared_ptr<SystemType> const &system) {
  std::vector<std::string> first_component_names = {"S"};

  // name param_composition components "a", "b", ... for each independent
  // composition axis
  composition::CompositionConverter const &composition_converter =
      get_composition_converter(*system);
  std::vector<std::string> second_component_names;
  for (Index i = 0; i < composition_converter.independent_compositions(); ++i) {
    second_component_names.push_back(composition_converter.comp_var(i));
  }
  return make_covariance_f(
      "param_thermochem_susc",
      "Thermo-chemical susceptibility (per unit cell) = "
      "cov(potential_energy, param_composition)*n_unitcells/(kB*T)",
      "potential_energy", "param_composition", first_component_names,
      second_component_names,
      make_susc_normalization_constant_f("param_thermochem_susc"));
}

}  // namespace clexmonte
}  // namespace CASM

#endif