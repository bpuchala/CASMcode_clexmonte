#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/methods/occupation_metropolis.hh"
#include "casm/clexmonte/monte_calculator/BaseMonteCalculator.hh"
#include "casm/clexmonte/monte_calculator/MonteCalculator.hh"
#include "casm/clexmonte/monte_calculator/analysis_functions.hh"
#include "casm/clexmonte/monte_calculator/modifying_functions.hh"
#include "casm/clexmonte/monte_calculator/sampling_functions.hh"
#include "casm/clexmonte/run/functions.hh"
#include "casm/clexmonte/state/enforce_composition.hh"
#include "casm/configuration/io/json/Configuration_json_io.hh"
#include "casm/monte/events/OccEventProposal.hh"
#include "casm/monte/sampling/RequestedPrecisionConstructor.hh"

namespace CASM {
namespace clexmonte {

/// \brief Propose and apply canonical events
class CanonicalEventGenerator {
 public:
  typedef BaseMonteCalculator::engine_type engine_type;

  /// \brief Constructor
  ///
  /// Notes:
  /// - `_canonical_swaps` should have size != 0
  ///
  /// \param _canonical_swaps Site swap types for canonical Monte Carlo events.
  ///     If size > 0, only these events are proposed.
  CanonicalEventGenerator(std::vector<monte::OccSwap> const &_canonical_swaps)
      : state(nullptr),
        occ_location(nullptr),
        canonical_swaps(_canonical_swaps) {
    if (canonical_swaps.size() == 0) {
      throw std::runtime_error(
          "Error in CanonicalEventGenerator: canonical_swaps.size() == 0");
    }
  }

  /// \brief The current state for which events are proposed and applied. Can be
  ///     nullptr, but must be set for use.
  state_type *state;

  /// Occupant tracker
  monte::OccLocation *occ_location;

  /// \brief Swap types for canonical Monte Carlo events
  std::vector<monte::OccSwap> canonical_swaps;

  /// \brief The current proposed event
  monte::OccEvent occ_event;

 public:
  /// \brief Set the current Monte Carlo state and occupant locations
  ///
  /// Notes:
  /// - Must be called before `propose` or `apply`
  ///
  /// \param _state The current state for which events are proposed and applied.
  ///     Throws if nullptr.
  /// \param _occ_location An occupant location tracker, which enables efficient
  ///     event proposal. It must already be initialized with the input state.
  ///     Throws if nullptr.
  void set(state_type *_state, monte::OccLocation *_occ_location) {
    this->state = throw_if_null(_state,
                                "Error in CanonicalEventGenerator::set: "
                                "_state==nullptr");
    this->occ_location = throw_if_null(_occ_location,
                                       "Error in CanonicalEventGenerator::set: "
                                       "_occ_location==nullptr");
  }

  /// \brief Propose a Monte Carlo occupation event, returning a reference
  ///
  /// Notes:
  /// - Must call `set` before `propose` or `apply`
  ///
  /// \param random_number_generator A random number generator
  monte::OccEvent const &propose(
      monte::RandomNumberGenerator<engine_type> &random_number_generator) {
    return monte::propose_canonical_event(this->occ_event, *this->occ_location,
                                          this->canonical_swaps,
                                          random_number_generator);
  }

  /// \brief Update the occupation of the current state using the provided event
  void apply(monte::OccEvent const &e) {
    this->occ_location->apply(e, get_occupation(*this->state));
  }
};

class CanonicalPotential : public BaseMontePotential {
 public:
  CanonicalPotential(std::shared_ptr<StateData> _state_data)
      : BaseMontePotential(_state_data),
        state(*state_data->state),
        n_unitcells(state_data->n_unitcells),
        occupation(get_occupation(state)),
        convert(*state_data->convert),
        composition_calculator(
            get_composition_calculator(*this->state_data->system)),
        composition_converter(
            get_composition_converter(*this->state_data->system)),
        param_composition(
            get_param_composition(*this->state_data->system, state.conditions)),
        formation_energy_clex(
            get_clex(*state_data->system, state, "formation_energy")) {
    if (param_composition.size() !=
        composition_converter.independent_compositions()) {
      throw std::runtime_error(
          "Error in CanonicalPotential: param_composition size error");
    }
  }

  // --- Data used in the potential calculation: ---

  state_type const &state;
  Index n_unitcells;
  Eigen::VectorXi const &occupation;
  monte::Conversions const &convert;
  composition::CompositionCalculator const &composition_calculator;
  composition::CompositionConverter const &composition_converter;
  Eigen::VectorXd param_composition;
  std::shared_ptr<clexulator::ClusterExpansion> formation_energy_clex;

  /// \brief Calculate (per_supercell) potential value
  double per_supercell() override {
    return formation_energy_clex->per_supercell();
  }

  /// \brief Calculate (per_unitcell) potential value
  double per_unitcell() override {
    return formation_energy_clex->per_unitcell();
  }

  /// \brief Calculate change in (per_supercell) potential value due
  ///     to a series of occupation changes
  double occ_delta_per_supercell(std::vector<Index> const &linear_site_index,
                                 std::vector<int> const &new_occ) override {
    return formation_energy_clex->occ_delta_value(linear_site_index, new_occ);
  }
};

class CanonicalCalculator : public BaseMonteCalculator {
 public:
  using BaseMonteCalculator::engine_type;

  CanonicalCalculator()
      : BaseMonteCalculator("CanonicalCalculator",  // calculator_name
                            {},                     // required_basis_set,
                            {},                     // required_local_basis_set,
                            {"formation_energy"},   // required_clex,
                            {},                     // required_multiclex,
                            {},                     // required_local_clex,
                            {},                     // required_local_multiclex,
                            {},                     // required_dof_spaces,
                            {},                     // required_params,
                            {},                     // optional_params,
                            false,                  // time_sampling_allowed,
                            false,                  // update_atoms,
                            false,                  // track_unique_atom_ids,
                            false,                  // save_atom_info,
                            false                   // is_multistate_method,
        ) {}

  /// \brief Construct functions that may be used to sample various quantities
  ///     of the Monte Carlo calculation as it runs
  std::map<std::string, state_sampling_function_type>
  standard_sampling_functions(
      std::shared_ptr<MonteCalculator> const &calculation) const override {
    std::vector<state_sampling_function_type> functions =
        monte_calculator::common_sampling_functions(
            calculation, "potential_energy",
            "Potential energy of the state (normalized per primitive cell)");

    // Specific to canonical
    // (none)

    std::map<std::string, state_sampling_function_type> function_map;
    for (auto const &f : functions) {
      function_map.emplace(f.name, f);
    }
    return function_map;
  }

  /// \brief Construct functions that may be used to sample various quantities
  ///     of the Monte Carlo calculation as it runs
  std::map<std::string, json_state_sampling_function_type>
  standard_json_sampling_functions(
      std::shared_ptr<MonteCalculator> const &calculation) const override {
    std::vector<json_state_sampling_function_type> functions =
        monte_calculator::common_json_sampling_functions(calculation);

    std::map<std::string, json_state_sampling_function_type> function_map;
    for (auto const &f : functions) {
      function_map.emplace(f.name, f);
    }
    return function_map;
  }

  /// \brief Construct functions that may be used to analyze Monte Carlo
  ///     calculation results
  std::map<std::string, results_analysis_function_type>
  standard_analysis_functions(
      std::shared_ptr<MonteCalculator> const &calculation) const override {
    std::vector<results_analysis_function_type> functions = {
        monte_calculator::make_heat_capacity_f(calculation)};

    std::map<std::string, results_analysis_function_type> function_map;
    for (auto const &f : functions) {
      function_map.emplace(f.name, f);
    }
    return function_map;
  }

  /// \brief Construct functions that may be used to modify states
  StateModifyingFunctionMap standard_modifying_functions(
      std::shared_ptr<MonteCalculator> const &calculation) const override {
    std::vector<StateModifyingFunction> functions = {
        monte_calculator::make_match_composition_f(calculation),
        monte_calculator::make_enforce_composition_f(calculation)};

    StateModifyingFunctionMap function_map;
    for (auto const &f : functions) {
      function_map.emplace(f.name, f);
    }
    return function_map;
  }

  /// \brief Construct functions that may be used to collect selected event data
  std::optional<monte::SelectedEventFunctions>
  standard_selected_event_functions(
      std::shared_ptr<MonteCalculator> const &calculation) const override {
    return std::nullopt;
  }

  /// \brief Construct default SamplingFixtureParams
  sampling_fixture_params_type make_default_sampling_fixture_params(
      std::shared_ptr<MonteCalculator> const &calculation, std::string label,
      bool write_results, bool write_trajectory, bool write_observations,
      bool write_status, std::optional<std::string> output_dir,
      std::optional<std::string> log_file,
      double log_frequency_in_s) const override {
    monte::SamplingParams sampling_params;
    {
      auto &s = sampling_params;
      s.sampler_names = {"clex.formation_energy", "potential_energy",
                         "mol_composition", "param_composition"};
      std::string prefix = "order_parameter.";
      for (auto const &pair : calculation->system()->dof_spaces) {
        s.sampler_names.push_back(prefix + pair.first);
      }
      prefix = "order_parameter.";
      std::string suffix = ".subspace_magnitudes";
      for (auto const &pair : calculation->system()->dof_subspaces) {
        s.sampler_names.push_back(prefix + pair.first + suffix);
      }
      if (write_trajectory) {
        s.do_sample_trajectory = true;
      }
    }

    monte::CompletionCheckParams<statistics_type> completion_check_params;
    {
      auto &c = completion_check_params;
      c.equilibration_check_f = monte::default_equilibration_check;
      c.calc_statistics_f =
          monte::default_statistics_calculator<statistics_type>();

      converge(calculation->sampling_functions, completion_check_params)
          .set_abs_precision("potential_energy", 0.001);
    }

    std::vector<std::string> analysis_names = {"heat_capacity"};

    return clexmonte::make_sampling_fixture_params(
        label, calculation->sampling_functions,
        calculation->json_sampling_functions, calculation->analysis_functions,
        sampling_params, completion_check_params, analysis_names, write_results,
        write_trajectory, write_observations, write_status, output_dir,
        log_file, log_frequency_in_s);
  }

  /// \brief Validate the state's configuration
  ///
  /// Notes:
  /// - All configurations are valid (validate_state checks for consistency
  ///   with the composition conditions)
  Validator validate_configuration(state_type &state) const override {
    return Validator{};
  }

  /// \brief Validate state's conditions
  ///
  /// Notes:
  /// - requires scalar temperature
  /// - validate composition consistency
  /// - warnings if other conditions are present
  Validator validate_conditions(state_type &state) const override {
    // Validate system
    if (this->system == nullptr) {
      throw std::runtime_error(
          "Error in CanonicalCalculator::validate_conditions: system==nullptr");
    }

    // validate state.conditions
    monte::ValueMap const &conditions = state.conditions;
    Validator v;
    v.insert(validate_keys(conditions.scalar_values,
                           {"temperature"} /*required*/, {} /*optional*/,
                           "scalar", "condition", false /*throw_if_invalid*/));
    v.insert(
        validate_keys(conditions.vector_values, {} /*required*/,
                      {"param_composition", "mol_composition"} /*optional*/,
                      "vector", "condition", false /*throw_if_invalid*/));
    v.insert(validate_composition_consistency(
        state, get_composition_converter(*this->system),
        this->mol_composition_tol));
    return v;
  }

  /// \brief Validate state
  Validator validate_state(state_type &state) const override {
    Validator v;
    v.insert(this->validate_configuration(state));
    v.insert(this->validate_conditions(state));
    if (!v.valid()) {
      return v;
    }

    // check if configuration is consistent with conditions
    auto const &composition_calculator =
        get_composition_calculator(*this->system);
    auto const &composition_converter =
        get_composition_converter(*this->system);

    Eigen::VectorXd mol_composition =
        composition_calculator.mean_num_each_component(get_occupation(state));
    Eigen::VectorXd param_composition =
        composition_converter.param_composition(mol_composition);

    Eigen::VectorXd target_mol_composition =
        get_mol_composition(*this->system, state.conditions);
    Eigen::VectorXd target_param_composition =
        composition_converter.param_composition(target_mol_composition);

    if (!CASM::almost_equal(mol_composition, target_mol_composition,
                            this->mol_composition_tol)) {
      std::stringstream msg;
      msg << "***" << std::endl;
      msg << "Calculated composition is not consistent with conditions "
             "composition."
          << std::endl;
      msg << "Calculated composition:" << std::endl;
      msg << "- mol_composition: " << mol_composition.transpose() << std::endl;
      msg << "- param_composition: " << param_composition.transpose()
          << std::endl;
      msg << "Conditions:" << std::endl;
      msg << "- mol_composition: " << target_mol_composition.transpose()
          << std::endl;
      msg << "- param_composition: " << target_param_composition.transpose()
          << std::endl;
      msg << "***" << std::endl;
      v.error.insert(msg.str());
    }
    return v;
  }

  /// \brief Validate and set the current state, construct state_data, construct
  ///     potential
  ///
  /// \param state Conditions is required to have `mol_composition` or
  ///     `param_composition`. If both are present, they are checked for
  ///     consistency. If not consistent, a warning is printed and
  ///     param_composition is modified to match mol_composition.
  ///     If both are not present, one is set from the other.
  /// \param occ_location Optional occupation location tracking. Not required
  ///     for potential evaluation. Required for a Monte Carlo run.
  void set_state_and_potential(state_type &state,
                               monte::OccLocation *occ_location) override {
    // Validate system
    if (this->system == nullptr) {
      throw std::runtime_error(
          "Error in CanonicalCalculator::run: system==nullptr");
    }

    // Validate state
    Validator v = this->validate_state(state);
    print(CASM::log(), v);
    if (!v.valid()) {
      throw std::runtime_error(
          "Error in CanonicalCalculator::run: Invalid initial state");
    }

    // Make state data
    this->state_data =
        std::make_shared<StateData>(this->system, &state, occ_location);

    // Make potential calculator
    this->potential = std::make_shared<CanonicalPotential>(this->state_data);
  }

  /// \brief Set event data (includes calculating all rates), using current
  /// state data
  void set_event_data() override {
    throw std::runtime_error(
        "Error in CanonicalCalculator::set_event_data: not valid");
  }

  /// \brief Perform a single run, evolving current state
  ///
  /// Notes:
  /// - state and occ_location are evolved and end in modified states
  /// - if run_manager.engine != nullptr, set this->engine = run_manager.engine
  void run(state_type &state, monte::OccLocation &occ_location,
           run_manager_type<engine_type> &run_manager) override {
    // Set state data and construct potential calculator
    this->set_state_and_potential(state, &occ_location);

    // Random number generator
    if (run_manager.engine == nullptr) {
      throw std::runtime_error(
          "Error in CanonicalCalculator::run: run_manager.engine==nullptr");
    }
    this->engine = run_manager.engine;
    monte::RandomNumberGenerator<engine_type> random_number_generator(
        run_manager.engine);

    //    // Enforce composition
    //    this->enforce_composition(state, occ_location,
    //    random_number_generator);

    // Get temperature
    double temperature = state.conditions.scalar_values.at("temperature");

    // Make delta potential function
    auto potential_occ_delta_per_supercell_f =
        [=](monte::OccEvent const &event) {
          return this->potential->occ_delta_per_supercell(
              event.linear_site_index, event.new_occ);
        };

    // Make event generator
    auto event_generator = std::make_shared<CanonicalEventGenerator>(
        get_canonical_swaps(*this->system));
    event_generator->set(&state, &occ_location);

    // Make event proposal function
    auto propose_event_f =
        [=](monte::RandomNumberGenerator<engine_type> &random_number_generator)
        -> monte::OccEvent const & {
      return event_generator->propose(random_number_generator);
    };

    // Make event application function
    auto apply_event_f = [=](monte::OccEvent const &occ_event) -> void {
      return event_generator->apply(occ_event);
    };

    // Run Monte Carlo at a single condition
    clexmonte::occupation_metropolis_v2(
        state, occ_location, temperature, potential_occ_delta_per_supercell_f,
        propose_event_f, apply_event_f, run_manager);
  }

  /// \brief Perform a single run, evolving one or more states
  void run(int current_state, std::vector<state_type> &states,
           std::vector<monte::OccLocation> &occ_locations,
           run_manager_type<engine_type> &run_manager) override {
    throw std::runtime_error(
        "Error: CanonicalCalculator does not allow multi-state runs");
  }

  // --- Parameters ---
  int verbosity_level = 10;
  double mol_composition_tol = CASM::TOL;

  /// \brief Reset the derived Monte Carlo calculator
  ///
  /// Parameters:
  ///   verbosity: str or int, default=10
  ///       If integer, the allowed range is `[0,100]`. If string, then:
  ///       - "none" is equivalent to integer value 0
  ///       - "quiet" is equivalent to integer value 5
  ///       - "standard" is equivalent to integer value 10
  ///       - "verbose" is equivalent to integer value 20
  ///       - "debug" is equivalent to integer value 100
  void _reset() override {
    ParentInputParser parser{params};

    // "verbosity": str or int, default=10
    this->verbosity_level = parse_verbosity(parser);
    CASM::log().set_verbosity(this->verbosity_level);

    // "mol_composition_tol": float, default=CASM::TOL
    this->mol_composition_tol = CASM::TOL;
    parser.optional(this->mol_composition_tol, "mol_composition_tol");

    // TODO: enumeration

    std::stringstream ss;
    ss << "Error in CanonicalCalculator: error reading calculation "
          "parameters.";
    std::runtime_error error_if_invalid{ss.str()};
    report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);

    return;
  }

  /// \brief Clone the CanonicalCalculator
  CanonicalCalculator *_clone() const override {
    return new CanonicalCalculator(*this);
  }

  /// \brief Enforce composition conditions consistency
  ///
  /// - If both present and not consistent, set param_composition to be
  ///   consistent with mol_composition and print warning
  /// - If only one set, set the other to be consistent
  void enforce_conditions_consistency(state_type &state) {
    enforce_composition_consistency(state,
                                    get_composition_converter(*this->system),
                                    this->mol_composition_tol);
  }

  /// \brief Enforce composition conditions
  void enforce_composition(
      state_type &state, monte::OccLocation &occ_location,
      monte::RandomNumberGenerator<engine_type> &random_number_generator) {
    auto const &composition_calculator =
        get_composition_calculator(*this->system);
    auto const &composition_converter =
        get_composition_converter(*this->system);

    Eigen::VectorXd mol_composition =
        composition_calculator.mean_num_each_component(get_occupation(state));
    Eigen::VectorXd param_composition =
        composition_converter.param_composition(mol_composition);

    Eigen::VectorXd target_mol_composition =
        get_mol_composition(*this->system, state.conditions);
    Eigen::VectorXd target_param_composition =
        composition_converter.param_composition(target_mol_composition);

    auto &log = CASM::log();
    log.begin<Log::quiet>("Enforcing mol_composition conditions");
    log.indent() << "calculated:" << std::endl;
    log.indent() << "- mol_composition: " << mol_composition.transpose()
                 << std::endl;
    log.indent() << "- param_composition: " << param_composition.transpose()
                 << std::endl;
    log.indent() << "target: " << std::endl;
    log.indent() << "- mol_composition: " << target_mol_composition.transpose()
                 << std::endl;
    log.indent() << "- param_composition: "
                 << target_param_composition.transpose() << std::endl;

    log.indent() << "enforcing composition..." << std::endl;
    clexmonte::enforce_composition(get_occupation(state),
                                   target_mol_composition,
                                   get_composition_calculator(*this->system),
                                   get_semigrand_canonical_swaps(*this->system),
                                   occ_location, random_number_generator);
    log.indent() << "DONE" << std::endl;

    mol_composition =
        composition_calculator.mean_num_each_component(get_occupation(state));
    param_composition =
        composition_converter.param_composition(mol_composition);
    log.indent() << "calculated:" << std::endl;
    log.indent() << "- mol_composition: " << mol_composition.transpose()
                 << std::endl;
    log.indent() << "- param_composition: " << param_composition.transpose()
                 << std::endl
                 << std::endl;
  }
};

}  // namespace clexmonte
}  // namespace CASM

extern "C" {
/// \brief Returns a clexmonte::BaseMonteCalculator* owning a
/// CanonicalCalculator
CASM::clexmonte::BaseMonteCalculator *make_CanonicalCalculator() {
  return new CASM::clexmonte::CanonicalCalculator();
}
}
