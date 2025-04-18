#include "casm/clexmonte/monte_calculator/KineticCalculator.hh"

#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/events/io/json/event_data_json_io.hh"
#include "casm/clexmonte/misc/check_params.hh"
#include "casm/clexmonte/misc/to_json.hh"
#include "casm/clexmonte/monte_calculator/MonteEventData.hh"
#include "casm/clexmonte/monte_calculator/analysis_functions.hh"
#include "casm/clexmonte/monte_calculator/kinetic_events.hh"
#include "casm/clexmonte/monte_calculator/kinetic_sampling_functions.hh"
#include "casm/clexmonte/monte_calculator/modifying_functions.hh"
#include "casm/clexmonte/monte_calculator/sampling_functions.hh"
#include "casm/clexmonte/monte_calculator/selected_event_functions.hh"
#include "casm/clexmonte/run/functions.hh"
#include "casm/configuration/io/json/Configuration_json_io.hh"
#include "casm/misc/TypeInfo.hh"
#include "casm/monte/events/OccEventProposal.hh"
#include "casm/monte/run_management/RunManager.hh"
#include "casm/monte/sampling/RequestedPrecisionConstructor.hh"
#include "casm/monte/sampling/io/json/SelectedEventFunctions_json_io.hh"

namespace CASM {
namespace clexmonte {
namespace kinetic_2 {

namespace {

template <int verbosity_level, typename T, typename U>
void read_option_with_verbosity(ParentInputParser &parser, Log &log, T &param,
                                fs::path option, U const &default_value) {
  param = default_value;
  parser.optional(param, option);

  log.begin_section<verbosity_level>();
  log.indent() << option.string() << "=" << std::boolalpha << param
               << std::endl;
  log.end_section();
}

template <typename T, typename U>
void read_option(ParentInputParser &parser, Log &log, T &param, fs::path option,
                 U const &default_value) {
  read_option_with_verbosity<Log::standard>(parser, log, param, option,
                                            default_value);
}

}  // namespace

KineticPotential::KineticPotential(std::shared_ptr<StateData> _state_data)
    : BaseMontePotential(_state_data),
      state(*state_data->state),
      formation_energy_clex(
          get_clex(*state_data->system, state, "formation_energy")) {}

/// \brief Calculate (per_supercell) potential value
double KineticPotential::per_supercell() {
  return formation_energy_clex->per_supercell();
}

/// \brief Calculate (per_unitcell) potential value
double KineticPotential::per_unitcell() {
  return formation_energy_clex->per_unitcell();
}

/// \brief Calculate change in (per_supercell) potential value due
///     to a series of occupation changes
double KineticPotential::occ_delta_per_supercell(
    std::vector<Index> const &linear_site_index,
    std::vector<int> const &new_occ) {
  return formation_energy_clex->occ_delta_value(linear_site_index, new_occ);
}

KineticCalculator::KineticCalculator()
    : BaseMonteCalculator(
          "KineticCalculator",   // calculator_name
          {},                    // required_basis_set,
          {},                    // required_local_basis_set,
          {"formation_energy"},  // required_clex,
          {},                    // required_multiclex,
          {},                    // required_local_clex,
          {},                    // required_local_multiclex,
          {},                    // required_dof_spaces,
          {},                    // required_params,
          {"verbosity", "print_event_data_summary", "mol_composition_tol",
           "event_data_type", "event_selector_type", "abnormal_event_handling",
           "impact_table_type", "assign_allowed_events_only",
           "selected_event_data", "state_graph"},  // optional_params,
          true,                                    // time_sampling_allowed,
          true,                                    // update_atoms,
          false,                                   // save_atom_info,
          false                                    // is_multistate_method,
      ) {
  // this could go into base constructor
  this->selected_event = std::make_shared<SelectedEvent>();
  this->selected_event_functions =
      std::make_shared<monte::SelectedEventFunctions>();
  this->selected_event_data = std::make_shared<monte::SelectedEventData>();
}

/// \brief Construct functions that may be used to sample various quantities
///     of the Monte Carlo calculation as it runs
std::map<std::string, state_sampling_function_type>
KineticCalculator::standard_sampling_functions(
    std::shared_ptr<MonteCalculator> const &calculation) const {
  using namespace monte_calculator;

  std::vector<state_sampling_function_type> functions =
      common_sampling_functions(
          calculation, "potential_energy",
          "Potential energy of the state (normalized per primitive cell)");

  // Specific to kmc
  functions.push_back(make_mean_R_squared_collective_isotropic_f(calculation));
  functions.push_back(
      make_mean_R_squared_collective_anisotropic_f(calculation));
  functions.push_back(make_mean_R_squared_individual_isotropic_f(calculation));
  functions.push_back(
      make_mean_R_squared_individual_anisotropic_f(calculation));
  functions.push_back(make_L_isotropic_f(calculation));
  functions.push_back(make_L_anisotropic_f(calculation));
  functions.push_back(make_D_tracer_isotropic_f(calculation));
  functions.push_back(make_D_tracer_anisotropic_f(calculation));
  functions.push_back(make_jumps_per_atom_by_type_f(calculation));
  functions.push_back(make_jumps_per_event_by_type_f(calculation));
  functions.push_back(make_jumps_per_atom_per_event_by_type_f(calculation));
  functions.push_back(make_selected_event_count_by_type_f(calculation));
  functions.push_back(make_selected_event_fraction_by_type_f(calculation));
  functions.push_back(
      make_selected_event_count_by_equivalent_index_f(calculation));
  functions.push_back(
      make_selected_event_fraction_by_equivalent_index_f(calculation));
  functions.push_back(
      make_selected_event_count_by_equivalent_index_and_direction_f(
          calculation));
  functions.push_back(
      make_selected_event_fraction_by_equivalent_index_and_direction_f(
          calculation));

  for (auto f : make_selected_event_count_by_equivalent_index_per_event_type_f(
           calculation)) {
    functions.push_back(f);
  }
  for (auto f :
       make_selected_event_fraction_by_equivalent_index_per_event_type_f(
           calculation)) {
    functions.push_back(f);
  }

  std::map<std::string, state_sampling_function_type> function_map;
  for (auto const &f : functions) {
    function_map.emplace(f.name, f);
  }
  return function_map;
}

/// \brief Construct functions that may be used to sample various quantities
///     of the Monte Carlo calculation as it runs
std::map<std::string, json_state_sampling_function_type>
KineticCalculator::standard_json_sampling_functions(
    std::shared_ptr<MonteCalculator> const &calculation) const {
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
KineticCalculator::standard_analysis_functions(
    std::shared_ptr<MonteCalculator> const &calculation) const {
  std::vector<results_analysis_function_type> functions = {
      monte_calculator::make_heat_capacity_f(calculation)};

  std::map<std::string, results_analysis_function_type> function_map;
  for (auto const &f : functions) {
    function_map.emplace(f.name, f);
  }
  return function_map;
}

/// \brief Construct functions that may be used to modify states
StateModifyingFunctionMap KineticCalculator::standard_modifying_functions(
    std::shared_ptr<MonteCalculator> const &calculation) const {
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
KineticCalculator::standard_selected_event_functions(
    std::shared_ptr<MonteCalculator> const &calculation) const {
  using namespace monte_calculator;
  monte::SelectedEventFunctions functions;

  // Event type data:
  functions.insert(make_selected_event_by_type_f(calculation));
  functions.insert(make_selected_event_by_equivalent_index_f(calculation));
  functions.insert(
      make_selected_event_by_equivalent_index_and_direction_f(calculation));
  for (auto f :
       make_selected_event_by_equivalent_index_per_event_type_f(calculation)) {
    functions.insert(f);
  }
  for (auto f : make_local_orbit_composition_f(calculation)) {
    functions.insert(f);
  }

  // Event state data:
  functions.insert(make_dE_activated_by_type_f(calculation));
  functions.insert(make_dE_activated_by_equivalent_index_f(calculation));

  return functions;
}

/// \brief Construct default SamplingFixtureParams
sampling_fixture_params_type
KineticCalculator::make_default_sampling_fixture_params(
    std::shared_ptr<MonteCalculator> const &calculation, std::string label,
    bool write_results, bool write_trajectory, bool write_observations,
    bool write_status, std::optional<std::string> output_dir,
    std::optional<std::string> log_file, double log_frequency_in_s) const {
  monte::SamplingParams sampling_params;
  {
    auto &s = sampling_params;
    s.sampler_names = {"clex.formation_energy",
                       "potential_energy",
                       "mol_composition",
                       "param_composition",
                       "mean_R_squared_collective_isotropic",
                       "mean_R_squared_individual_isotropic",
                       "L_isotropic",
                       "D_tracer_isotropic",
                       "mean_R_squared_collective_anisotropic",
                       "mean_R_squared_individual_anisotropic",
                       "L_anisotropic",
                       "D_tracer_anisotropic",
                       "jumps_per_atom_by_type",
                       "jumps_per_event_by_type",
                       "jumps_per_atom_per_event_by_type"};
    s.do_sample_time = true;
    std::string prefix;
    prefix = "order_parameter.";
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

    completion_check_params.cutoff_params.max_count = 100;
  }

  std::vector<std::string> analysis_names = {"heat_capacity"};

  return clexmonte::make_sampling_fixture_params(
      label, calculation->sampling_functions,
      calculation->json_sampling_functions, calculation->analysis_functions,
      sampling_params, completion_check_params, analysis_names, write_results,
      write_trajectory, write_observations, write_status, output_dir, log_file,
      log_frequency_in_s);
}

/// \brief Validate the state's configuration
///
/// Notes:
/// - All configurations are valid (validate_state checks for consistency
///   with the composition conditions)
Validator KineticCalculator::validate_configuration(state_type &state) const {
  return Validator{};
}

/// \brief Validate state's conditions
///
/// Notes:
/// - requires scalar temperature
/// - validate composition consistency
/// - warnings if other conditions are present
Validator KineticCalculator::validate_conditions(state_type &state) const {
  // Validate system
  if (this->system == nullptr) {
    throw std::runtime_error(
        "Error in KineticCalculator::validate_conditions: system==nullptr");
  }

  // validate state.conditions
  monte::ValueMap const &conditions = state.conditions;
  Validator v;
  v.insert(validate_keys(conditions.scalar_values, {"temperature"} /*required*/,
                         {} /*optional*/, "scalar", "condition",
                         false /*throw_if_invalid*/));
  v.insert(validate_keys(conditions.vector_values, {} /*required*/,
                         {"param_composition", "mol_composition"} /*optional*/,
                         "vector", "condition", false /*throw_if_invalid*/));
  v.insert(validate_composition_consistency(
      state, get_composition_converter(*this->system),
      this->mol_composition_tol));
  return v;
}

/// \brief Validate state
Validator KineticCalculator::validate_state(state_type &state) const {
  Validator v;
  v.insert(this->validate_configuration(state));
  v.insert(this->validate_conditions(state));
  if (!v.valid()) {
    return v;
  }

  // check if configuration is consistent with conditions
  auto const &composition_calculator =
      get_composition_calculator(*this->system);
  auto const &composition_converter = get_composition_converter(*this->system);

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
/// \param state State to set
/// \param occ_location Pointer to OccLocation to use, or may be nullptr
void KineticCalculator::set_state_and_potential(
    state_type &state, monte::OccLocation *occ_location) {
  // Validate system
  if (this->system == nullptr) {
    throw std::runtime_error(
        "Error in KineticCalculator::run: system==nullptr");
  }

  // Validate state
  Validator v = this->validate_state(state);
  clexmonte::print(CASM::log(), v);
  if (!v.valid()) {
    throw std::runtime_error(
        "Error in KineticCalculator::run: Invalid initial state");
  }

  // Make state data
  this->state_data =
      std::make_shared<StateData>(this->system, &state, occ_location);

  // Make potential calculator
  this->potential = std::make_shared<KineticPotential>(this->state_data);
}

/// \brief Set event data (includes calculating all rates), using current
/// state data
///
/// Notes:
/// - Validates this->state_data is not null
/// - Validates this->state_data->occ_location is not null
/// - Uses this->engine for the event selector
/// - Resets the `n_encountered_abnormal` and `n_selected_abnormal`
///   counters.
void KineticCalculator::set_event_data() {
  if (this->state_data == nullptr) {
    throw std::runtime_error(
        "Error in KineticCalculator::set_event_data: "
        "this->state_data==nullptr");
  }
  if (this->state_data->occ_location == nullptr) {
    throw std::runtime_error(
        "Error in KineticCalculator::set_event_data: "
        "this->state_data->occ_location==nullptr");
  }
  if (this->event_data == nullptr) {
    throw std::runtime_error(
        "Error in KineticCalculator::set_event_data: "
        "this->event_data==nullptr");
  }

  // Currently, event_filters are only set at _reset() by reading from params
  this->event_data->update(this->state_data, event_filters, this->engine);
}

/// \brief Perform a single run, evolving current state
void KineticCalculator::run(state_type &state, monte::OccLocation &occ_location,
                            run_manager_type<engine_type> &run_manager) {
  Log &log = CASM::log();

  if (run_manager.engine == nullptr) {
    throw std::runtime_error(
        "Error in KineticCalculator::run: run_manager.engine==nullptr");
  }
  this->engine = run_manager.engine;

  if (run_manager.sampling_fixtures.size() == 0) {
    throw std::runtime_error(
        "Error in KineticCalculator::run: "
        "run_manager.sampling_fixtures.size()==0");
  }

  // Set state and potential
  // - Validates this->system is not null
  // - Validates state
  // - Throw if validation fails
  // - Constructs this->state_data
  // - Constructs this->potential
  this->set_state_and_potential(state, &occ_location);

  // Set event data
  // - Set or re-set state in event_state_calculators
  // - Update event_data if supercell has changed
  // - Throw if this->state_data is null
  // - Constructs this->event_data->event_selector
  // - Calculates all rates
  log.begin_section<Log::standard>();
  log.indent() << "Setting event data ... " << std::endl;
  this->set_event_data();
  log.indent() << "Setting event data ... DONE" << std::endl << std::endl;

  // Construct EventDataSummary
  if (this->print_event_data_summary) {
    log.indent() << "Generating event data summary ... " << std::endl;
    MonteEventData monte_event_data(this->event_data, nullptr);
    double energy_bin_width = 0.1;
    double freq_bin_width = 0.1;
    double rate_bin_width = 0.1;

    EventDataSummary event_data_summary(this->state_data, monte_event_data,
                                        energy_bin_width, freq_bin_width,
                                        rate_bin_width);
    log.indent() << "Generating event data summary ... DONE" << std::endl
                 << std::endl;
    log.end_section();
    print<Log::standard>(log, event_data_summary);

    if (event_data_summary.n_events_allowed == 0) {
      throw std::runtime_error("Error: Cannot run. No allowed events.");
    }
  }

  //  if (event_data_summary.n_events_allowed == 0) {
  //    throw std::runtime_error("Error: Cannot run. No allowed events.");
  //  }

  // Construct KMCData
  this->kmc_data = std::make_shared<kmc_data_type>();

  // Update atom_name_index_list
  // -- These do not change (no atoms moving to/from reservoirs) --
  auto event_system = get_event_system(*this->system);

  // Optional: Manages constructing histogram data structures
  // and collecting selected event data

  std::optional<monte::SelectedEventDataCollector> collector;
  if (this->selected_event_function_params) {
    if (this->selected_event_functions == nullptr) {
      throw std::runtime_error(
          "Error in KineticCalculator::run: "
          "this->selected_event_functions==nullptr");
    }
    collector = monte::SelectedEventDataCollector(
        *selected_event_functions, *selected_event_function_params,
        selected_event_data);
  }

  // Check this->selected_event is not null
  if (this->selected_event == nullptr) {
    throw std::runtime_error(
        "Error in KineticCalculator::run: this->selected_event==nullptr");
  }

  this->event_data->run(state, occ_location, *this->selected_event, collector,
                        run_manager, this->kmc_data, event_system);

  //  // Function to set selected event
  //  bool requires_event_state =
  //      collector.has_value() && collector->requires_event_state;
  //  auto set_selected_event_f = [=](SelectedEvent &selected_event) {
  //    this->event_data->select_event(selected_event, requires_event_state);
  //  };
  //
  //  // Run Kinetic Monte Carlo at a single condition
  //  kinetic_monte_carlo_v2<EventID>(state, occ_location, *this->kmc_data,
  //                                  *this->selected_event,
  //                                  set_selected_event_f, collector,
  //                                  run_manager, event_system);

  // Warn if abnormal events were encountered or selected
  check_n_encountered_abnormal(this->event_data->n_encountered_abnormal);
  check_n_selected_abnormal(this->event_data->n_selected_abnormal);
}

/// \brief Print a warning to std::cerr if abnormal events were
///     encountered
void KineticCalculator::check_n_encountered_abnormal(
    std::map<std::string, Index> const &n_encountered_abnormal) const {
  if (n_encountered_abnormal.empty()) {
    return;
  }
  Log &log = CASM::err_log();
  log << "## WARNING: ENCOUNTERED ABNORMAL EVENTS #############\n"
         "#                                                   #\n"
         "# Number encountered by type:                       #\n";
  for (auto const &pair : n_encountered_abnormal) {
    log << "  - " << pair.first << ": " << pair.second << "\n";
  }
  log << "#                                                   #\n"
         "#####################################################\n"
      << std::endl;
}

/// \brief Print a warning to std::cerr if abnormal events were
///     encountered
void KineticCalculator::check_n_selected_abnormal(
    std::map<std::string, Index> const &n_selected_abnormal) const {
  if (n_selected_abnormal.empty()) {
    return;
  }

  Log &log = CASM::err_log();
  log << "## WARNING: SELECTED ABNORMAL EVENTS ################\n"
         "#                                                   #\n"
         "# Number selected by type:                          #\n";
  for (auto const &pair : n_selected_abnormal) {
    log << "  - " << pair.first << ": " << pair.second << "\n";
  }
  log << "#                                                   #\n"
         "#####################################################\n"
      << std::endl;
}

/// \brief Perform a single run, evolving one or more states
void KineticCalculator::run(int current_state, std::vector<state_type> &states,
                            std::vector<monte::OccLocation> &occ_locations,
                            run_manager_type<engine_type> &run_manager) {
  throw std::runtime_error(
      "Error: KineticCalculator does not allow multi-state runs");
}

template <bool DebugMode>
void KineticCalculator::make_complete_event_data_impl() {
  if constexpr (DebugMode) {
    auto &log = CASM::log();
    log << "!! make_complete_event_data_impl !! " << std::endl;
    log << "!! DebugMode=" << DebugMode << " !! " << std::endl;
    log << std::endl;
  }

  this->event_data = std::make_shared<CompleteKineticEventData<DebugMode>>(
      system, event_filters, event_data_options);
}

template <bool DebugMode>
void KineticCalculator::make_allowed_event_data_impl() {
  if constexpr (DebugMode) {
    auto &log = CASM::log();
    log << "!! make_allowed_event_data_impl !! " << std::endl;
    log << "!! DebugMode=" << DebugMode << " !! " << std::endl;
    log << std::endl;
  }

  typedef AllowedEventCalculator<DebugMode> event_calculator_type;

  if (this->event_selector_type ==
      kinetic_event_selector_type::vector_sum_tree) {
    typedef vector_sum_tree_event_selector_type<event_calculator_type>
        event_selector_type;

    this->event_data = std::make_shared<
        AllowedKineticEventData<event_selector_type, DebugMode>>(
        system, event_data_options);

  } else if (this->event_selector_type ==
             kinetic_event_selector_type::sum_tree) {
    typedef sum_tree_event_selector_type<event_calculator_type>
        event_selector_type;

    this->event_data = std::make_shared<
        AllowedKineticEventData<event_selector_type, DebugMode>>(
        system, event_data_options);

  } else if (this->event_selector_type ==
             kinetic_event_selector_type::direct_sum) {
    typedef direct_sum_event_selector_type<event_calculator_type>
        event_selector_type;

    this->event_data = std::make_shared<
        AllowedKineticEventData<event_selector_type, DebugMode>>(
        system, event_data_options);

  } else {
    throw std::runtime_error(
        "Error in KineticCalculator: "
        "invalid event_selector_type for event_data_type");
  }
}

/// \brief Reset the derived Monte Carlo calculator
///
/// Parameters:
///
///   verbosity: str or int, default=10
///       If integer, the allowed range is `[0,100]`. If string, then:
///       - "none" is equivalent to integer value 0
///       - "quiet" is equivalent to integer value 5
///       - "standard" is equivalent to integer value 10
///       - "verbose" is equivalent to integer value 20
///       - "debug" is equivalent to integer value 100
void KineticCalculator::_reset() {
  // -- Parsing ----------------------------

  // `params` is BaseMonteCalculator::params
  // - originally "calculation_options" from the run params input file
  // - `params` argument of MonteCalculator Python constructor

  ParentInputParser parser{params};

  // "verbosity": str or int, default=10
  this->verbosity_level = parse_verbosity(parser);
  CASM::log().set_verbosity(this->verbosity_level);

  auto &log = CASM::log();

  // system check:
  {
    auto it = this->system->clex_data.find("entropy");
    bool found_entropy = (it != this->system->clex_data.end());

    log.read<Log::standard>("System parameters");
    if (!found_entropy) {
      log.verbatim(
          "No \"entropy\" clex found:                                        \n"
          "- Event rates will be calculated without a change in entropy.     \n"
          "- A local cluster expansion will be used for the attempt          \n"
          "  frequency.\n");
    } else {
      log.verbatim(
          "Found \"entropy\" clex:                                           \n"
          "- Event rates will be calculated using the change in entropy.     \n"
          "- A local cluster expansion will be used for the kinetically-     \n"
          "  resolved activation entropy.\n");
    }
    log << std::endl;
  }

  // read other parameters

  log.read<Log::standard>("KineticCalculator parameters");
  log.indent() << "verbosity=" << this->verbosity_level << std::endl;

  // "print_event_data_summary": bool, default=false
  this->print_event_data_summary = false;
  parser.optional(this->print_event_data_summary, "print_event_data_summary");
  log.indent() << "print_event_data_summary=" << std::boolalpha
               << this->print_event_data_summary << std::endl;

  // "mol_composition_tol": float, default=CASM::TOL
  this->mol_composition_tol = CASM::TOL;
  parser.optional(this->mol_composition_tol, "mol_composition_tol");
  log.indent() << "mol_composition_tol=" << this->mol_composition_tol
               << std::endl;

  // TODO: enumeration

  // TODO: Read event_filters from params
  this->event_filters = std::nullopt;
  log.indent() << "event_filters=" << qto_json(this->event_filters)
               << std::endl;

  // Read selected event data params
  this->selected_event_function_params.reset();
  if (parser.self.contains("selected_event_data")) {
    auto selected_event_data_subparser =
        parser.subparse<monte::SelectedEventFunctionParams>(
            "selected_event_data");
    if (selected_event_data_subparser->valid()) {
      this->selected_event_function_params =
          std::move(selected_event_data_subparser->value);
      log.indent() << "selected_event_data="
                   << qto_json(*this->selected_event_function_params)
                   << std::endl;
    } else {
      log.indent() << "selected_event_data=invalid" << std::endl;
    }
  } else {
    log.indent() << "selected_event_data=null" << std::endl;
  }

  // Read "event_data_type"
  // - "high_memory": complete event list,
  // - "default" (default): allowed event list w/ vector index, or
  // - "low_memory": allowed event list w/ map index
  std::string event_data_type_str = "default";
  parser.optional(event_data_type_str, "event_data_type");
  if (event_data_type_str == "high_memory") {
    this->event_data_type = kinetic_event_data_type::high_memory;
    log.indent() << "event_data_type="
                 << "\"high_memory\"" << std::endl;
  } else if (event_data_type_str == "default") {
    this->event_data_type = kinetic_event_data_type::default_memory;
    log.indent() << "event_data_type="
                 << "\"default\"" << std::endl;

    if (event_filters.has_value()) {
      parser.insert_error(
          "event_data_type",
          "event_filters are not supported by event_data_type 'default'");
    }

  } else {
    parser.insert_error("event_data_type",
                        "Invalid event_data_type: " + event_data_type_str);
  }

  // Read "event_selector_type"
  // - "vector_sum_tree" (default): binary sum tree (log complexity) built
  //   using std::vector
  // - "sum_tree": binary sum tree (log complexity) built using std::map and
  //   linked lists
  // - "direct_sum": direct sum of rate vector (linear complexity)
  std::string event_selector_type_str = "vector_sum_tree";
  parser.optional(event_selector_type_str, "event_selector_type");
  if (event_selector_type_str == "vector_sum_tree") {
    this->event_selector_type = kinetic_event_selector_type::vector_sum_tree;
    log.indent() << "event_selector_type="
                 << "\"vector_sum_tree\"" << std::endl;
  } else if (event_selector_type_str == "sum_tree") {
    this->event_selector_type = kinetic_event_selector_type::sum_tree;
    log.indent() << "event_selector_type="
                 << "\"sum_tree\"" << std::endl;
  } else if (event_selector_type_str == "direct_sum") {
    this->event_selector_type = kinetic_event_selector_type::direct_sum;
    log.indent() << "event_selector_type="
                 << "\"direct_sum\"" << std::endl;
  } else {
    parser.insert_error("event_selector_type", "Invalid event_selector_type: " +
                                                   event_selector_type_str);
  }

  // -- Abnormal event handling
  fs::path base("abnormal_event_handling");
  check_params(params[base], {} /*required_params*/,
               {"output_dir", "tol", "encountered_events",
                "selected_events"} /*optional_params*/,
               base);
  read_option(parser, log, this->event_data_options.output_dir,
              base / "output_dir", fs::path("output"));

  read_option(parser, log, this->event_data_options.local_corr_compare_tol,
              base / "tol", CASM::TOL);

  // -- Encountered events
  check_params(params[base]["encountered_events"], {} /*required_params*/,
               {"warn", "throw", "n_write", "disallow"} /*optional_params*/,
               base / "encountered_events");
  read_option(parser, log,
              this->event_data_options.warn_if_encountered_event_is_abnormal,
              base / "encountered_events" / "warn", true);

  read_option(parser, log,
              this->event_data_options.throw_if_encountered_event_is_abnormal,
              base / "encountered_events" / "throw", true);

  read_option(parser, log,
              this->event_data_options.n_write_if_encountered_event_is_abnormal,
              base / "encountered_events" / "n_write", 100);

  read_option(
      parser, log,
      this->event_data_options.disallow_if_encountered_event_is_abnormal,
      base / "encountered_events" / "disallow", false);

  // -- Selected events
  check_params(params[base]["selected_events"], {} /*required_params*/,
               {"warn", "throw", "n_write", "disallow"} /*optional_params*/,
               base / "selected_events");
  read_option(parser, log,
              this->event_data_options.warn_if_selected_event_is_abnormal,
              base / "selected_events" / "warn", true);

  read_option(parser, log,
              this->event_data_options.throw_if_selected_event_is_abnormal,
              base / "selected_events" / "throw", true);

  read_option(parser, log,
              this->event_data_options.n_write_if_encountered_event_is_abnormal,
              base / "selected_events" / "n_write", 100);

  // Read "impact_table_type" (optional, only takes affect for "low_memory")
  this->event_data_options.use_neighborlist_impact_table = true;
  std::string impact_table_type_str = "neighborlist";
  parser.optional(impact_table_type_str, "impact_table_type");
  if (impact_table_type_str == "neighborlist") {
    this->event_data_options.use_neighborlist_impact_table = true;
    log.indent() << "impact_table_type="
                 << "\"neighborlist\"" << std::endl;
  } else if (impact_table_type_str == "relative") {
    this->event_data_options.use_neighborlist_impact_table = false;
    log.indent() << "impact_table_type="
                 << "\"relative\"" << std::endl;
  } else {
    parser.insert_error("impact_table_type",
                        "Invalid impact_table_type: " + impact_table_type_str);
  }

  // Read "assign_allowed_events_only"
  //
  /// \brief If true (default) check if potentially impacted events are allowed
  ///     and only assign them to the event list if they are (adds an
  ///     additional check, but may reduce the size of the event list).
  ///     Otherwise, assign all potentially impacted events to the event list
  ///     (whether they are allowed will still be checked during the rate
  ///     calculation).
  read_option(parser, log, this->event_data_options.assign_allowed_events_only,
              "assign_allowed_events_only", true);

  // -- State saving
  if (parser.self.contains("state_graph")) {
    fs::path base("state_graph");
    check_params(params[base], {} /*required_params*/,
                 {"state_selection", "event_selection", "n_recent_events",
                  "n_states"} /*optional_params*/,
                 base);

    this->event_data_options.state_graph_options = state_graph::Options();
    state_graph::Options &opt = *this->event_data_options.state_graph_options;

    // "state_graph"/"state_selection":
    // - "n_jump_first", "basin_jump_first", "n_look_first", "basin_look_first"
    std::string state_selection_str = "n_jump_first";
    parser.optional(state_selection_str, "state_selection");
    if (state_selection_str == "n_jump_first") {
      opt.state_selection_method =
          state_graph::StateSelectionMethod::N_JUMP_FIRST;
      log.indent() << "state_graph/state_selection="
                   << "\"n_jump_first\"" << std::endl;
    } else if (state_selection_str == "basin_jump_first") {
      opt.state_selection_method =
          state_graph::StateSelectionMethod::BASIN_JUMP_FIRST;
      log.indent() << "state_graph/state_selection="
                   << "\"basin_jump_first\"" << std::endl;
    } else if (state_selection_str == "n_look_first") {
      opt.state_selection_method =
          state_graph::StateSelectionMethod::N_LOOK_FIRST;
      log.indent() << "state_graph/state_selection="
                   << "\"n_look_first\"" << std::endl;
    } else if (state_selection_str == "basin_look_first") {
      opt.state_selection_method =
          state_graph::StateSelectionMethod::BASIN_LOOK_FIRST;
      log.indent() << "state_graph/state_selection="
                   << "\"basin_look_first\"" << std::endl;
    } else {
      parser.insert_error(base / "state_selection",
                          "Invalid state_selection: " + state_selection_str);
    }

    // "state_graph"/"event_selection"
    // - "first_passage_time_analysis", "mean_rate_method", "memory_only"
    std::string event_selection_str = "memory_only";
    parser.optional(event_selection_str, "memory_only");
    if (event_selection_str == "first_passage_time_analysis") {
      opt.event_selection_method =
          state_graph::EventSelectionMethod::FIRST_PASSAGE_TIME_ANALYSIS;
      log.indent() << "state_graph/event_selection="
                   << "\"first_passage_time_analysis\"" << std::endl;
    } else if (event_selection_str == "mean_rate_method") {
      opt.event_selection_method =
          state_graph::EventSelectionMethod::MEAN_RATE_METHOD;
      log.indent() << "state_graph/event_selection="
                   << "\"mean_rate_method\"" << std::endl;
    } else if (event_selection_str == "memory_only") {
      opt.event_selection_method =
          state_graph::EventSelectionMethod::MEMORY_ONLY;
      log.indent() << "state_graph/event_selection="
                   << "\"memory_only\"" << std::endl;
    } else {
      parser.insert_error(base / "event_selection",
                          "Invalid event_selection: " + event_selection_str);
    }

    // "state_graph" / "n_recent_events"
    read_option(parser, log, opt.n_recent_events, base / "n_recent_events",
                100);

    // "state_graph" / "n_states"
    read_option(parser, log, opt.n_states, base / "n_states", 20);
  } else {
    log.indent() << "state_graph=null" << std::endl;
  }

  log << std::endl;
  log.end_section();

  std::stringstream ss;
  ss << "Error in KineticCalculator: error reading calculation "
        "parameters.";
  std::runtime_error error_if_invalid{ss.str()};
  report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);

  // -- After parsing ----------------------------

  // Make event data
  if (this->event_data_type == kinetic_event_data_type::high_memory) {
    if (this->event_selector_type != kinetic_event_selector_type::sum_tree) {
      throw std::runtime_error(
          "Error in KineticCalculator: "
          "for event_data_type=\"high_memory\", "
          "only event_selector_type=\"sum_tree\" is allowed.");
    }

    if (this->verbosity_level == 100) {
      const bool DebugMode = true;
      this->make_complete_event_data_impl<DebugMode>();
    } else {
      const bool DebugMode = false;
      this->make_complete_event_data_impl<DebugMode>();
    }
  } else {
    // event_data_type == default_memory or low_memory
    if (this->verbosity_level == 100) {
      const bool DebugMode = true;
      this->make_allowed_event_data_impl<DebugMode>();
    } else {
      const bool DebugMode = false;
      this->make_allowed_event_data_impl<DebugMode>();
    }
  }

  return;
}

/// \brief Clone the KineticCalculator
KineticCalculator *KineticCalculator::_clone() const {
  return new KineticCalculator(*this);
}

}  // namespace kinetic_2
}  // namespace clexmonte
}  // namespace CASM

extern "C" {
/// \brief Returns a clexmonte::BaseMonteCalculator* owning a KineticCalculator
CASM::clexmonte::BaseMonteCalculator *make_KineticCalculator() {
  return new CASM::clexmonte::kinetic_2::KineticCalculator();
}
}
