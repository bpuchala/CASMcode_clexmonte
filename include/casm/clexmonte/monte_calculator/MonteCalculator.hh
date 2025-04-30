#ifndef CASM_clexmonte_MonteCalculator
#define CASM_clexmonte_MonteCalculator

#include "casm/clexmonte/monte_calculator/BaseMonteCalculator.hh"
#include "casm/clexmonte/monte_calculator/MonteEventData.hh"

namespace CASM {
namespace clexmonte {

/// \brief Implements a potential
class MontePotential {
 public:
  MontePotential(std::shared_ptr<BaseMontePotential> _pot,
                 std::shared_ptr<RuntimeLibrary> _lib)
      : m_pot(_pot), m_lib(_lib) {}

  ~MontePotential() {
    // ensure BaseMontePotential is deleted before library
    m_pot.reset();
  }

  /// State data for sampling functions, for the current state
  std::shared_ptr<StateData> state_data() { return m_pot->state_data; }

  /// \brief Calculate (per_supercell) potential value
  double per_supercell() { return m_pot->per_supercell(); }

  /// \brief Calculate (per_unitcell) potential value
  double per_unitcell() { return m_pot->per_unitcell(); }

  /// \brief Calculate change in (per_supercell) semi-grand potential value due
  ///     to a series of occupation changes
  double occ_delta_per_supercell(std::vector<Index> const &linear_site_index,
                                 std::vector<int> const &new_occ) {
    return m_pot->occ_delta_per_supercell(linear_site_index, new_occ);
  }

 private:
  std::shared_ptr<BaseMontePotential> m_pot;
  std::shared_ptr<RuntimeLibrary> m_lib;
};

/// \brief Wrapper for Monte Carlo calculations implementations
class MonteCalculator {
 public:
  typedef default_engine_type engine_type;
  typedef monte::KMCData<config_type, statistics_type, engine_type>
      kmc_data_type;

  /// \brief Constructor.
  ///
  /// Note: For most uses it is recommended to construct
  /// a std::shared_ptr<MonteCalculator> using the `make_monte_calculator`
  /// factory function.
  ///
  /// \param _base_calculator The underlying implementation
  /// \param _lib If the `base_calculator` is from a runtime library, it should
  ///     be provided to ensure the lifetime of the library. Otherwise, give
  ///     nullptr.
  MonteCalculator(
      std::unique_ptr<clexmonte::BaseMonteCalculator> _base_calculator,
      std::shared_ptr<RuntimeLibrary> _lib)
      : m_calc(_base_calculator), m_lib(_lib) {}

  ~MonteCalculator() {
    // ensure BaseMonteCalculator is deleted before library
    m_calc.reset();
  }

  // --- Set at construction: ---

  /// Calculator name
  std::string const &calculator_name() const { return m_calc->calculator_name; }

  /// Method allows time-based sampling?
  bool time_sampling_allowed() const { return m_calc->time_sampling_allowed; }

  /// Method tracks atom locations? (like in KMC)
  bool update_atoms() const { return m_calc->update_atoms; }

  /// Method tracks all atoms added / removed from supercell? (KMC with
  /// deposition / dissolution)
  bool track_unique_atom_ids() const { return m_calc->track_unique_atom_ids; }

  /// Method saves atom initial / final position info? (KMC with deposition /
  /// dissolution)
  bool save_atom_info() const { return m_calc->save_atom_info; }

  // --- Set at `reset`: ---

  /// \brief Set parameters, check for required system data, and reset derived
  /// Monte Carlo calculator
  void reset(jsonParser const &_params, std::shared_ptr<system_type> system) {
    m_calc->reset(_params, system);
  }

  /// Calculator method parameters
  jsonParser const &params() const { return m_calc->params; }

  /// System data
  std::shared_ptr<system_type> system() const { return m_calc->system; }

  // --- Set by user after `reset`, before `run`: ---

  /// State sampling functions
  std::map<std::string, state_sampling_function_type> sampling_functions;

  /// JSON State sampling functions
  std::map<std::string, json_state_sampling_function_type>
      json_sampling_functions;

  /// Results analysis functions
  std::map<std::string, results_analysis_function_type> analysis_functions;

  /// State modifying functions
  StateModifyingFunctionMap modifying_functions;

  // --- Set when `set_state_and_potential` or `run` is called: ---

  /// State data for sampling functions, for the current state
  std::shared_ptr<StateData> state_data() {
    if (m_calc->state_data == nullptr) {
      throw std::runtime_error(
          "Error in MonteCalculator::state_data: State data is not "
          "yet constructed. To use outside of the `run` method, call "
          "`set_state_and_potential` first.");
    }
    return m_calc->state_data;
  }

  /// \brief Potential calculator
  MontePotential potential() {
    if (m_calc->potential == nullptr) {
      throw std::runtime_error(
          "Error in MonteCalculator::potential: Potential calculator is not "
          "yet constructed. To use outside of the `run` method, call "
          "`set_state_and_potential` first.");
    }
    return MontePotential(m_calc->potential, m_lib);
  }

  /// \brief Validate the state's configuration
  Validator validate_configuration(state_type &state) const {
    return m_calc->validate_configuration(state);
  }

  /// \brief Validate the state's conditions
  Validator validate_conditions(state_type &state) const {
    return m_calc->validate_conditions(state);
  }

  /// \brief Validate the state
  Validator validate_state(state_type &state) const {
    return m_calc->validate_state(state);
  }

  /// \brief Validate and set the current state, construct state_data, construct
  ///     potential
  void set_state_and_potential(state_type &state,
                               monte::OccLocation *occ_location) {
    m_calc->set_state_and_potential(state, occ_location);
  }

  /// \brief Construct and initialize a new occupant location list for the
  ///     current state
  std::shared_ptr<monte::OccLocation> make_occ_location() {
    if (m_calc->state_data == nullptr) {
      throw std::runtime_error(
          "Error in MonteCalculator::make_occ_location: State data is not "
          "yet constructed. To make an occupant location list, call "
          "`set_state_and_potential` first.");
    }
    auto const &state = *m_calc->state_data->state;

    monte::Conversions const &convert =
        get_index_conversions(*m_calc->system, state);
    monte::OccCandidateList const &occ_candidate_list =
        get_occ_candidate_list(*m_calc->system, state);

    auto occ_location = std::make_shared<monte::OccLocation>(
        convert, occ_candidate_list, m_calc->update_atoms,
        m_calc->save_atom_info);

    double time = 0.0;
    occ_location->initialize(get_occupation(state), time);

    m_calc->state_data->owned_occ_location = occ_location;
    m_calc->state_data->occ_location = occ_location.get();

    return occ_location;
  }

  // --- Constructed if applicable ---

  /// Selected event functions (if applicable)
  /// - Constructed by constructor
  std::shared_ptr<monte::SelectedEventFunctions> selected_event_functions() {
    return m_calc->selected_event_functions;
  }

  /// Selected event data collection parameters (if applicable)
  /// - Constructed by `reset`
  std::shared_ptr<monte::SelectedEventFunctionParams>
  selected_event_function_params() {
    return m_calc->selected_event_function_params;
  }

  /// Event data access (if applicable)
  /// - Constructed by `reset`
  /// - Updated by `set_event_data` or `run`
  MonteEventData event_data() {
    if (m_calc->event_data == nullptr) {
      throw std::runtime_error(
          "Error in MonteCalculator::event_data: Event data does not "
          "exist.");
    }
    return MonteEventData(m_calc->event_data, m_lib);
  }

  std::vector<clexmonte::PrimEventData> const &prim_event_list() {
    if (m_calc->event_data == nullptr) {
      throw std::runtime_error(
          "Error in MonteCalculator::prim_event_list: Event data does not "
          "exist.");
    }
    return m_calc->event_data->prim_event_list;
  }

  /// Set selected event data collection paramters
  void set_selected_event_function_params(
      std::shared_ptr<monte::SelectedEventFunctionParams>
          selected_event_function_params) {
    m_calc->selected_event_function_params = selected_event_function_params;
  }

  /// Selected event data access (if applicable)
  std::shared_ptr<monte::SelectedEventData> selected_event_data() {
    return m_calc->selected_event_data;
  }

  /// Selected event access (if applicable)
  std::shared_ptr<SelectedEvent> selected_event() {
    return m_calc->selected_event;
  }

  /// \brief Set event data (includes calculating all rates), using current
  /// state data
  void set_event_data() {
    if (m_calc->state_data == nullptr) {
      throw std::runtime_error(
          "Error in MonteCalculator::set_event_data: State data is not "
          "yet constructed. To construct event data, call "
          "`set_state_and_potential` with an occupant location list first.");
    }
    if (m_calc->state_data->occ_location == nullptr) {
      throw std::runtime_error(
          "Error in MonteCalculator::set_event_data: State data is "
          "constructed, but there is no occupant location list. To construct "
          "event data, call `set_state_and_potential` with an occupant "
          "location list first.");
    }
    m_calc->set_event_data();
  }

  // --- Set backup random number engine: ---

  /// \brief Set the random number engine
  /// \param engine The new random number engine. An exception is raised if
  ///     `engine` is null.
  void set_engine(std::shared_ptr<engine_type> engine) {
    if (engine == nullptr) {
      throw std::runtime_error(
          "Error in MonteCalculator::set_engine: engine is null.");
    }
    m_calc->engine = engine;
  }

  /// Random number engine access (not null)
  std::shared_ptr<engine_type> engine() {
    if (m_calc->engine == nullptr) {
      throw std::runtime_error(
          "Error in MonteCalculator::engine: engine is null");
    }
    return m_calc->engine;
  }

  // --- Set when `run` is called: ---

  /// Run manager access
  std::shared_ptr<run_manager_type<engine_type>> run_manager() {
    if (m_calc->shared_run_manager == nullptr) {
      throw std::runtime_error(
          "Error in MonteCalculator::run_manager: Run manager is not "
          "yet set.");
    }
    return m_calc->shared_run_manager;
  }

  /// KMC data for sampling functions, for the current state (if applicable)
  std::shared_ptr<kmc_data_type> kmc_data() {
    if (m_calc->kmc_data == nullptr) {
      throw std::runtime_error(
          "Error in MonteCalculator::kmc_data: KMC data is not "
          "yet constructed.");
    }
    return m_calc->kmc_data;
  }

  /// \brief Perform a single run, evolving current state
  void run(state_type &state, monte::OccLocation &occ_location,
           std::shared_ptr<run_manager_type<engine_type>> run_manager) {
    if (!run_manager) {
      throw std::runtime_error(
          "Error in MonteCalculator::run: Run manager is not set.");
    }
    m_calc->shared_run_manager = run_manager;
    m_calc->run(state, occ_location, *m_calc->shared_run_manager);
  }

  /// \brief Construct functions that may be used to sample various quantities
  /// of
  ///     the Monte Carlo calculation as it runs
  std::map<std::string, state_sampling_function_type>
  standard_sampling_functions(
      std::shared_ptr<MonteCalculator> const &calculation) const {
    return m_calc->standard_sampling_functions(calculation);
  }

  /// \brief Construct functions that may be used to sample various quantities
  ///     of the Monte Carlo calculation as it runs
  std::map<std::string, json_state_sampling_function_type>
  standard_json_sampling_functions(
      std::shared_ptr<MonteCalculator> const &calculation) const {
    return m_calc->standard_json_sampling_functions(calculation);
  }

  /// \brief Construct functions that may be used to analyze Monte Carlo
  ///     calculation results
  std::map<std::string, results_analysis_function_type>
  standard_analysis_functions(
      std::shared_ptr<MonteCalculator> const &calculation) const {
    return m_calc->standard_analysis_functions(calculation);
  }

  /// \brief Construct functions that may be used to modify states
  StateModifyingFunctionMap standard_modifying_functions(
      std::shared_ptr<MonteCalculator> const &calculation) const {
    return m_calc->standard_modifying_functions(calculation);
  }

  /// \brief Construct functions that may be used to collect selected event data
  std::optional<monte::SelectedEventFunctions>
  standard_selected_event_functions(
      std::shared_ptr<MonteCalculator> const &calculation) const {
    return m_calc->standard_selected_event_functions(calculation);
  }

  /// \brief Construct default SamplingFixtureParams
  sampling_fixture_params_type make_default_sampling_fixture_params(
      std::shared_ptr<MonteCalculator> const &calculation, std::string label,
      bool write_results = true, bool write_trajectory = false,
      bool write_observations = false, bool write_status = true,
      std::optional<std::string> output_dir = std::nullopt,
      std::optional<std::string> log_file = std::nullopt,
      double log_frequency_in_s = 600.0) const {
    return m_calc->make_default_sampling_fixture_params(
        calculation, label, write_results, write_trajectory, write_observations,
        write_status, output_dir, log_file, log_frequency_in_s);
  }

  // --- Experimental, to support multi-state methods: ---

  /// \brief Check if a multi-state method
  bool is_multistate_method() const { return m_calc->is_multistate_method; }

  /// \brief Number of states, for multi-state methods
  int n_states() const { return m_calc->multistate_data.size(); }

  /// \brief Current state index
  int current_state() const { return m_calc->current_state; }

  /// \brief State data for sampling functions, for specified state
  std::shared_ptr<StateData> multistate_state_data(int state_index) {
    return m_calc->multistate_data.at(state_index);
  }

  /// \brief Potential calculator, for specified state
  MontePotential multistate_potential(int state_index) {
    return MontePotential(m_calc->multistate_potential.at(state_index), m_lib);
  }

  /// \brief Perform a single run, evolving one or more states
  void run(int current_state, std::vector<state_type> &states,
           std::vector<monte::OccLocation> &occ_locations,
           std::shared_ptr<run_manager_type<engine_type>> run_manager) {
    if (!run_manager) {
      throw std::runtime_error(
          "Error in MonteCalculator::run: Run manager is not set.");
    }
    m_calc->shared_run_manager = run_manager;
    m_calc->run(current_state, states, occ_locations,
                *m_calc->shared_run_manager);
  }

 private:
  notstd::cloneable_ptr<BaseMonteCalculator> m_calc;
  std::shared_ptr<RuntimeLibrary> m_lib;
};

/// \brief MonteCalculator factory function
std::shared_ptr<MonteCalculator> make_monte_calculator(
    jsonParser const &params, std::shared_ptr<system_type> system,
    std::shared_ptr<MonteCalculator::engine_type> engine,
    std::unique_ptr<BaseMonteCalculator> base_calculator,
    std::shared_ptr<RuntimeLibrary> lib);

/// \brief MonteCalculator factory function, from source
std::shared_ptr<MonteCalculator> make_monte_calculator_from_source(
    fs::path dirpath, std::string calculator_name,
    std::shared_ptr<system_type> system, jsonParser const &params,
    std::shared_ptr<MonteCalculator::engine_type> engine,
    std::string compile_options, std::string so_options);

Eigen::VectorXd scalar_conditions(
    std::shared_ptr<MonteCalculator> const &calculation, std::string key);

Eigen::VectorXd vector_conditions(
    std::shared_ptr<MonteCalculator> const &calculation, std::string key);

Eigen::VectorXd matrix_conditions(
    std::shared_ptr<MonteCalculator> const &calculation, std::string key);

system_type const &get_system(
    std::shared_ptr<MonteCalculator> const &calculation);

state_type const &get_state(
    std::shared_ptr<MonteCalculator> const &calculation);

std::vector<PrimEventData> const &get_prim_event_list(
    std::shared_ptr<MonteCalculator> const &calculation);

/// \brief Make temporary monte::OccLocation if necessary
void make_temporary_if_necessary(state_type const &state,
                                 monte::OccLocation *&occ_location,
                                 std::unique_ptr<monte::OccLocation> &tmp,
                                 MonteCalculator const &calculation);

}  // namespace clexmonte
}  // namespace CASM

#endif
