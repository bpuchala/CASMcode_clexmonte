#include "casm/clexmonte/methods/occupation_metropolis.hh"
#include "casm/clexmonte/monte_calculator/BaseMonteCalculator.hh"
#include "casm/clexmonte/monte_calculator/MonteCalculator.hh"
#include "casm/clexmonte/monte_calculator/analysis_functions.hh"
#include "casm/clexmonte/monte_calculator/sampling_functions.hh"
#include "casm/clexmonte/run/functions.hh"
#include "casm/configuration/io/json/Configuration_json_io.hh"
#include "casm/monte/events/OccEventProposal.hh"
#include "casm/monte/sampling/RequestedPrecisionConstructor.hh"

namespace CASM {
namespace clexmonte {

/// \brief Propose and apply semi-grand canonical events
class SemiGrandCanonicalEventGenerator {
 public:
  typedef BaseMonteCalculator::engine_type engine_type;

  /// \brief Constructor
  ///
  /// Notes:
  /// - One and only one of `_semigrand_canonical_swaps` and
  /// `_semigrand_canonical_multiswaps` should have size != 0
  ///
  /// \param _semigrand_canonical_swaps Single site swap types for semi-grand
  ///     canonical Monte Carlo events. If size > 0, only these events are
  ///     proposed.
  /// \param _semigrand_canonical_multiswaps Multiple site swap types for
  ///     semi-grand canonical Monte Carlo events, such as charge neutral
  ///     events. These events are only proposed if no single swaps are
  ///     provided.
  SemiGrandCanonicalEventGenerator(
      std::vector<monte::OccSwap> const &_semigrand_canonical_swaps,
      std::vector<monte::MultiOccSwap> const &_semigrand_canonical_multiswaps)
      : state(nullptr),
        occ_location(nullptr),
        semigrand_canonical_swaps(_semigrand_canonical_swaps),
        semigrand_canonical_multiswaps(_semigrand_canonical_multiswaps),
        use_multiswaps(semigrand_canonical_swaps.size() == 0) {
    if (semigrand_canonical_swaps.size() == 0 &&
        semigrand_canonical_multiswaps.size() == 0) {
      throw std::runtime_error(
          "Error in SemiGrandCanonicalEventGenerator: "
          "semigrand_canonical_swaps.size() == 0 && "
          "semigrand_canonical_multiswaps.size() == 0");
    }
    if (semigrand_canonical_swaps.size() != 0 &&
        semigrand_canonical_multiswaps.size() != 0) {
      throw std::runtime_error(
          "Error in SemiGrandCanonicalEventGenerator: "
          "semigrand_canonical_swaps.size() != 0 && "
          "semigrand_canonical_multiswaps.size() != 0");
    }
  }

  /// \brief The current state for which events are proposed and applied. Can be
  ///     nullptr, but must be set for use.
  state_type *state;

  /// Occupant tracker
  monte::OccLocation *occ_location;

  /// \brief Single swap types for semi-grand canonical Monte Carlo events
  std::vector<monte::OccSwap> semigrand_canonical_swaps;

  /// \brief Multiple swap types for semi-grand canonical Monte Carlo events
  std::vector<monte::MultiOccSwap> semigrand_canonical_multiswaps;

  /// \brief If true, propose events from multiswaps, else propose events from
  /// single swaps
  bool use_multiswaps;

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
    this->state =
        throw_if_null(_state,
                      "Error in SemiGrandCanonicalEventGenerator::set: "
                      "_state==nullptr");
    this->occ_location =
        throw_if_null(_occ_location,
                      "Error in SemiGrandCanonicalEventGenerator::set: "
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
    if (this->use_multiswaps) {
      return monte::propose_semigrand_canonical_multiswap_event(
          this->occ_event, *this->occ_location,
          this->semigrand_canonical_multiswaps, random_number_generator);
    } else {
      return monte::propose_semigrand_canonical_event(
          this->occ_event, *this->occ_location, this->semigrand_canonical_swaps,
          random_number_generator);
    }
  }

  /// \brief Update the occupation of the current state using the provided event
  void apply(monte::OccEvent const &e) {
    this->occ_location->apply(e, get_occupation(*this->state));
  }
};

class SemiGrandCanonicalPotential : public BaseMontePotential {
 public:
  SemiGrandCanonicalPotential(std::shared_ptr<StateData> _state_data)
      : BaseMontePotential(_state_data),
        state(*state_data->state),
        n_unitcells(state_data->n_unitcells),
        occupation(get_occupation(state)),
        convert(*state_data->convert),
        composition_calculator(
            get_composition_calculator(*this->state_data->system)),
        composition_converter(
            get_composition_converter(*this->state_data->system)),
        param_chem_pot(state.conditions.vector_values.at("param_chem_pot")),
        formation_energy_clex(
            get_clex(*state_data->system, state, "formation_energy")) {
    if (param_chem_pot.size() !=
        composition_converter.independent_compositions()) {
      throw std::runtime_error(
          "Error in SemiGrandCanonicalPotential: param_chem_pot size error");
    }

    matrix_R_transpose = composition_converter.dparam_dmol();

    delta_N.resize(composition_converter.components().size());
  }

  // --- Data used in the potential calculation: ---

  state_type const &state;
  Index n_unitcells;
  Eigen::VectorXi const &occupation;
  monte::Conversions const &convert;
  composition::CompositionCalculator const &composition_calculator;
  composition::CompositionConverter const &composition_converter;
  std::shared_ptr<clexulator::ClusterExpansion> formation_energy_clex;

  Eigen::VectorXd const &param_chem_pot;

  // :math:`R^{\mathsf{T}}`, where
  //
  //     \vec{x} = R^{\mathsf{T}} ( \vec{n} - \vec{n}_0 )
  //
  Eigen::MatrixXd matrix_R_transpose;

  // :math:`\vec{N}`, the change in number of each component per supercell
  Eigen::VectorXd delta_N;

  /// \brief Calculate (per_supercell) potential value
  double per_supercell() override {
    Eigen::VectorXd mol_composition =
        composition_calculator.mean_num_each_component(occupation);
    Eigen::VectorXd param_composition =
        composition_converter.param_composition(mol_composition);

    return formation_energy_clex->per_supercell() -
           n_unitcells * param_chem_pot.dot(param_composition);
  }

  /// \brief Calculate (per_unitcell) potential value
  double per_unitcell() override { return this->per_supercell() / n_unitcells; }

  /// \brief Calculate change in (per_supercell) semi-grand potential value due
  ///     to a series of occupation changes
  double occ_delta_per_supercell(std::vector<Index> const &linear_site_index,
                                 std::vector<int> const &new_occ) override {
    // Epot = (E - N_unit mu_x.T x)
    // dEpot = (E_f - N_unit * mu_x.T * x_f) - (E_i - N * mu_x.T * x_i)
    //      = (E_f - E_i) - N_unit * mu_x.T * (x_f - x_i)
    //      = dE - N * mu_x.T * dx
    //
    // x = R.T * (n - n_0)
    // dx = R.T * dn
    //
    // dEpot = dE - mu_x.T * R.T * dN

    double delta_formation_energy =
        formation_energy_clex->occ_delta_value(linear_site_index, new_occ);

    delta_N.setZero();
    for (Index i = 0; i < linear_site_index.size(); ++i) {
      Index l = linear_site_index[i];
      Index asym = convert.l_to_asym(l);
      Index curr_species = convert.species_index(asym, occupation(l));
      Index new_species = convert.species_index(asym, new_occ[i]);
      delta_N[curr_species] += -1.0;
      delta_N[new_species] += 1.0;
    }

    return delta_formation_energy -
           param_chem_pot.dot(matrix_R_transpose * delta_N);
  }
};

class SemiGrandCanonicalCalculator : public BaseMonteCalculator {
 public:
  using BaseMonteCalculator::engine_type;

  SemiGrandCanonicalCalculator()
      : BaseMonteCalculator("SemiGrandCanonicalCalculator",  // calculator_name
                            {},                    // required_basis_set,
                            {},                    // required_local_basis_set,
                            {"formation_energy"},  // required_clex,
                            {},                    // required_multiclex,
                            {},                    // required_local_clex,
                            {},                    // required_local_multiclex,
                            {},                    // required_dof_spaces,
                            {},                    // required_params,
                            {},                    // optional_params,
                            false,                 // time_sampling_allowed,
                            false,                 // update_atoms,
                            false,                 // track_unique_atom_ids,
                            false,                 // save_atom_info,
                            false                  // is_multistate_method,
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

    // Specific to semi-grand canonical
    functions.push_back(monte_calculator::make_param_chem_pot_f(calculation));

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
        monte_calculator::make_heat_capacity_f(calculation),
        monte_calculator::make_mol_susc_f(calculation),
        monte_calculator::make_param_susc_f(calculation),
        monte_calculator::make_mol_thermochem_susc_f(calculation),
        monte_calculator::make_param_thermochem_susc_f(calculation)};

    std::map<std::string, results_analysis_function_type> function_map;
    for (auto const &f : functions) {
      function_map.emplace(f.name, f);
    }
    return function_map;
  }

  /// \brief Construct functions that may be used to modify states
  StateModifyingFunctionMap standard_modifying_functions(
      std::shared_ptr<MonteCalculator> const &calculation) const override {
    return StateModifyingFunctionMap();
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

      converge(calculation->sampling_functions, completion_check_params)
          .set_abs_precision("potential_energy", 0.001)
          .set_abs_precision("param_composition", 0.001);
    }

    std::vector<std::string> analysis_names = {
        "heat_capacity", "mol_susc", "param_susc", "mol_thermochem_susc",
        "param_thermochem_susc"};

    return clexmonte::make_sampling_fixture_params(
        label, calculation->sampling_functions,
        calculation->json_sampling_functions, calculation->analysis_functions,
        sampling_params, completion_check_params, analysis_names, write_results,
        write_trajectory, write_observations, write_status, output_dir,
        log_file, log_frequency_in_s);
  }

  /// \brief Validate the state's configuration (all are valid)
  Validator validate_configuration(state_type &state) const override {
    return Validator{};
  }

  /// \brief Validate state's conditions
  ///
  /// Notes:
  /// - requires scalar temperature
  /// - requires vector param_chem_pot
  /// - warnings if other conditions are present
  Validator validate_conditions(state_type &state) const override {
    // validate state.conditions
    monte::ValueMap const &conditions = state.conditions;
    Validator v;
    v.insert(validate_keys(conditions.scalar_values,
                           {"temperature"} /*required*/, {} /*optional*/,
                           "scalar", "condition", false /*throw_if_invalid*/));
    v.insert(validate_keys(conditions.vector_values,
                           {"param_chem_pot"} /*required*/, {} /*optional*/,
                           "vector", "condition", false /*throw_if_invalid*/));

    return v;
  }

  /// \brief Validate state
  Validator validate_state(state_type &state) const override {
    Validator v;
    v.insert(this->validate_configuration(state));
    v.insert(this->validate_conditions(state));
    return v;
  }

  /// \brief Validate and set the current state, construct state_data, construct
  ///     potential
  void set_state_and_potential(state_type &state,
                               monte::OccLocation *occ_location) override {
    // Validate system
    if (this->system == nullptr) {
      throw std::runtime_error(
          "Error in SemiGrandCanonicalCalculator::run: system==nullptr");
    }

    // Validate state
    Validator v = this->validate_state(state);
    print(CASM::log(), v);
    if (!v.valid()) {
      throw std::runtime_error(
          "Error in SemiGrandCanonicalCalculator::run: Invalid initial state");
    }

    // Make state data
    this->state_data =
        std::make_shared<StateData>(this->system, &state, occ_location);

    // Make potential calculator
    this->potential =
        std::make_shared<SemiGrandCanonicalPotential>(this->state_data);
  }

  /// \brief Set event data (includes calculating all rates), using current
  /// state data
  void set_event_data() override {
    throw std::runtime_error(
        "Error in SemiGrandCanonicalCalculator::set_event_data: not valid");
  }

  /// \brief Perform a single run, evolving current state
  void run(state_type &state, monte::OccLocation &occ_location,
           run_manager_type<engine_type> &run_manager) override {
    this->set_state_and_potential(state, &occ_location);

    // Get temperature
    double temperature = state.conditions.scalar_values.at("temperature");

    auto potential_occ_delta_per_supercell_f =
        [=](monte::OccEvent const &event) {
          return this->potential->occ_delta_per_supercell(
              event.linear_site_index, event.new_occ);
        };

    // Random number generator
    if (run_manager.engine == nullptr) {
      throw std::runtime_error(
          "Error in SemiGrandCanonicalCalculator::run: "
          "run_manager.engine==nullptr");
    }
    this->engine = run_manager.engine;
    monte::RandomNumberGenerator<engine_type> random_number_generator(
        run_manager.engine);

    // Make event generator
    auto event_generator = std::make_shared<SemiGrandCanonicalEventGenerator>(
        get_semigrand_canonical_swaps(*this->system),
        get_semigrand_canonical_multiswaps(*this->system));
    event_generator->set(&state, &occ_location);

    auto propose_event_f =
        [=](monte::RandomNumberGenerator<engine_type> &random_number_generator)
        -> monte::OccEvent const & {
      return event_generator->propose(random_number_generator);
    };

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
        "Error: SemiGrandCanonicalCalculator does not allow multi-state runs");
  }

  // --- Parameters ---
  int verbosity_level = 10;

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
  void _reset() override {
    ParentInputParser parser{params};

    // "verbosity": str or int, default=10
    this->verbosity_level = parse_verbosity(parser);
    CASM::log().set_verbosity(this->verbosity_level);

    // TODO: enumeration

    std::stringstream ss;
    ss << "Error in SemiGrandCanonicalCalculator: error reading calculation "
          "parameters.";
    std::runtime_error error_if_invalid{ss.str()};
    report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);

    return;
  }

  /// \brief Clone the SemiGrandCanonicalCalculator
  SemiGrandCanonicalCalculator *_clone() const override {
    return new SemiGrandCanonicalCalculator(*this);
  }
};

}  // namespace clexmonte
}  // namespace CASM

extern "C" {
/// \brief Returns a clexmonte::BaseMonteCalculator* owning a
/// SemiGrandCanonicalCalculator
CASM::clexmonte::BaseMonteCalculator *make_SemiGrandCanonicalCalculator() {
  return new CASM::clexmonte::SemiGrandCanonicalCalculator();
}
}
