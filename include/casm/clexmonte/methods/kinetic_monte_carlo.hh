/// An implementation of a kinetic Monte Carlo main loop
/// that makes use of the RunManager provided by
/// casm/monte/run_management to implement sampling
/// fixtures and results data structures and input/output
/// methods and a data structure that allows sampling
/// atomic displacements for kinetic coefficient
/// calculations.

#ifndef CASM_clexmonte_methods_kinetic_monte_carlo
#define CASM_clexmonte_methods_kinetic_monte_carlo

// logging
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clexmonte/definitions.hh"
#include "casm/configuration/occ_events/OccSystem.hh"
#include "casm/monte/Conversions.hh"
#include "casm/monte/MethodLog.hh"
#include "casm/monte/checks/io/json/CompletionCheck_json_io.hh"
#include "casm/monte/events/OccLocation.hh"
#include "casm/monte/methods/kinetic_monte_carlo.hh"
#include "casm/monte/run_management/State.hh"
#include "casm/monte/sampling/SelectedEventFunctions.hh"

namespace CASM {
namespace clexmonte {

/// \brief Construct a list of atom names corresponding to OccLocation atoms
std::vector<Index> make_atom_name_index_list(
    monte::OccLocation const &occ_location,
    occ_events::OccSystem const &occ_system);

template <bool DebugMode, typename ConfigType,
          typename SetSelectedEventFunction,
          typename ApplySelectedEventFunction, typename StatisticsType,
          typename EngineType>
void kinetic_monte_carlo_v2(
    state_type &state, monte::OccLocation &occ_location,
    monte::KMCData<ConfigType, StatisticsType, EngineType> &kmc_data,
    SelectedEvent &selected_event,
    SetSelectedEventFunction set_selected_event_f,
    ApplySelectedEventFunction apply_selected_event_f,
    std::optional<monte::SelectedEventDataCollector> &collector,
    monte::RunManager<ConfigType, StatisticsType, EngineType> &run_manager,
    std::shared_ptr<occ_events::OccSystem> event_system);

// --- Implementation ---

template <bool Debug, int size = 50>
void begin_section(const std::string &message) {
  if constexpr (Debug) {
    Log &log = CASM::log();
    log.indent() << "## " << message << " "
                 << std::string(size - 4 - message.size(), '#') << std::endl;
  }
}

template <bool Debug, int size = 50>
void end_section() {
  if constexpr (Debug) {
    Log &log = CASM::log();
    log.indent() << std::string(size, '#') << std::endl << std::endl;
  }
}

template <bool Debug, int size = 50>
void debug_collect(monte::SelectedEventDataCollector const &collector) {
  if constexpr (Debug) {
    Log &log = CASM::log();

    jsonParser json;

    // Discrete vector int data
    log.custom("Discrete vector int functions");
    json["discrete_vector_int_functions"].put_obj();
    for (auto const &f : collector.discrete_vector_int_f) {
      log.indent() << "- function: " << f.name << std::endl;
      jsonParser &tjson = json["discrete_vector_int_functions"][f.name];
      log.increase_indent();

      // Check if function requires event state
      bool requires_event_state = f.requires_event_state;
      log.indent() << "- requires_event_state: " << requires_event_state
                   << std::endl;
      tjson["requires_event_state"] = requires_event_state;

      // Check if function has value in current state
      bool has_value = f.has_value();
      log.indent() << "- has_value: " << has_value << std::endl;
      tjson["has_value"] = has_value;

      // If function has a value, evaluate it
      if (has_value) {
        Eigen::VectorXl value = f.function();
        to_json(f.shape, tjson["shape"]);
        to_json(f.component_names, tjson["component_names"]);
        to_json(value, tjson["value"], jsonParser::as_array());
        log.indent() << "- shape: " << tjson["shape"] << std::endl;
        log.indent() << "- component_names: " << tjson["component_names"]
                     << std::endl;
        log.indent() << "- value: " << tjson["value"] << std::endl;

        if (f.value_labels.has_value()) {
          auto const &value_to_label = f.value_labels.value();
          if (value_to_label.count(value)) {
            to_json(value_to_label.at(value), tjson["label"]);
            log.indent() << "- label: " << tjson["label"] << std::endl;
          }
        }
      }
      log.decrease_indent();
    }
    if (collector.discrete_vector_int_f.empty()) {
      log.indent() << "- no functions" << std::endl;
    }
    log << std::endl;

    // Discrete vector float data
    log.custom("Discrete vector float functions");
    json["discrete_vector_float_functions"].put_obj();
    for (auto const &f : collector.discrete_vector_float_f) {
      log.indent() << "- function: " << f.name << std::endl;
      jsonParser &tjson = json["discrete_vector_float_functions"][f.name];
      log.increase_indent();

      // Check if function requires event state
      bool requires_event_state = f.requires_event_state;
      log.indent() << "- requires_event_state: " << requires_event_state
                   << std::endl;
      tjson["requires_event_state"] = requires_event_state;

      // Check if function has value in current state
      bool has_value = f.has_value();
      log.indent() << "- has_value: " << has_value << std::endl;
      tjson["has_value"] = has_value;

      // If function has a value, evaluate it
      if (has_value) {
        Eigen::VectorXd value = f.function();
        to_json(f.shape, tjson["shape"]);
        to_json(f.component_names, tjson["component_names"]);
        to_json(value, tjson["value"], jsonParser::as_array());
        log.indent() << "- shape: " << tjson["shape"] << std::endl;
        log.indent() << "- component_names: " << tjson["component_names"]
                     << std::endl;
        log.indent() << "- value: " << tjson["value"] << std::endl;

        if (f.value_labels.has_value()) {
          auto const &value_to_label = f.value_labels.value();
          if (value_to_label.count(value)) {
            to_json(value_to_label.at(value), tjson["label"]);
            log.indent() << "- label: " << tjson["label"] << std::endl;
          }
        }
      }
      log.decrease_indent();
    }
    if (collector.discrete_vector_float_f.empty()) {
      log.indent() << "- no functions" << std::endl;
    }
    log << std::endl;

    // Continuous 1d data
    log.custom("Continuous 1d functions");
    json["continuous_1d_functions"].put_obj();
    for (auto const &f : collector.continuous_1d_f) {
      log.indent() << "- function: " << f.name << std::endl;
      jsonParser &tjson = json["continuous_1d_functions"][f.name];
      log.increase_indent();

      // Check if function requires event state
      bool requires_event_state = f.requires_event_state;
      log.indent() << "- requires_event_state: " << requires_event_state
                   << std::endl;
      tjson["requires_event_state"] = requires_event_state;

      // partition index
      int partition = f.partition();
      log.indent() << "- partition index: " << partition << std::endl;
      tjson["partition_index"] = partition;

      // validate partition index
      if (partition < 0 || partition >= f.partition_names.size()) {
        std::stringstream ss;
        ss << "Error: partition index (=" << partition << ") out of range.";
        throw std::runtime_error(ss.str());
      }

      // partition name
      std::string partition_name = f.partition_names[partition];
      log.indent() << "- partition name: " << partition_name << std::endl;
      tjson["partition_name"] = partition_name;

      // value
      double value = f.function();
      log.indent() << "- value: " << value << std::endl;
      tjson["value"] = value;

      log.decrease_indent();
    }
    if (collector.continuous_1d_f.empty()) {
      log.indent() << "- no functions" << std::endl;
    }
    log << std::endl;

    // Generic functions
    log.custom("Generic functions");
    for (auto const &f : collector.generic_f) {
      log.indent() << "- function: " << f.name << std::endl;
      jsonParser &tjson = json["generic_functions"][f.name];
      log.increase_indent();

      // Check function evaluation order
      Index order = f.order;
      log.indent() << "- order: " << order << std::endl;
      tjson["order"] = order;

      // Check if function requires event state
      bool requires_event_state = f.requires_event_state;
      log.indent() << "- requires_event_state: " << requires_event_state
                   << std::endl;
      tjson["requires_event_state"] = requires_event_state;

      // Check if function has value in current state
      bool has_value = f.has_value();
      log.indent() << "- has_value: " << has_value << std::endl;
      tjson["has_value"] = has_value;

      // If function has a value, evaluate it
      if (has_value) {
        bool completed = true;
        log.indent() << "- completed: " << true << std::endl;
        tjson["completed"] = true;
      }

      log.decrease_indent();
    }
    if (collector.generic_f.empty()) {
      log.indent() << "- no functions" << std::endl;
    }
    log << std::endl;

    log.custom("JSON Summary");
    log.indent() << json << std::endl << std::endl;
  }
}

/// \brief Construct a list of atom names corresponding to OccLocation atoms
///
/// Notes:
/// - If atoms are conserved, then the order of this list will remain unchanged
///   during the course of a calculation
/// - Values are set to -1 if atom is no longer in supercell
inline std::vector<Index> make_atom_name_index_list(
    monte::OccLocation const &occ_location,
    occ_events::OccSystem const &occ_system) {
  // sanity check:
  monte::Conversions const &convert = occ_location.convert();
  if (convert.species_size() != occ_system.orientation_name_list.size()) {
    throw std::runtime_error(
        "Error in CASM::clexmonte::kinetic::make_snapshot_for_conserved_atoms: "
        "mismatch between monte::Conversions and occ_events::OccSystem.");
  }

  // collect atom name indices
  std::vector<Index> atom_name_index_list(occ_location.atom_size(), -1);
  for (Index i = 0; i < occ_location.mol_size(); ++i) {
    monte::Mol const &mol = occ_location.mol(i);
    Index b = convert.l_to_b(mol.l);
    Index occupant_index =
        occ_system.orientation_to_occupant_index[b][mol.species_index];
    Index atom_position_index = 0;
    for (Index atom_id : mol.component) {
      Index atom_name_index =
          occ_system.atom_position_to_name_index[b][occupant_index]
                                                [atom_position_index];
      atom_name_index_list.at(atom_id) = atom_name_index;
      ++atom_position_index;
    }
  }
  return atom_name_index_list;
}

/// \brief Run a kinetic Monte Carlo calculation
///
/// \param state The state. Consists of both the initial
///     configuration and conditions. Conditions must include `temperature`
///     and any others required by `potential`.
/// \param occ_location An occupant location tracker, which enables efficient
///     event proposal. It must already be initialized with the input state.
/// \param kmc_data Stores data to be made available to the sampling functions
///     along with the current state.
/// \param selected_event The selected event data object which the selected
///     event data functions are expecting to be set with the selected event
///     data.
/// \param set_selected_event_f A function that can set `selected_event` with
///     signature `set_selected_event_f(SelectedEvent &selected_event, bool
///     requires_event_state)`.
/// \param apply_selected_event_f A function that can apply the selected event
///     and prepare to select the next event, with
///     signature `apply_selected_event_f(state_type &state,
///     monte::OccLocation &occ_location, SelectedEvent &selected_event)`.
/// \param collector Collects selected event data
/// \param run_manager Contains sampling fixtures and after completion holds
///     final results
/// \param event_system Defines the system for OccPosition / OccTrajectory /
///     OccEvent. Used in particular to determine
///     `kmc_data.atom_name_index_list`.
///
/// \returns A Results<ConfigType> instance with run results.
///
/// Required interface for `State<ConfigType>`:
/// - `Eigen::VectorXi &get_occupation(State<ConfigType> const &state)`
/// - `Eigen::Matrix3l const &get_transformation_matrix_to_super(
///        State<ConfigType> const &state)`
///
/// State properties that are set:
/// - None
///
template <bool DebugMode, typename ConfigType,
          typename SetSelectedEventFunction,
          typename ApplySelectedEventFunction, typename StatisticsType,
          typename EngineType>
void kinetic_monte_carlo_v2(
    state_type &state, monte::OccLocation &occ_location,
    monte::KMCData<ConfigType, StatisticsType, EngineType> &kmc_data,
    SelectedEvent &selected_event,
    SetSelectedEventFunction set_selected_event_f,
    ApplySelectedEventFunction apply_selected_event_f,
    std::optional<monte::SelectedEventDataCollector> &collector,
    monte::RunManager<ConfigType, StatisticsType, EngineType> &run_manager,
    std::shared_ptr<occ_events::OccSystem> event_system) {
  // Validate existence of sampling fixtures
  if (run_manager.sampling_fixtures.empty()) {
    throw std::runtime_error(
        "Error in clexmonte::kinetic_monte_carlo_v2: no sampling fixtures");
  }

  // Initialize atom positions & time
  kmc_data.sampling_fixture_label.clear();
  kmc_data.sampling_fixture = nullptr;
  kmc_data.total_rate = 0.0;
  kmc_data.time = 0.0;
  kmc_data.unique_atom_id = occ_location.unique_atom_id();
  kmc_data.atom_positions_cart = occ_location.atom_positions_cart();
  kmc_data.prev_time.clear();
  kmc_data.prev_atom_positions_cart.clear();
  kmc_data.prev_unique_atom_id.clear();
  for (auto &fixture_ptr : run_manager.sampling_fixtures) {
    kmc_data.prev_time.emplace(fixture_ptr->label(), kmc_data.time);
    kmc_data.prev_atom_positions_cart.emplace(fixture_ptr->label(),
                                              kmc_data.atom_positions_cart);
    kmc_data.prev_unique_atom_id.emplace(fixture_ptr->label(),
                                         kmc_data.unique_atom_id);
  }

  // Pre- and post- sampling actions:

  // notes: it is important this uses
  // - the total_rate obtained for the state before applying the selected event
  auto pre_sample_action =
      [&](monte::SamplingFixture<ConfigType, StatisticsType, EngineType>
              &fixture,
          monte::State<ConfigType> const &state) {
        // if save states:
        // - determine the state at sample time

        // set data that can be used in sampling functions
        kmc_data.sampling_fixture_label = fixture.label();
        kmc_data.sampling_fixture = &fixture;
        kmc_data.unique_atom_id = occ_location.unique_atom_id();
        kmc_data.atom_name_index_list =
            make_atom_name_index_list(occ_location, *event_system);
        kmc_data.atom_positions_cart = occ_location.atom_positions_cart();
        kmc_data.total_rate = selected_event.total_rate;
        if (fixture.params().sampling_params.sample_mode ==
            monte::SAMPLE_MODE::BY_TIME) {
          kmc_data.time = fixture.next_sample_time();
        }
      };

  auto post_sample_action =
      [&](monte::SamplingFixture<ConfigType, StatisticsType, EngineType>
              &fixture,
          monte::State<ConfigType> const &state) {
        // set data that can be used in sampling functions
        kmc_data.prev_time[fixture.label()] = kmc_data.time;
        kmc_data.prev_atom_positions_cart[fixture.label()] =
            kmc_data.atom_positions_cart;
        kmc_data.prev_unique_atom_id[fixture.label()] = kmc_data.unique_atom_id;

        // if save states:
        // - select next event
      };

  // Main loop
  double event_time;
  bool collect_selected_event_data = collector.has_value();
  selected_event.reset();
  begin_section<DebugMode>("Initialize RunManager");
  run_manager.initialize(occ_location.mol_size());
  run_manager.update_next_sampling_fixture();
  end_section<DebugMode>();

  // Sample data, if a sample is due by count
  // - This location correctly handles sampling at count=0 and count!=0
  // - If the sample count is n, then the state after the n-th step/pass is
  //   sampled.
  begin_section<DebugMode>("Sample by count, if due");
  run_manager.template sample_data_by_count_if_due<DebugMode>(
      state, pre_sample_action, post_sample_action);
  end_section<DebugMode>();

  begin_section<DebugMode, 80>("Check for completion");
  while (!run_manager.is_complete()) {
    // Debug - Write completion check results
    if constexpr (DebugMode) {
      Log &log = CASM::log();
      jsonParser json;
      for (auto &fixture_ptr : run_manager.sampling_fixtures) {
        to_json(fixture_ptr->completion_check_results(),
                json[fixture_ptr->label()]);
      }
      log << json << std::endl << std::endl;
    }
    end_section<DebugMode, 80>();

    begin_section<DebugMode>("Write status, if due");
    run_manager.write_status_if_due();
    end_section<DebugMode>();

    // --

    // Select an event. This function:
    // - If no event groups:
    //   - Select event from event selector
    // - If event groups:
    //   - Select event by finding the group with minimum next event time
    // - Sets the `selected_event` object with the next event and time
    begin_section<DebugMode>("Select an event");
    set_selected_event_f(selected_event);
    end_section<DebugMode>();
    event_time = selected_event.time;

    // Sample data, if a sample is due by time
    // - This location correctly handles sampling at sample times >= 0
    // - Samples are taken using the state before the event for each
    //   sample time <= the event time
    // - If the sample time is exactly equal to the event time (should be
    //   vanishingly rare), then the state before the event occurs is sampled.
    // - While next event time > next sample time:
    //   - If no event groups:
    //     - Set current time to sample time
    //     - Sample
    //   - If event groups: (TODO)
    //     - Resolve state at current time
    //     - Sample
    //     - Calculate next event & time for each group
    //     - Select next event
    begin_section<DebugMode>("Sample by time, if due");
    run_manager.template sample_data_by_time_if_due<DebugMode>(
        event_time, state, pre_sample_action, post_sample_action);
    end_section<DebugMode>();

    // Set time -- for all fixtures and kmc_data
    begin_section<DebugMode>("Update time");
    run_manager.set_time(event_time);
    kmc_data.time = event_time;
    if constexpr (DebugMode) {
      Log &log = CASM::log();
      log.indent() << "- time: " << kmc_data.time << std::endl << std::endl;
    }
    end_section<DebugMode>();

    // Increment count -- for all fixtures
    begin_section<DebugMode>("Update step / pass / count");
    run_manager.increment_step();
    run_manager.increment_n_accept();
    if constexpr (DebugMode) {
      Log &log = CASM::log();

      for (auto &fixture_ptr : run_manager.sampling_fixtures) {
        auto const &counter = fixture_ptr->counter();
        log.indent() << "- sampling fixture: " << fixture_ptr->label()
                     << std::endl;
        log.increase_indent();
        log.indent() << "- step: " << counter.step << std::endl;
        log.indent() << "- pass: " << counter.pass << std::endl;
        log.indent() << "- count: " << counter.count << std::endl;
        log.indent() << "- time: " << counter.time << std::endl;
        log.indent() << "- n_accept: " << counter.n_accept << std::endl;
        log.indent() << "- n_reject: " << counter.n_reject << std::endl;
        log.decrease_indent();
      }
      log.indent() << "- time: " << kmc_data.time << std::endl << std::endl;
    }
    end_section<DebugMode>();

    // If using event groups:
    // - Restore the selected group to exit state
    // - Set the impacted events & groups
    // - Evolve impacted groups to the event time
    // - Regroup impacted events as necessary
    //   - If a new group, add initial state
    //   - Add current site and event information for existing states
    // - Clear the impacted events
    // Evaluate selected event functions
    // Apply event
    // If not using event groups:
    // - Set the impacted events
    // - Update event rates
    // - Clear the impacted events
    // If using event groups:
    // - Set the impacted events
    // - Update event rates
    // - Clear the impacted events
    // - Save new state(s) & transition(s) / Remove old state(s) & transition(s)
    // - Select next event for new and impacted groups
    begin_section<DebugMode>("Apply selected event");
    apply_selected_event_f(state, occ_location, selected_event);
    end_section<DebugMode>();

    // Sample data, if a sample is due by count
    // - This location correctly handles sampling at count=0 and count!=0
    // - If the sample count is n, then the state after the n-th step/pass is
    //   sampled.
    // - If count == next sample count:
    //   - If event groups: Resolve state at current time
    //   - Sample
    //   - If event groups: Calculate next event & time for each group
    begin_section<DebugMode>("Sample by count, if due");
    run_manager.template sample_data_by_count_if_due<DebugMode>(
        state, pre_sample_action, post_sample_action);
    end_section<DebugMode>();

    begin_section<DebugMode, 80>("Check for completion");
  }
  // Debug - Write completion check results
  if constexpr (DebugMode) {
    Log &log = CASM::log();
    jsonParser json;
    for (auto &fixture_ptr : run_manager.sampling_fixtures) {
      to_json(fixture_ptr->completion_check_results(),
              json[fixture_ptr->label()]);
    }
    log << json << std::endl << std::endl;
  }
  end_section<DebugMode, 80>();

  begin_section<DebugMode>("Finalize RunManager");
  run_manager.finalize(state);
  end_section<DebugMode>();
}

}  // namespace clexmonte
}  // namespace CASM

#endif
