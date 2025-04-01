#include "casm/clexmonte/monte_calculator/kinetic_events.hh"

#include "casm/casm_io/SafeOfstream.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/events/io/json/event_data_json_io.hh"
#include "casm/clexmonte/events/io/stream/EventState_stream_io.hh"
#include "casm/clexmonte/misc/to_json.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/state/io/json/State_json_io.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/configuration/Configuration.hh"
#include "casm/crystallography/io/UnitCellCoordIO.hh"
#include "casm/monte/events/OccLocation.hh"
#include "casm/monte/run_management/RunManager.hh"
#include "casm/monte/run_management/State.hh"

// this must be included after the classes are defined for its application here
#include "casm/clexmonte/events/state_graph_impl.hh"
#include "casm/clexmonte/methods/kinetic_monte_carlo.hh"

namespace CASM {
namespace clexmonte {
namespace kinetic_2 {

namespace {

template <bool DebugMode>
bool check_requires_event_state(
    std::optional<monte::SelectedEventDataCollector> &collector,
    bool selected_abnormal_event_handling_on) {
  bool requires_event_state =
      (collector.has_value() && collector->requires_event_state) ||
      selected_abnormal_event_handling_on;

  if constexpr (DebugMode) {
    Log &log = CASM::log();
    log.custom(
        "Check if event selection requires re-calculating the event state:");

    // Check 1a: Selected event functions exist
    log.indent() << "- Selected event functions exist=" << std::boolalpha
                 << collector.has_value() << std::endl;

    // Check 1b: Selected event functions require event state
    if (collector.has_value()) {
      log.indent() << "- Selected event functions require event state="
                   << std::boolalpha << collector->requires_event_state
                   << std::endl;
    }

    // Check 2:
    log.indent() << "- selected_abnormal_event_handling_on=" << std::boolalpha
                 << selected_abnormal_event_handling_on << std::endl;

    // Result:
    log.indent() << "- Event selection requires re-calculating the event state="
                 << std::boolalpha << requires_event_state << std::endl
                 << std::endl;
    log.end_section();
  }

  return requires_event_state;
}

}  // namespace

// -- CompleteKineticEventData --

template <bool DebugMode>
CompleteEventCalculator<DebugMode>::CompleteEventCalculator(
    std::vector<PrimEventData> const &_prim_event_list,
    std::vector<EventStateCalculator> const &_prim_event_calculators,
    std::map<EventID, EventData> const &_event_list,
    bool _abnormal_event_handling_on,
    AbnormalEventHandlingFunction &_handling_f,
    std::map<std::string, Index> &_n_encountered_abnormal)
    : prim_event_list(_prim_event_list),
      prim_event_calculators(_prim_event_calculators),
      event_list(_event_list),
      abnormal_event_handling_on(_abnormal_event_handling_on),
      handling_f(_handling_f),
      n_encountered_abnormal(_n_encountered_abnormal) {}

/// \brief Update `event_state` for event `id` in the current state and
/// return the event rate

template <bool DebugMode>
double CompleteEventCalculator<DebugMode>::calculate_rate(EventID const &id) {
  EventData const &event_data = event_list.at(id);
  PrimEventData const &prim_event_data =
      prim_event_list.at(id.prim_event_index);
  // Note: to keep all event state calculations, uncomment this:
  // EventState &event_state = event_data.event_state;
  prim_event_calculators.at(id.prim_event_index)
      .calculate_event_state(event_state, event_data.unitcell_index,
                             event_data.event.linear_site_index,
                             prim_event_data);

  // ---
  // can check event state and handle non-normal event states here
  // ---
  if (abnormal_event_handling_on) {
    if (event_state.is_allowed && !event_state.is_normal) {
      if constexpr (DebugMode) {
        Log &log = CASM::log();
        log.custom("Handle encountered abnormal event...");
        log.indent() << "- event_type_name=" << prim_event_data.event_type_name
                     << std::endl;
        log.indent() << "Handling encountered abnormal event..." << std::endl;
      }
      Index &n = n_encountered_abnormal[prim_event_data.event_type_name];
      n += 1;
      handling_f(n, event_state, event_data, prim_event_data,
                 *prim_event_calculators.at(id.prim_event_index).state());

      if constexpr (DebugMode) {
        Log &log = CASM::log();
        log.indent() << "Handling encountered abnormal event... DONE"
                     << std::endl;
        log.end_section();
      }
    }
  }

  return event_state.rate;
}

template <bool DebugMode>
CompleteKineticEventData<DebugMode>::CompleteKineticEventData(
    std::shared_ptr<system_type> _system,
    std::optional<std::vector<EventFilterGroup>> _event_filters,
    EventDataOptions _options)
    : options(_options),
      transformation_matrix_to_super(Eigen::Matrix3l::Zero(3, 3)) {
  if constexpr (DebugMode) {
    Log &log = CASM::log();
    log.custom("Construct CompleteKineticEventData");
    log.end_section();
  }

  system = _system;
  if (!is_clex_data(*system, "formation_energy")) {
    throw std::runtime_error(
        "Error constructing CompleteKineticEventData: no 'formation_energy' "
        "clex.");
  }

  prim_event_list = clexmonte::make_prim_event_list(*system);
  if (prim_event_list.empty()) {
    throw std::runtime_error(
        "Error constructing CompleteKineticEventData: "
        "prim event list is empty.");
  }
  if constexpr (DebugMode) {
    Log &log = CASM::log();
    log.custom("Prim event list");
    log.indent() << qto_json(prim_event_list) << std::endl << std::endl;
  }

  prim_impact_info_list = clexmonte::make_prim_impact_info_list(
      *system, prim_event_list, {"formation_energy"});

  if (_event_filters.has_value()) {
    event_filters = _event_filters.value();
  }

  set_encountered_abnormal_event_handling(BasicAbnormalEventHandler(
      "encountered" /*std::string _event_kind*/,
      options.throw_if_encountered_event_is_abnormal /*bool _do_throw*/,
      options.warn_if_encountered_event_is_abnormal /*bool _do_warn*/,
      options.disallow_if_encountered_event_is_abnormal /*bool _disallow*/,
      options.n_write_if_encountered_event_is_abnormal /*int _n_write*/,
      options.output_dir /*fs::path _output_dir*/,
      options.local_corr_compare_tol /*double _tol*/));

  set_selected_abnormal_event_handling(BasicAbnormalEventHandler(
      "selected" /*std::string _event_kind*/,
      options.throw_if_selected_event_is_abnormal /*bool _do_throw*/,
      options.warn_if_selected_event_is_abnormal /*bool _do_warn*/,
      false /*bool _disallow*/,
      options.n_write_if_selected_event_is_abnormal /*int _n_write*/,
      options.output_dir /*fs::path _output_dir*/,
      options.local_corr_compare_tol /*double _tol*/));

  /// Construct the state graph if requested
  if (_options.state_graph_options.has_value()) {
    this->state_graph = std::make_shared<state_graph::StateGraph>(
        *_options.state_graph_options);
  }
}

/// \brief Update for given state, conditions, occupants, event filters
///
/// Notes:
/// - This constructs the complete event list and impact table, and constructs
///   the event selector, which calculates all event rates.
/// - If there are no event filters and the supercell remains unchanged from
///   the previous update, then the event list and impact table are not
///   re-constructed, but the event rates are still re-calculated.
/// - Resets the `n_encountered_abnormal` and `n_selected_abnormal`
///   counters.
template <bool DebugMode>
void CompleteKineticEventData<DebugMode>::update(
    std::shared_ptr<StateData> _state_data,
    std::optional<std::vector<EventFilterGroup>> _event_filters,
    std::shared_ptr<engine_type> engine) {
  // Current state info
  state_data = _state_data;
  state_type const &state = *state_data->state;
  monte::OccLocation const &occ_location = *state_data->occ_location;

  // if same supercell && no event filters
  // -> just re-set state & avoid re-constructing event list
  if (this->transformation_matrix_to_super ==
          get_transformation_matrix_to_super(state) &&
      !_event_filters.has_value()) {
    for (auto &event_state_calculator : prim_event_calculators) {
      event_state_calculator.set(&state);
    }
  } else {
    if (_event_filters.has_value()) {
      event_filters = _event_filters.value();
    }

    // These are constructed/re-constructed so cluster expansions point
    // at the current state
    prim_event_calculators.clear();
    for (auto const &prim_event_data : prim_event_list) {
      prim_event_calculators.emplace_back(system,
                                          prim_event_data.event_type_name);
      prim_event_calculators.back().set(&state);

      // Set a custom event state calculation function if it exists:
      auto it = custom_event_state_calculation_f.find(
          prim_event_data.event_type_name);
      if (it != custom_event_state_calculation_f.end()) {
        prim_event_calculators.back().set_custom_event_state_calculation(
            it->second);
      }
    }

    // Construct CompleteEventList
    event_list = clexmonte::make_complete_event_list(
        prim_event_list, prim_impact_info_list, occ_location, event_filters);

    // Reset "not normal" event counters
    n_encountered_abnormal.clear();
    n_selected_abnormal.clear();

    // Construct CompleteEventCalculator
    if (encountered_abnormal_event_handling_on == true &&
        encountered_abnormal_event_handling_f == nullptr) {
      throw std::runtime_error(
          "Error in CompleteKineticEventData::update: "
          "encountered_abnormal_event_handling_on == true && "
          "encountered_abnormal_event_handling_f == nullptr");
    }
    if (selected_abnormal_event_handling_on == true &&
        selected_abnormal_event_handling_f == nullptr) {
      throw std::runtime_error(
          "Error in CompleteKineticEventData::update: "
          "selected_abnormal_event_handling_on == true && "
          "selected_abnormal_event_handling_f == nullptr");
    }
    event_calculator = std::make_shared<CompleteEventCalculator<DebugMode>>(
        prim_event_list, prim_event_calculators, event_list.events,
        encountered_abnormal_event_handling_on,
        encountered_abnormal_event_handling_f, n_encountered_abnormal);

    transformation_matrix_to_super = get_transformation_matrix_to_super(state);
  }

  Index n_unitcells = transformation_matrix_to_super.determinant();

  // Make event selector
  // - This calculates all rates at construction
  event_selector =
      std::make_shared<CompleteKineticEventData::event_selector_type>(
          event_calculator,
          clexmonte::make_complete_event_id_list(n_unitcells, prim_event_list),
          event_list.impact_table,
          std::make_shared<lotto::RandomGeneratorT<engine_type>>(engine));
}

template <bool DebugMode>
void CompleteKineticEventData<DebugMode>::run(
    state_type &state, monte::OccLocation &occ_location,
    kmc_data_type &kmc_data, SelectedEvent &selected_event,
    std::optional<monte::SelectedEventDataCollector> &collector,
    run_manager_type &run_manager,
    std::shared_ptr<occ_events::OccSystem> event_system) {
  // Function to set selected event
  bool requires_event_state = check_requires_event_state<DebugMode>(
      collector, this->selected_abnormal_event_handling_on);
  auto set_selected_event_f = [=](SelectedEvent &selected_event) {
    this->select_event(selected_event, requires_event_state);
  };

  auto set_impacted_events_f = [=](SelectedEvent &selected_event) {
    // Set impacted events
    this->event_selector->set_impacted_events(selected_event.event_id);
  };

  // Run Kinetic Monte Carlo at a single condition
  kinetic_monte_carlo_v2<DebugMode>(
      state, occ_location, kmc_data, selected_event, set_selected_event_f,
      set_impacted_events_f, collector, run_manager, event_system);
}

/// \brief Update for given state, conditions, occupants, event filters
template <bool DebugMode>
void CompleteKineticEventData<DebugMode>::select_event(
    SelectedEvent &selected_event, bool requires_event_state) {
  // This function:
  // - Updates rates of events impacted by the *last* selected event (if there
  //   was a previous selection)
  // - Updates the total rate
  // - Chooses an event and time increment (does not apply event)
  //
  // It does not apply the event or set the impacted events.

  std::tie(selected_event.event_id, selected_event.time_increment) =
      event_selector->only_select_event();
  selected_event.total_rate = event_selector->total_rate();
  EventID const &event_id = selected_event.event_id;
  EventData const &event_data = event_list.events.at(event_id);
  PrimEventData const &prim_event_data =
      prim_event_list[event_id.prim_event_index];
  selected_event.event_data = &event_data;
  selected_event.prim_event_data = &prim_event_data;

  if (requires_event_state) {
    EventStateCalculator &prim_event_calculator =
        prim_event_calculators.at(event_id.prim_event_index);
    prim_event_calculator.calculate_event_state(
        m_event_state, event_data.unitcell_index,
        event_data.event.linear_site_index, prim_event_data);
    selected_event.event_state = &m_event_state;

    if (selected_abnormal_event_handling_on && !m_event_state.is_normal) {
      if constexpr (DebugMode) {
        Log &log = CASM::log();
        log.custom("Handle selected abnormal event...");
        log.indent() << "- event_type_name=" << prim_event_data.event_type_name
                     << std::endl;
        log.indent() << "Handling selected abnormal event ..." << std::endl;
      }
      Index &n = n_selected_abnormal[prim_event_data.event_type_name];
      n += 1;
      selected_abnormal_event_handling_f(n, m_event_state, event_data,
                                         prim_event_data,
                                         *prim_event_calculator.state());

      if constexpr (DebugMode) {
        Log &log = CASM::log();
        log.indent() << "Handling selected abnormal event... DONE" << std::endl;
        log.end_section();
      }
    }
  }
}

// -- AllowedKineticEventData --

template <bool DebugMode>
AllowedEventCalculator<DebugMode>::AllowedEventCalculator(
    std::vector<PrimEventData> const &_prim_event_list,
    std::vector<EventStateCalculator> const &_prim_event_calculators,
    AllowedEventList &_event_list, bool _abnormal_event_handling_on,
    AbnormalEventHandlingFunction &_handling_f,
    std::map<std::string, Index> &_n_encountered_abnormal)
    : prim_event_list(_prim_event_list),
      prim_event_calculators(_prim_event_calculators),
      event_list(_event_list),
      abnormal_event_handling_on(_abnormal_event_handling_on),
      handling_f(_handling_f),
      n_encountered_abnormal(_n_encountered_abnormal) {}

/// \brief Update `event_state` for event `event_index` in the current state
/// and return the event rate; if the event is no longer allowed, free the
/// event.
template <bool DebugMode>
double AllowedEventCalculator<DebugMode>::calculate_rate(Index event_index) {
  AllowedEventData const &allowed_event_data =
      event_list.allowed_event_map.events()[event_index];
  // EventID original_event_id = allowed_event_data.event_id;
  if (!allowed_event_data.is_assigned) {
    event_state.is_allowed = false;
    event_state.rate = 0.0;
  } else {
    this->calculate_rate(allowed_event_data.event_id);

    // free event from AllowedEventList if not allowed
    if (!event_state.is_allowed) {
      event_list.allowed_event_map.free(allowed_event_data.event_id);
    }
  }

  return event_state.rate;
}

/// \brief Update `event_state` for any event `event_id` in the current state
/// and return the event rate
template <bool DebugMode>
double AllowedEventCalculator<DebugMode>::calculate_rate(
    EventID const &event_id) {
  Index prim_event_index = event_id.prim_event_index;
  PrimEventData const &prim_event_data =
      this->prim_event_list[prim_event_index];
  event_data.unitcell_index = event_id.unitcell_index;

  // set linear_site_index
  set_event_linear_site_index(
      event_data.event.linear_site_index, event_data.unitcell_index,
      event_list.neighbor_index[prim_event_index], *event_list.supercell_nlist);

  // calculate event state
  prim_event_calculators.at(prim_event_index)
      .calculate_event_state(event_state, event_data.unitcell_index,
                             event_data.event.linear_site_index,
                             prim_event_data);

  // ---
  // can check event state and handle non-normal event states here
  // ---
  if (abnormal_event_handling_on) {
    if (event_state.is_allowed && !event_state.is_normal) {
      if constexpr (DebugMode) {
        Log &log = CASM::log();
        log.custom("Handle encountered abnormal event...");
        log.indent() << "- event_type_name=" << prim_event_data.event_type_name
                     << std::endl;
        log.indent() << "Handling encountered abnormal event..." << std::endl;
      }
      Index &n = n_encountered_abnormal[prim_event_data.event_type_name];
      n += 1;
      handling_f(n, event_state, event_data, prim_event_data,
                 *prim_event_calculators.at(prim_event_index).state());

      if constexpr (DebugMode) {
        Log &log = CASM::log();
        log.indent() << "Handling encountered abnormal event... DONE"
                     << std::endl;
        log.end_section();
      }
    }
  }

  return event_state.rate;
}

/// \brief Set `event_data` for event `event_index`, returning a reference
/// which is valid until the next call to this method
template <bool DebugMode>
EventData const &AllowedEventCalculator<DebugMode>::set_event_data(
    Index event_index) {
  return set_event_data(event_list.allowed_event_map.event_id(event_index));
}

/// \brief Set `event_data` for any event `event_id`, returning a reference
/// which is valid until the next call to this method
template <bool DebugMode>
EventData const &AllowedEventCalculator<DebugMode>::set_event_data(
    EventID const &event_id) {
  Index prim_event_index = event_id.prim_event_index;
  PrimEventData const &prim_event_data =
      this->prim_event_list[prim_event_index];
  Index unitcell_index = event_id.unitcell_index;

  // set this->event_data.unitcell_index
  this->event_data.unitcell_index = unitcell_index;

  // set this->event_data.event
  set_event(this->event_data.event, prim_event_data, unitcell_index,
            event_list.occ_location,
            event_list.neighbor_index[prim_event_index],
            *event_list.supercell_nlist);

  return this->event_data;
}

template <typename EventSelectorType, bool DebugMode>
AllowedKineticEventData<EventSelectorType, DebugMode>::AllowedKineticEventData(
    std::shared_ptr<system_type> _system, EventDataOptions _options)
    : options(_options) {
  if constexpr (DebugMode) {
    Log &log = CASM::log();
    log.custom("Construct AllowedKineticEventData");
    log.indent() << "Event data and selection:" << std::endl;
    log.indent() << "- impact_table_type="
                 << (options.use_neighborlist_impact_table
                         ? std::string("\"neighborlist\"")
                         : std::string("\"relative\""))
                 << std::endl;
    log.indent() << "- event_selector_type=\""
                 << this->event_selector_type_str() << "\"" << std::endl;
    log.indent() << "- assign_allowed_events_only=" << std::boolalpha
                 << options.assign_allowed_events_only << std::endl;
    log.indent() << std::endl;
    log.end_section();
  }

  system = _system;
  if (!is_clex_data(*system, "formation_energy")) {
    throw std::runtime_error(
        "Error constructing AllowedKineticEventData: no 'formation_energy' "
        "clex.");
  }

  prim_event_list = clexmonte::make_prim_event_list(*system);
  if (prim_event_list.empty()) {
    throw std::runtime_error(
        "Error constructing AllowedKineticEventData: "
        "prim event list is empty.");
  }
  if constexpr (DebugMode) {
    Log &log = CASM::log();
    log.custom("Prim event list");
    log.indent() << qto_json(prim_event_list) << std::endl << std::endl;
  }

  prim_impact_info_list = clexmonte::make_prim_impact_info_list(
      *system, prim_event_list, {"formation_energy"});

  BasicAbnormalEventHandler encountered_abnormal_event_handling_f(
      "encountered", options.throw_if_encountered_event_is_abnormal,
      options.warn_if_encountered_event_is_abnormal,
      options.disallow_if_encountered_event_is_abnormal,
      options.n_write_if_encountered_event_is_abnormal, options.output_dir,
      options.local_corr_compare_tol);
  set_encountered_abnormal_event_handling(
      encountered_abnormal_event_handling_f);
  this->encountered_abnormal_event_handling_on =
      encountered_abnormal_event_handling_f.handling_on();

  BasicAbnormalEventHandler selected_abnormal_event_handling_f(
      "selected", options.throw_if_selected_event_is_abnormal,
      options.warn_if_selected_event_is_abnormal, false,
      options.n_write_if_selected_event_is_abnormal, options.output_dir,
      options.local_corr_compare_tol);
  set_selected_abnormal_event_handling(selected_abnormal_event_handling_f);
  this->selected_abnormal_event_handling_on =
      selected_abnormal_event_handling_f.handling_on();

  /// Construct the state graph if requested
  if (_options.state_graph_options.has_value()) {
    this->state_graph = std::make_shared<state_graph::StateGraph>(
        *_options.state_graph_options);
  }

  if constexpr (DebugMode) {
    Log &log = CASM::log();
    log.indent() << "Construct AllowedKineticEventData: DONE" << std::endl
                 << std::endl;
  }
}

/// \brief Update for given state, conditions, occupants, event filters
///
/// Notes:
/// - This constructs the complete event list and impact table, and constructs
///   the event selector, which calculates all event rates.
/// - Event filters are ignored (with a warning). This is a TODO feature.
/// - Resets the `n_encountered_abnormal` and `n_selected_abnormal`
///   counters.
template <typename EventSelectorType, bool DebugMode>
void AllowedKineticEventData<EventSelectorType, DebugMode>::update(
    std::shared_ptr<StateData> _state_data,
    std::optional<std::vector<EventFilterGroup>> _event_filters,
    std::shared_ptr<engine_type> engine) {
  random_generator =
      std::make_shared<lotto::RandomGeneratorT<engine_type>>(engine);
  state_data = _state_data;

  // Warning if event_filters:
  if (_event_filters.has_value()) {
    std::cerr << "#############################################" << std::endl;
    std::cerr << "Warning: Event filters are being ignored. Use" << std::endl;
    std::cerr << "the \"high_memory\" event data type to apply " << std::endl;
    std::cerr << "event filters.                               " << std::endl;
    std::cerr << "#############################################" << std::endl;
  }

  // Current state info
  state_type const &state = *state_data->state;
  monte::OccLocation const &occ_location = *state_data->occ_location;

  if constexpr (DebugMode) {
    Log &log = CASM::log();
    log.custom("Monte Carlo State");
    log.indent() << qto_json(state) << std::endl << std::endl;
  }

  // These are constructed/re-constructed so cluster expansions point
  // at the current state
  prim_event_calculators.clear();
  for (auto const &prim_event_data : prim_event_list) {
    prim_event_calculators.emplace_back(system,
                                        prim_event_data.event_type_name);
    prim_event_calculators.back().set(&state);

    // Set a custom event state calculation function if it exists:
    auto it =
        custom_event_state_calculation_f.find(prim_event_data.event_type_name);
    if (it != custom_event_state_calculation_f.end()) {
      prim_event_calculators.back().set_custom_event_state_calculation(
          it->second);
    }
  }

  // Construct AllowedEventList
  event_list = std::make_shared<clexmonte::AllowedEventList>(
      prim_event_list, prim_impact_info_list, get_dof_values(state),
      occ_location, get_prim_neighbor_list(*system),
      get_supercell_neighbor_list(*system, state), options.use_map_index,
      options.use_neighborlist_impact_table,
      options.assign_allowed_events_only);

  if constexpr (DebugMode) {
    Log &log = CASM::log();
    log.custom("Event list summary");
    log.indent() << "- Event list container size: "
                 << event_list->allowed_event_map.n_total() << std::endl;
    log.indent() << "- Number of events: "
                 << event_list->allowed_event_map.n_assigned() << std::endl;
    log << std::endl;
    log.end_section();
  }

  // Reset "not normal" event counters
  n_encountered_abnormal.clear();
  n_selected_abnormal.clear();

  // Construct AllowedEventCalculator
  if (encountered_abnormal_event_handling_on == true &&
      encountered_abnormal_event_handling_f == nullptr) {
    throw std::runtime_error(
        "Error in AllowedKineticEventData::update: "
        "encountered_abnormal_event_handling_on == true && "
        "encountered_abnormal_event_handling_f == nullptr");
  }
  if (selected_abnormal_event_handling_on == true &&
      selected_abnormal_event_handling_f == nullptr) {
    throw std::runtime_error(
        "Error in AllowedKineticEventData::update: "
        "selected_abnormal_event_handling_on == true && "
        "selected_abnormal_event_handling_f == nullptr");
  }
  event_calculator = std::make_shared<AllowedEventCalculator<DebugMode>>(
      prim_event_list, prim_event_calculators, *event_list,
      encountered_abnormal_event_handling_on,
      encountered_abnormal_event_handling_f, n_encountered_abnormal);

  // Make event selector
  // - This calculates all rates at construction
  this->make_event_selector();
  this->event_list->allowed_event_map.clear_has_new_events();
}

template <typename EventSelectorType, bool DebugMode>
void AllowedKineticEventData<EventSelectorType, DebugMode>::run(
    state_type &state, monte::OccLocation &occ_location,
    kmc_data_type &kmc_data, SelectedEvent &selected_event,
    std::optional<monte::SelectedEventDataCollector> &collector,
    run_manager_type &run_manager,
    std::shared_ptr<occ_events::OccSystem> event_system) {
  // Function to set selected event
  bool requires_event_state = check_requires_event_state<DebugMode>(
      collector, this->selected_abnormal_event_handling_on);
  auto set_selected_event_f = [=](SelectedEvent &selected_event) {
    this->select_event(selected_event, requires_event_state);
  };

  auto set_impacted_events_f = [=](SelectedEvent &selected_event) {
    // Set impacted events
    this->event_selector->set_impacted_events(selected_event.event_index);
  };

  // Run Kinetic Monte Carlo at a single condition
  kinetic_monte_carlo_v2<DebugMode>(
      state, occ_location, kmc_data, selected_event, set_selected_event_f,
      set_impacted_events_f, collector, run_manager, event_system);
}

// -- EventSelectorType specializations --
namespace {

/// \brief Template class to specialize the implementation of the
///     `make_event_selector` function and `type_str` function
template <typename EventSelectorType, bool DebugMode>
struct event_selector_impl;

/// "sum_tree" event selector
template <bool DebugMode>
struct event_selector_impl<
    sum_tree_event_selector_type<AllowedEventCalculator<DebugMode>>,
    DebugMode> {
  typedef default_engine_type engine_type;
  typedef AllowedEventCalculator<DebugMode> event_calculator_type;
  typedef sum_tree_event_selector_type<event_calculator_type>
      event_selector_type;

  /// \brief Return "sum_tree"
  static std::string type_str() { return "sum_tree"; }

  /// \brief Construct the "sum_tree" event selector
  static std::shared_ptr<event_selector_type> make_event_selector(
      std::shared_ptr<event_calculator_type> event_calculator,
      std::shared_ptr<AllowedEventList> event_list,
      std::shared_ptr<lotto::RandomGeneratorT<engine_type>> random_generator) {
    return std::make_shared<event_selector_type>(
        event_calculator, event_list->allowed_event_map.event_index_list(),
        GetImpactFromAllowedEventList(event_list), random_generator);
  }
};

/// "vector_sum_tree" event selector
template <bool DebugMode>
struct event_selector_impl<
    vector_sum_tree_event_selector_type<AllowedEventCalculator<DebugMode>>,
    DebugMode> {
  typedef default_engine_type engine_type;
  typedef AllowedEventCalculator<DebugMode> event_calculator_type;
  typedef vector_sum_tree_event_selector_type<event_calculator_type>
      event_selector_type;

  /// \brief Return "vector_sum_tree"
  static std::string type_str() { return "vector_sum_tree"; }

  /// \brief Construct the "sum_tree" event selector
  static std::shared_ptr<event_selector_type> make_event_selector(
      std::shared_ptr<event_calculator_type> event_calculator,
      std::shared_ptr<AllowedEventList> event_list,
      std::shared_ptr<lotto::RandomGeneratorT<engine_type>> random_generator) {
    return std::make_shared<event_selector_type>(
        event_calculator, event_list->allowed_event_map.events().size(),
        GetImpactFromAllowedEventList(event_list), random_generator);
  }
};

/// "direct_sum" event selector
template <bool DebugMode>
struct event_selector_impl<
    direct_sum_event_selector_type<AllowedEventCalculator<DebugMode>>,
    DebugMode> {
  typedef default_engine_type engine_type;
  typedef AllowedEventCalculator<DebugMode> event_calculator_type;
  typedef direct_sum_event_selector_type<event_calculator_type>
      event_selector_type;

  /// \brief Return "direct_sum"
  static std::string type_str() { return "direct_sum"; }

  /// \brief Construct the "direct_sum" event selector
  static std::shared_ptr<event_selector_type> make_event_selector(
      std::shared_ptr<event_calculator_type> event_calculator,
      std::shared_ptr<AllowedEventList> event_list,
      std::shared_ptr<lotto::RandomGeneratorT<engine_type>> random_generator) {
    return std::make_shared<event_selector_type>(
        event_calculator, event_list->allowed_event_map.events().size(),
        GetImpactFromAllowedEventList(event_list), random_generator);
  }
};

template <typename EventSelectorType, typename EngineType, bool DebugMode>
std::shared_ptr<EventSelectorType> make_event_selector_impl(
    std::shared_ptr<AllowedEventCalculator<DebugMode>> event_calculator,
    std::shared_ptr<AllowedEventList> event_list,
    std::shared_ptr<lotto::RandomGeneratorT<EngineType>> random_generator) {
  return event_selector_impl<EventSelectorType, DebugMode>::make_event_selector(
      event_calculator, event_list, random_generator);
}

}  // namespace
// -- end EventSelectorType specializations --

template <typename EventSelectorType, bool DebugMode>
std::string AllowedKineticEventData<
    EventSelectorType, DebugMode>::event_selector_type_str() const {
  return event_selector_impl<EventSelectorType, DebugMode>::type_str();
}

/// \brief Constructs `event_selector`; must be called after `update`
template <typename EventSelectorType, bool DebugMode>
void AllowedKineticEventData<EventSelectorType,
                             DebugMode>::make_event_selector() {
  if constexpr (DebugMode) {
    Log &log = CASM::log();
    log.custom("Make event selector");
    log.indent() << "- event_selector_type=\""
                 << this->event_selector_type_str() << "\"" << std::endl;
    log.indent() << "- Event list container size: "
                 << event_list->allowed_event_map.n_total() << std::endl;
    log.indent() << "- Number of events: "
                 << event_list->allowed_event_map.n_assigned() << std::endl;
    log.indent() << "- Constructing event selector..." << std::endl;
  }

  // Make event selector
  // - This calculates all rates at construction
  event_selector =
      make_event_selector_impl<EventSelectorType, engine_type, DebugMode>(
          event_calculator, event_list, random_generator);

  if constexpr (DebugMode) {
    Log &log = CASM::log();
    log.indent() << "- Constructing event selector... DONE" << std::endl;
    log.indent() << "- total_rate=" << this->event_selector->total_rate()
                 << std::endl;
    log << std::endl;
    log.end_section();
  }
}

/// \brief Update for given state, conditions, occupants, event filters
template <typename EventSelectorType, bool DebugMode>
void AllowedKineticEventData<EventSelectorType, DebugMode>::select_event(
    SelectedEvent &selected_event, bool requires_event_state) {
  // If updating the event list with impacted events after the previous step
  // caused the event list to increase in size, then it needs to be
  // re-constructed.
  if (this->event_list->allowed_event_map.has_new_events()) {
    if constexpr (DebugMode) {
      Log &log = CASM::log();
      log.custom("Select event requires re-constructing event selector");
      log << std::endl;
      CASM::log().increase_indent();
    }

    this->make_event_selector();
    this->event_list->allowed_event_map.clear_has_new_events();

    if constexpr (DebugMode) {
      CASM::log().decrease_indent();
    }
  }

  if constexpr (DebugMode) {
    Log &log = CASM::log();
    log.custom("Event list summary");
    log.indent() << "- Event list container size: "
                 << event_list->allowed_event_map.n_total() << std::endl;
    log.indent() << "- Number of events: "
                 << event_list->allowed_event_map.n_assigned() << std::endl;
    log << std::endl;
    log.end_section();
  }

  // The function `only_select_event` does the following:
  // - Updates rates of events impacted by the *last* selected event (if there
  //   was a previous selection)
  // - Updates the total rate
  // - If saving states, saves the state
  // - Chooses an event and time increment
  //
  // It does not apply the event or set the impacted events.
  Index selected_event_index;

  if (!this->state_graph) {
    std::tie(selected_event_index, selected_event.time_increment) =
        event_selector->only_select_event();
  } else {
    // update impacted event rates
    event_selector->update_impacted_event_rates();

    // save state

    // determine next exit state & event

    // set state to exit state (TBD: checking if time is past next sample time)
  }
  selected_event.total_rate = event_selector->total_rate();

  EventID const &event_id =
      this->event_list->allowed_event_map.event_id(selected_event_index);
  EventData const &event_data =
      event_calculator->set_event_data(selected_event_index);
  PrimEventData const &prim_event_data =
      prim_event_list[event_id.prim_event_index];

  selected_event.event_id = event_id;
  selected_event.event_index = selected_event_index;
  selected_event.event_data = &event_data;
  selected_event.prim_event_data = &prim_event_data;

  if constexpr (DebugMode) {
    Log &log = CASM::log();
    log.custom("Selected event");

    // get ijk
    auto const &unitcell_index_converter =
        state_data->occ_location->convert().unitcell_index_converter();
    auto ijk = unitcell_index_converter(event_id.unitcell_index);
    jsonParser ijk_json;
    to_json(ijk, ijk_json, jsonParser::as_array());

    // format output
    {
      jsonParser tjson;
      to_json(selected_event, tjson, *get_event_system(*system));

      // Add event position
      jsonParser ajson;
      ajson["unitcell_index"] = event_data.unitcell_index;
      ajson["unitcell_ijk"] = ijk_json;

      ajson["event_sites_relative"] = jsonParser::array();
      for (auto const &site : prim_event_data.sites) {
        ajson["event_sites_relative"].push_back(qto_json(site));
      }

      ajson["event_sites_absolute"] = jsonParser::array();
      for (auto const &site : prim_event_data.sites) {
        ajson["event_sites_absolute"].push_back(qto_json(site + ijk));
      }
      tjson["event_position"] = ajson;

      // Add Monte Carlo state
      tjson["state"] = *state_data->state;

      log << tjson << std::endl << std::endl;
    }
  }

  if (requires_event_state) {
    if constexpr (DebugMode) {
      Log &log = CASM::log();
      log.indent() << "- Selected event state calculation required=true"
                   << std::endl;
      log.indent() << "- Event state calculation..." << std::endl;
    }

    EventStateCalculator &prim_event_calculator =
        prim_event_calculators.at(event_id.prim_event_index);
    prim_event_calculator.calculate_event_state(
        m_event_state, event_data.unitcell_index,
        event_data.event.linear_site_index, prim_event_data);
    selected_event.event_state = &m_event_state;

    if constexpr (DebugMode) {
      Log &log = CASM::log();
      log.indent() << "- Event state calculation... DONE" << std::endl
                   << std::endl;

      jsonParser event_json;
      to_json(m_event_state, event_json["event_state"]);
      event_json["unitcell_index"] = event_data.unitcell_index;
      event_json["linear_site_index"] = event_data.event.linear_site_index;
      to_json(prim_event_data, event_json["prim_event_data"]);
      log << event_json << std::endl << std::endl;
    }

    if (selected_abnormal_event_handling_on && !m_event_state.is_normal) {
      if constexpr (DebugMode) {
        Log &log = CASM::log();
        log.custom("Handle selected abnormal event...");
        log.indent() << "- event_type_name=" << prim_event_data.event_type_name
                     << std::endl;
        log.indent() << "Handling selected abnormal event ..." << std::endl;
      }
      Index &n = n_selected_abnormal[prim_event_data.event_type_name];
      n += 1;
      selected_abnormal_event_handling_f(n, m_event_state, event_data,
                                         prim_event_data,
                                         *prim_event_calculator.state());

      if constexpr (DebugMode) {
        Log &log = CASM::log();
        log.indent() << "Handling selected abnormal event... DONE" << std::endl;
        log.end_section();
      }
    }
  } else {
    if constexpr (DebugMode) {
      Log &log = CASM::log();
      log.indent() << "- Selected event state calculation required=false"
                   << std::endl;
    }
  }

  if constexpr (DebugMode) {
    Log &log = CASM::log();
    log.end_section();
  }
}

// Explicit instantiation:

// DebugMode=false
template class CompleteKineticEventData<false>;
template class AllowedKineticEventData<
    vector_sum_tree_event_selector_type<AllowedEventCalculator<false>>, false>;
template class AllowedKineticEventData<
    sum_tree_event_selector_type<AllowedEventCalculator<false>>, false>;
template class AllowedKineticEventData<
    direct_sum_event_selector_type<AllowedEventCalculator<false>>, false>;

// DebugMode=true
template class CompleteKineticEventData<true>;
template class AllowedKineticEventData<
    vector_sum_tree_event_selector_type<AllowedEventCalculator<true>>, true>;
template class AllowedKineticEventData<
    sum_tree_event_selector_type<AllowedEventCalculator<true>>, true>;
template class AllowedKineticEventData<
    direct_sum_event_selector_type<AllowedEventCalculator<true>>, true>;

}  // namespace kinetic_2
}  // namespace clexmonte
}  // namespace CASM
