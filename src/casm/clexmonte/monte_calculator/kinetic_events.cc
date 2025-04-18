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
    //    this->state_graph = std::make_shared<state_graph::StateGraph>(
    //        *_options.state_graph_options);
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
    SelectedEvent &selected_event,
    std::optional<monte::SelectedEventDataCollector> &collector,
    run_manager_type &run_manager, std::shared_ptr<kmc_data_type> _kmc_data,
    std::shared_ptr<occ_events::OccSystem> event_system) {
  if (_kmc_data == nullptr) {
    throw std::runtime_error(
        "Error in CompleteKineticEventData::run: _kmc_data==nullptr");
  }
  this->kmc_data = _kmc_data;
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
      state, occ_location, *_kmc_data, selected_event, set_selected_event_f,
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
  selected_event.time = this->kmc_data->time + selected_event.time_increment;
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
    AllowedEventList &_event_list,
    std::vector<std::shared_ptr<event_group::EventGroup<DebugMode>>>
        &_event_group,
    bool _abnormal_event_handling_on,
    AbnormalEventHandlingFunction &_handling_f,
    std::map<std::string, Index> &_n_encountered_abnormal)
    : prim_event_list(_prim_event_list),
      prim_event_calculators(_prim_event_calculators),
      event_list(_event_list),
      event_group(_event_group),
      abnormal_event_handling_on(_abnormal_event_handling_on),
      handling_f(_handling_f),
      n_encountered_abnormal(_n_encountered_abnormal) {}

/// \brief Update `event_state` for event `event_index` in the current state
/// and return the event rate; if the event is no longer allowed and not
/// included in any event group, free the event.
///
/// Calculate event rates:
/// - If event is ungrouped (group=0):
///   - add rate to event selector
/// - If event is grouped (group!=0),
///   but index_in_group is not set (index_in_group=-1):
///   - add ref_rate and ref_dE_activated to EventGroup,
///   - set index_in_group
///   - add rate and dE_activated to group's current state
///   - set 0.0 rate in event selector
/// - If event is grouped (group!=0),
///   and index_in_group is set (index_in_group!=-1):
///   - set rate and dE_activated to group's current state
///   - set 0.0 rate in event selector
///
/// \param event_index Linear index of event in event_list.allowed_event_map
///
/// \return rate, The rate of the event in the current state
template <bool DebugMode>
double AllowedEventCalculator<DebugMode>::calculate_rate(Index event_index) {
  AllowedEventData const &allowed_event_data =
      event_list.allowed_event_map.events()[event_index];

  auto &log = CASM::log();
  log.increase_indent();
  if (!allowed_event_data.is_assigned) {
    log.indent() << "- calculate " << event_index << " / {"
                 << allowed_event_data.event_id.unitcell_index << ","
                 << allowed_event_data.event_id.prim_event_index
                 << "}: (not assigned) " << std::endl;

    event_state.is_allowed = false;
    event_state.rate = 0.0;

    // this value goes into the ungrouped event selector:
    log.decrease_indent();
    return event_state.rate;
  } else {
    this->calculate_rate(allowed_event_data.event_id);

    if (allowed_event_data.group == 0) {
      // Ungrouped events:
      // - If the event is not allowed, we need to free the event from the
      //   AllowedEventMap

      if (!event_state.is_allowed) {
        log.indent() << "- calculate " << event_index << " / {"
                     << allowed_event_data.event_id.unitcell_index << ","
                     << allowed_event_data.event_id.prim_event_index
                     << "}: (group 0, not allowed) rate=" << event_state.rate
                     << std::endl;
        event_list.allowed_event_map.free(allowed_event_data.event_id);
      } else {
        log.indent() << "- calculate " << event_index << " / {"
                     << allowed_event_data.event_id.unitcell_index << ","
                     << allowed_event_data.event_id.prim_event_index
                     << "}: (group 0, is allowed) rate=" << event_state.rate
                     << std::endl;
      }

      // This value goes into the ungrouped event selector:
      log.decrease_indent();
      return event_state.rate;
    } else {
      // Grouped events:
      // - If the event is not allowed, we still have to keep the event in the
      //   AllowedEventMap
      // - If the event is newly added to a group (index_in_group == -1), we
      //   need to store the reference rate and dE_activated in the group
      // - If the event is already in a group (index_in_group != -1), we need to
      //   update the current state rate and dE_activated in the group

      auto &group = *this->event_group[allowed_event_data.group];
      if (allowed_event_data.index_in_group == -1) {
        Index index_in_group = group.event.size();
        log.indent() << "- calculate " << event_index << " / {"
                     << allowed_event_data.event_id.unitcell_index << ","
                     << allowed_event_data.event_id.prim_event_index
                     << "}: (group " << allowed_event_data.group
                     << ", index_in_group " << index_in_group
                     << ", new to group) rate=" << event_state.rate
                     << std::endl;

        group.event.emplace_back(event_index, event_state.rate,
                                 event_state.dE_activated);
        event_list.allowed_event_map.set_event_index_in_group(event_index,
                                                              index_in_group);
        group.new_state.rate.push_back(event_state.rate);
        group.new_state.dE_activated.push_back(event_state.dE_activated);
      } else {
        log.indent() << "- calculate " << event_index << " / {"
                     << allowed_event_data.event_id.unitcell_index << ","
                     << allowed_event_data.event_id.prim_event_index
                     << "}: (group " << allowed_event_data.group
                     << ", index_in_group " << allowed_event_data.index_in_group
                     << ", existing in group) rate=" << event_state.rate
                     << std::endl;

        group.new_state.rate[allowed_event_data.index_in_group] =
            event_state.rate;
        group.new_state.dE_activated[allowed_event_data.index_in_group] =
            event_state.dE_activated;
      }

      // This value goes into the ungrouped event selector:
      log.decrease_indent();
      return 0.0;
    }
  }
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

  /// Set use_event_groups:
  this->use_event_groups = options.state_graph_options.has_value();

  /// Construct the state graph if requested
  if (_options.state_graph_options.has_value()) {
    //    this->state_graph = std::make_shared<state_graph::StateGraph>(
    //        *_options.state_graph_options);
    //    if constexpr (DebugMode) {
    //      Log &log = CASM::log();
    //      log.custom("Construct state graph");
    //      log.indent() << "- n_recent_events = "
    //                   << this->state_graph->opt.n_recent_events << std::endl;
    //      log.indent() << "- n_state = " << this->state_graph->opt.n_states
    //                   << std::endl
    //                   << std::endl;
    //    }
  } else {
    if constexpr (DebugMode) {
      Log &log = CASM::log();
      log.custom("Construct state graph");
      log.indent() << "- No state graph" << std::endl << std::endl;
    }
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
      prim_event_list, prim_event_calculators, *event_list, event_group,
      encountered_abnormal_event_handling_on,
      encountered_abnormal_event_handling_f, n_encountered_abnormal);

  // Make event selector
  // - This calculates all rates at construction
  this->make_event_selector();
}

template <typename EventSelectorType, bool DebugMode>
void AllowedKineticEventData<EventSelectorType, DebugMode>::run(
    state_type &state, monte::OccLocation &occ_location,
    SelectedEvent &selected_event,
    std::optional<monte::SelectedEventDataCollector> &collector,
    run_manager_type &run_manager, std::shared_ptr<kmc_data_type> _kmc_data,
    std::shared_ptr<occ_events::OccSystem> event_system) {
  if (_kmc_data == nullptr) {
    throw std::runtime_error(
        "Error in AllowedKineticEventData::run: _kmc_data==nullptr");
  }
  this->kmc_data = _kmc_data;

  // Function to set selected event
  bool requires_event_state = check_requires_event_state<DebugMode>(
      collector, this->selected_abnormal_event_handling_on);
  auto set_selected_event_f = [=](SelectedEvent &selected_event) {
    this->select_event(selected_event, requires_event_state);
  };

  // Function to set impacted events and handle the consequences of applying
  // the last selected event
  auto set_impacted_events_f = [=](SelectedEvent &selected_event) {
    this->set_impacted_events(selected_event);
  };

  // Run Kinetic Monte Carlo at a single condition
  this->kmc_data = _kmc_data;
  kinetic_monte_carlo_v2<DebugMode>(
      state, occ_location, *_kmc_data, selected_event, set_selected_event_f,
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

/// \brief Constructs `event_selector` from the current `event_calculator`,
///     `event_list`, and `random_generator`
///
/// - This is called by `update`
/// - This should be called if the AllowedEventMap is resized in order to
///   reconstruct the event selector
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

  this->event_list->allowed_event_map.clear_has_been_resized();

  // If using event groups:
  // - If no default group, create it
  // - build default group and select first event / time
  if (this->use_event_groups && this->event_group.size() == 0) {
    Index group = 0;
    monte::TimeType group_time = 0.0;
    if (this->kmc_data != nullptr) {
      group_time = this->kmc_data->time;
    }
    this->event_group.emplace_back(
        std::make_shared<event_group::EventGroup<DebugMode>>(
            this->state_data, group_time, group));
    this->current_groups.insert(group);
    this->select_next_event_for(this->current_groups);
  }
}

/// \brief Reconstruct the event selector if updating the allowed event list
///     caused it to increase in size
template <typename EventSelectorType, bool DebugMode>
void AllowedKineticEventData<EventSelectorType,
                             DebugMode>::make_event_selector_if_resized() {
  if (this->event_list->allowed_event_map.has_been_resized()) {
    if constexpr (DebugMode) {
      Log &log = CASM::log();
      log.custom("Reconstructing the event selector");
      log << std::endl;
      CASM::log().increase_indent();
    }

    this->make_event_selector();

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
}

template <typename EventSelectorType, bool DebugMode>
void AllowedKineticEventData<EventSelectorType, DebugMode>::set_impacted_events(
    SelectedEvent &selected_event) {
  Log &log = CASM::log();

  // Set impacted events
  this->event_selector->set_impacted_events(selected_event.event_index);

  if (!this->use_event_groups) {
    log.indent() << "- use_event_groups=false" << std::endl;
    this->make_event_selector_if_resized();
    this->event_selector->update_impacted_event_rates();
    this->event_selector->clear_impacted_events();
  } else {
    log.indent() << "- use_event_groups=true" << std::endl;
    this->regroup_impacted_events();
  }
}

/// \brief Update for given state, conditions, occupants, event filters
///
/// - If there are no event groups, just select from the event selector
/// - If there are event groups, determine which event group moves next
template <typename EventSelectorType, bool DebugMode>
void AllowedKineticEventData<EventSelectorType, DebugMode>::select_event(
    SelectedEvent &selected_event, bool requires_event_state) {
  Index selected_event_index;

  if (!this->use_event_groups) {
    //    Log &log = CASM::log();
    //    log.indent() << "- use_event_groups=false" << std::endl;
    if (this->kmc_data == nullptr) {
      throw std::runtime_error(
          "Error in AllowedKineticEventData::select_event: "
          "Cannot select event, kmc_data==nullptr");
    }

    // The function `only_select_event` does the following:
    // - Updates rates of events impacted by the *last* selected event if that
    //   hasn't been done yet (by checking if the impacted events ptr is set)
    // - Updates the total rate
    // - Chooses an event and time increment
    //
    // It does not apply the event or set the impacted events.

    std::tie(selected_event_index, selected_event.time_increment) =
        event_selector->only_select_event();
    selected_event.time = this->kmc_data->time + selected_event.time_increment;
    selected_event.total_rate = event_selector->total_rate();

  } else {
    Log &log = CASM::log();
    log.custom("Select event, using event groups");
    log.indent() << "- use_event_groups=true..." << std::endl;

    // Set `next_event_group` by finding which group moves next
    this->set_next_event_group();

    // Set selected event
    selected_event.group = this->next_event_group;
    log.indent() << "- next_event_group=" << this->next_event_group
                 << std::endl;
    auto const &g = *this->event_group[this->next_event_group];
    log.indent() << "- g.next_event_index=" << g.next_event_index << std::endl;
    log.indent() << "- g.next_time=" << g.next_time << std::endl;
    log.indent() << "- g.next_time_increment=" << g.next_time_increment
                 << std::endl;
    log.indent() << "- g.next_state=" << g.next_state << std::endl;
    selected_event_index = g.next_event_index;
    selected_event.time = g.next_time;
    selected_event.time_increment = g.next_time_increment;
    selected_event.group_state = g.next_state;

    if (selected_event_index == -1) {
      throw std::runtime_error(
          "Error in AllowedKineticEventData::select_event: "
          "selected_event_index == -1");
    }
    if (selected_event_index >= this->event_list->allowed_event_map.n_total()) {
      throw std::runtime_error(
          "Error in AllowedKineticEventData::select_event: "
          "selected_event_index >= n_total");
    }
    if (!this->event_list->allowed_event_map.event_data(selected_event_index)
             .is_assigned) {
      throw std::runtime_error(
          "Error in AllowedKineticEventData::select_event: "
          "selected_event_index is not assigned");
    }

    log.indent() << "- selected_group=" << this->next_event_group << std::endl;
    log.indent() << "- selected_event_index=" << selected_event_index
                 << std::endl;
    log.indent() << "- next_state=" << selected_event.group_state << std::endl;
    log.indent() << "- time=" << selected_event.time << std::endl;
    log.indent() << "- time_increment=" << selected_event.time_increment
                 << std::endl;
    log.indent() << "- find next event overall... DONE" << std::endl;
    log << std::endl;
    log.end_section();
  }

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

  {
    Log &log = CASM::log();
    log.indent() << "- selected={" << prim_event_data.event_type_name << ","
                 << prim_event_data.equivalent_index << ","
                 << prim_event_data.is_forward << "," << event_id.unitcell_index
                 << "} / {" << event_id.unitcell_index << ","
                 << event_id.prim_event_index << "}" << std::endl
                 << std::endl;
    log.end_section();
  }

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

/// \brief Find and regroup impacted events
///
/// - Checks if `event_selector` has impacted events
/// - Evolve impacted EventGroup to the current time
/// - Regroup impacted events
/// - Set the next event and time for the impacted & new EventGroup
/// - Determine the next event overall
///
template <typename EventSelectorType, bool DebugMode>
void AllowedKineticEventData<EventSelectorType,
                             DebugMode>::regroup_impacted_events() {
  Log &log = CASM::log();

  log.indent() << "- regrouping impacted events..." << std::endl;
  log.indent() << "- checking event impact..." << std::endl;
  if (!this->event_selector->has_impacted_events()) {
    //    log.indent() << "- has_impacted_events=false" << std::endl;
    log.indent() << "- regrouping impacted events... DONE" << std::endl;
    return;
  }

  // Groups (including "ungrouped" group 0) impacted by the last event:
  static std::set<Index> impacted_groups;
  impacted_groups.clear();

  // Groups (excluding "ungrouped" group 0) impacted by the last event:
  static std::set<Index> impacted_groups_excluding_group_zero;
  impacted_groups_excluding_group_zero.clear();

  // Groups to update:
  // - Groups that need the next selected event to be chosen
  static std::set<Index> new_and_modified_groups;
  new_and_modified_groups.clear();

  // Groups that are merged into another group and need to be deleted:
  static std::set<Index> groups_to_delete;
  groups_to_delete.clear();

  // Check last event impact
  //  log.indent() << "- has_impacted_events=true" << std::endl;
  auto &allowed_event_map = this->event_list->allowed_event_map;
  for (Index event_index : this->event_selector->get_impacted_events()) {
    Index group = allowed_event_map.event_group(event_index);
    impacted_groups.insert(group);
    if (group != 0) {
      impacted_groups_excluding_group_zero.insert(group);
    }
  }
  //  log.indent() << "- impacted_groups.size=" << impacted_groups.size()
  //               << std::endl;
  //  log.indent() << "- impacted_groups=" << qto_json(impacted_groups)
  //               << std::endl;

  // Resolve which state impacted groups are in at the time the event
  if (this->kmc_data == nullptr) {
    throw std::runtime_error(
        "Error in AllowedKineticEventData::regroup_impacted_events: "
        "Cannot evolve impacted groups to the current time, "
        "kmc_data==nullptr");
  }
  this->resolve_state_for(this->kmc_data->time, impacted_groups);

  // Regroup impacted events if necessary

  // ----------------------------------- //
  // case 1: 0 existing groups (excluding group 0) impacted
  // - create a new group (up to max group size)
  // - add all to new group, set index_in_group to -1
  if (impacted_groups_excluding_group_zero.size() == 0) {
    //    log.indent() << "- 0 existing groups impacted" << std::endl;
    if (this->n_groups() >= 10) {
      //      log.indent() << "- max group size reached" << std::endl;
    } else {
      //      log.indent() << "- max group size not reached" << std::endl;
      Index new_group = this->add_group();
      this->add_events_to_group(this->event_selector->get_impacted_events(),
                                new_group);
      new_and_modified_groups = impacted_groups;
      new_and_modified_groups.insert(new_group);
    }

  }
  // ----------------------------------- //
  // case 2: 1 existing group impacted
  // - add ungrouped impacted events to existing group
  else if (impacted_groups_excluding_group_zero.size() == 1) {
    Index existing_group = *impacted_groups_excluding_group_zero.begin();

    //    log.indent() << "- 1 existing group impacted (group " <<
    //    existing_group
    //                 << ")" << std::endl;
    //    log.increase_indent();
    for (Index event_index : this->event_selector->get_impacted_events()) {
      Index group = allowed_event_map.event_group(event_index);
      if (group == 0) {
        //        log.indent() << "- move event " << event_index << " from group
        //        "
        //                     << group << " to group " << existing_group <<
        //                     std::endl;
        allowed_event_map.set_event_group(event_index, existing_group);
        allowed_event_map.set_event_index_in_group(event_index, -1);
      }
    }
    log.decrease_indent();
    new_and_modified_groups = impacted_groups;
  }
  // ----------------------------------- //
  // case 3: >1 existing group impacted
  // - merge existing groups
  // - pick one existing group to keep,
  // - add all impacted events to the existing group
  // - schedule other existing groups for deletion
  else {
    //    log.indent() << "- >1 existing group impacted" << std::endl;
    Index existing_group = *impacted_groups_excluding_group_zero.begin();
    //    log.indent() << "- keeping group " << existing_group << std::endl;
    //    log.increase_indent();
    for (Index event_index : this->event_selector->get_impacted_events()) {
      Index group = allowed_event_map.event_group(event_index);
      if (group != existing_group) {
        //        log.indent() << "- move event " << event_index << " from group
        //        "
        //                     << group << " to group " << existing_group <<
        //                     std::endl;
        allowed_event_map.set_event_group(event_index, existing_group);
        allowed_event_map.set_event_index_in_group(event_index, -1);
      }
    }
    //    log.decrease_indent();
    for (Index group : impacted_groups) {
      if (group == existing_group) {
        new_and_modified_groups.insert(existing_group);
      } else if (group != 0) {
        groups_to_delete.insert(group);
      }
    }
  }
  // ----------------------------------- //

  log.indent() << "- new_and_modified_groups="
               << qto_json(new_and_modified_groups) << std::endl;
  for (Index group : new_and_modified_groups) {
    log.increase_indent();
    log.indent() << "- group=" << group
                 << " #events=" << this->event_group[group]->event.size()
                 << std::endl;
    log.decrease_indent();
  }
  log.indent() << "- groups_to_delete=" << qto_json(groups_to_delete)
               << std::endl;
  log.indent() << "- checking event impact... DONE" << std::endl;
  log << std::endl;

  this->make_event_selector_if_resized();
  this->event_selector->update_impacted_event_rates();
  this->event_selector->clear_impacted_events();

  // If saving states, saves state
  {
    //    log.indent() << "- saving states..." << std::endl;
    //    log.increase_indent();
    for (Index group : new_and_modified_groups) {
      auto &g = *this->event_group[group];
      if (group != 0) {
        //        log.indent() << "- save state (group=" << group << ")" <<
        //        std::endl;
        g.save_state(*this->state_data);
      }
      //      else {
      //        log.indent() << "- do not save state (group=" << group << ")"
      //                     << std::endl;
      //      }
    }
    //    log.indent() << "- saving states... DONE" << std::endl;
    //    log << std::endl;
    //    log.decrease_indent();
  }

  //  if (selected_event.group_state != -1) {
  //    log.indent() << "- restore state..." << std::endl;
  //    // Restore selected state
  //    this->event_group[selected_event.group]->restore_state(
  //        selected_event.group_state, *this->state_data);
  //    log.indent() << "- restore state... DONE" << std::endl;
  //    log << std::endl;
  //  }

  // Select the next event for each new or modified group:
  this->select_next_event_for(new_and_modified_groups);

  // Delete groups scheduled for deletion
  this->delete_groups(groups_to_delete);

  // Finish
  log.indent() << "- regrouping impacted events... DONE" << std::endl
               << std::endl;
}

/// \brief Get the current number of groups (includes group 0)
template <typename EventSelectorType, bool DebugMode>
Index AllowedKineticEventData<EventSelectorType, DebugMode>::n_groups() {
  Index n_groups = 0;
  for (Index group = 0; group < this->event_group.size(); ++group) {
    if (this->event_group[group] != nullptr) {
      ++n_groups;
    }
  }
  return n_groups;
}

/// \brief Construct and add a new event group
///
/// \return new_group The index of the new group
template <typename EventSelectorType, bool DebugMode>
Index AllowedKineticEventData<EventSelectorType, DebugMode>::add_group() {
  Log &log = CASM::log();
  Index new_group = 1;
  while (new_group < this->event_group.size()) {
    if (this->event_group[new_group] == nullptr) {
      break;
    }
    ++new_group;
  }
  log.indent() << "- adding group " << new_group << std::endl;
  if (this->kmc_data == nullptr) {
    throw std::runtime_error(
        "Error in AllowedKineticEventData::regroup_impacted_events: "
        "Cannot evolve impacted groups to the current time, "
        "kmc_data==nullptr");
  }

  // Determine the new group's current time - Use 0.0 as a default
  monte::TimeType new_group_time = 0.0;
  if (this->kmc_data != nullptr) {
    new_group_time = this->kmc_data->time;
  }
  if (new_group == this->event_group.size()) {
    this->event_group.emplace_back(
        std::make_shared<event_group::EventGroup<DebugMode>>(
            this->state_data, new_group_time, new_group));
  } else {
    this->event_group[new_group] =
        std::make_shared<event_group::EventGroup<DebugMode>>(
            this->state_data, new_group_time, new_group);
  }
  this->current_groups.insert(new_group);
  return new_group;
}

/// \brief Add events to the specified event group
///
/// - This currently does not remove events from an existing group. It's
///   expected existing groups other than group 0 will be deleted afterwards.
///
/// \param event_indices Indices of events to move into the group
/// \param group Index of the event group
template <typename EventSelectorType, bool DebugMode>
void AllowedKineticEventData<EventSelectorType, DebugMode>::add_events_to_group(
    std::vector<Index> const &event_indices, Index group) {
  Log &log = CASM::log();
  log.increase_indent();
  AllowedEventMap &allowed_event_map = this->event_list->allowed_event_map;
  for (Index event_index : event_indices) {
    Index initial_group = allowed_event_map.event_group(event_index);
    log.indent() << "- move event " << event_index << " from group "
                 << initial_group << " to group " << group << std::endl;
    allowed_event_map.set_event_group(event_index, group);
    allowed_event_map.set_event_index_in_group(event_index, -1);
  }
  log.decrease_indent();
}

/// \brief Select the next event for specified event groups
///
/// \param groups Groups to update by selecting the next event to occur. May
///     include group 0.
template <typename EventSelectorType, bool DebugMode>
void AllowedKineticEventData<EventSelectorType, DebugMode>::
    select_next_event_for(std::set<Index> const &groups) {
  Log &log = CASM::log();
  log.indent() << "- Selecting next events..." << std::endl;
  log.increase_indent();

  // TODO: select events for each impacted / new group
  for (Index group : groups) {
    log.indent() << "- group=" << group << ":" << std::endl;
    log.increase_indent();
    auto &g = *this->event_group[group];
    if (group == 0) {
      if (this->event_selector->total_rate() == 0.0) {
        log.indent() << "- NO ALLOWED EVENTS" << std::endl;

        g.next_state = -1;
        g.next_event_index = -1;
        g.next_time_increment = std::numeric_limits<double>::max();
        g.next_time = std::numeric_limits<double>::max();
      } else {
        log.indent() << "- select an ungrouped event" << std::endl;
        g.next_state = -1;
        std::tie(g.next_event_index, g.next_time_increment) =
            this->event_selector->only_select_event();
        g.next_time = g.current_time + g.next_time_increment;
      }
      log.indent() << "- total_rate=" << this->event_selector->total_rate()
                   << std::endl;
    } else {
      g.select_next_event(*this->random_generator);
    }

    log.indent() << "- next_state=" << g.next_state << std::endl;
    log.indent() << "- next_event_index=" << g.next_event_index << std::endl;
    log.indent() << "- current_time=" << g.current_time << std::endl;
    log.indent() << "- next_time_increment=" << g.next_time_increment
                 << std::endl;
    log.indent() << "- next_time=" << g.next_time << std::endl;
    log.decrease_indent();
  }
  log.decrease_indent();
  log.indent() << "- Selecting next events... DONE" << std::endl;
  log << std::endl;
}

/// \brief Resolve which state all groups are in at the specified time,
///    under the assumption that they only transition between transient
///    states in the current chain
template <typename EventSelectorType, bool DebugMode>
void AllowedKineticEventData<EventSelectorType, DebugMode>::resolve_state(
    monte::TimeType time) {
  Log &log = CASM::log();
  log.indent() << "- resolving state of all groups..." << std::endl;
  log.increase_indent();
  for (Index group = 0; group < this->event_group.size(); ++group) {
    if (group == 0) {
      continue;
    }
    if (this->event_group[group] == nullptr) {
      log.indent() << "- group=" << group << " does not exist" << std::endl;
      continue;
    }
    auto &g = *this->event_group[group];
    log.indent() << "- group=" << group << ":" << std::endl;
    g.resolve_state(time, *this->random_generator);
  }

  log.indent() << "- resolving state of all groups... DONE" << std::endl;
  log << std::endl;
  log.decrease_indent();
}

/// \brief Resolve which state specified groups are in at the specified time,
///    under the assumption that they only transition between transient
///    states in the current chain
template <typename EventSelectorType, bool DebugMode>
void AllowedKineticEventData<EventSelectorType, DebugMode>::resolve_state_for(
    monte::TimeType time, std::set<Index> const &groups) {
  Log &log = CASM::log();
  log.indent() << "- resolving state for groups..." << std::endl;
  log.increase_indent();
  for (Index group : groups) {
    auto &g = *this->event_group[group];
    log.indent() << "- group=" << group << ":" << std::endl;
    g.resolve_state(time, *this->random_generator);
  }

  log.indent() << "- resolving state for groups... DONE" << std::endl;
  log << std::endl;
  log.decrease_indent();
}

/// \brief Delete specified event groups
///
/// \param groups Groups to delete. Groups are deleted by resetting the
///     corresponding `event_group` pointer to `nullptr`.
template <typename EventSelectorType, bool DebugMode>
void AllowedKineticEventData<EventSelectorType, DebugMode>::delete_groups(
    std::set<Index> const &groups) {
  Log &log = CASM::log();
  log.indent() << "- delete groups..." << std::endl;
  log.increase_indent();
  for (Index group : groups) {
    log.indent() << "- delete group " << group << std::endl;
    this->event_group[group].reset();
    this->current_groups.erase(group);
  }
  log.decrease_indent();
  log.indent() << "- delete groups... DONE" << std::endl;
}

/// \brief Set `next_event_group` by finding which group moves next
template <typename EventSelectorType, bool DebugMode>
void AllowedKineticEventData<EventSelectorType,
                             DebugMode>::set_next_event_group() {
  Log &log = CASM::log();
  log.indent() << "- find next event overall..." << std::endl;
  log.increase_indent();

  Index group = 0;
  this->next_event_group = group;
  double next_time;
  double min_time;
  {
    log.indent() << "- group=" << group << ":" << std::endl;
    log.increase_indent();

    next_time = this->event_group[group]->next_time;
    min_time = next_time;
    log.indent() << "- next_time=" << next_time << "(min=" << min_time
                 << ", selected_group=" << this->next_event_group << ")"
                 << std::endl;
    log.decrease_indent();
  }

  ++group;
  for (; group < this->event_group.size(); ++group) {
    log.indent() << "- group=" << group << ":" << std::endl;
    log.increase_indent();
    if (this->event_group[group] == nullptr) {
      log.indent() << "- group does not exist" << std::endl;
      log.decrease_indent();
      continue;
    }
    next_time = this->event_group[group]->next_time;
    if (next_time < min_time) {
      min_time = next_time;
      this->next_event_group = group;
    }
    log.indent() << "- next_time=" << next_time << "(min=" << min_time
                 << ", selected_group=" << this->next_event_group << ")"
                 << std::endl;
    log.decrease_indent();
  }
  log.decrease_indent();
  log.indent() << "- find next event overall... DONE" << std::endl << std::endl;
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
