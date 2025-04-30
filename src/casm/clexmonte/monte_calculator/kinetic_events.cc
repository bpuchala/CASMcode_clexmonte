#include "casm/clexmonte/monte_calculator/kinetic_events.hh"

#include "casm/casm_io/SafeOfstream.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/events/io/json/event_data_json_io.hh"
#include "casm/clexmonte/events/io/stream/EventState_stream_io.hh"
#include "casm/clexmonte/events/state_graph.hh"
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
    throw std::runtime_error(
        "Error constructing CompleteKineticEventData: state graph is not "
        "supported");
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

  auto apply_selected_event_f = [=](state_type &state,
                                    monte::OccLocation &occ_location,
                                    SelectedEvent &selected_event) {
    // Set impacted events
    occ_location.apply(selected_event.event_data->event, get_occupation(state));
    this->event_selector->set_impacted_events(selected_event.event_id);
  };

  // Run Kinetic Monte Carlo at a single condition
  kinetic_monte_carlo_v2<DebugMode>(
      state, occ_location, *_kmc_data, selected_event, set_selected_event_f,
      apply_selected_event_f, collector, run_manager, event_system);
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

// -- EventGroupManager --

/// \brief Save the occupation and location info for a state
void save(Index linear_site_index, SavedSiteData &site_data,
          StateData const &state_data) {
  Eigen::VectorXi const &occupation = get_occupation(*state_data.state);
  monte::OccLocation const &occ_location = *state_data.occ_location;

  // save occupation
  site_data.occ = occupation(linear_site_index);

  // save location info
  Index mol_id = occ_location.l_to_mol_id(linear_site_index);
  monte::Mol const &mol = occ_location.mol(mol_id);
  site_data.atom_id.clear();
  site_data.atom.clear();
  for (Index _atom_id : mol.component) {
    site_data.atom_id.push_back(_atom_id);
    site_data.atom.push_back(occ_location.atom(_atom_id));
  }
}

/// \brief Restore the occupation and location info for a state
void restore(Index linear_site_index, SavedSiteData const &site_data,
             StateData &state_data) {
  // Note: previously state_data was only planned as const access
  state_type &state = const_cast<state_type &>(*state_data.state);
  Eigen::VectorXi &occupation = get_occupation(state);
  monte::OccLocation &occ_location =
      const_cast<monte::OccLocation &>(*state_data.occ_location);

  occupation(linear_site_index) = site_data.occ;
  occ_location.replace_mol(linear_site_index, site_data.occ, site_data.atom_id,
                           site_data.atom);
}

/// \brief Check if the saved occupation matches the current occupation
///     (indistinguishable occupant comparison, not tracer / atom id equality)
bool is_equal(Index linear_site_index, SavedSiteData const &site_data,
              StateData const &state_data) {
  return site_data.occ == get_occupation(*state_data.state)(linear_site_index);
}

template <bool DebugMode>
AllowedEventCalculator<DebugMode>::AllowedEventCalculator(
    std::vector<PrimEventData> const &_prim_event_list,
    std::vector<EventStateCalculator> const &_prim_event_calculators,
    AllowedEventList &_event_list,
    std::shared_ptr<EventGroupManager<DebugMode>> _event_group_manager,
    bool _abnormal_event_handling_on,
    AbnormalEventHandlingFunction &_handling_f,
    std::map<std::string, Index> &_n_encountered_abnormal)
    : prim_event_list(_prim_event_list),
      prim_event_calculators(_prim_event_calculators),
      event_list(_event_list),
      event_group_manager(_event_group_manager),
      abnormal_event_handling_on(_abnormal_event_handling_on),
      handling_f(_handling_f),
      n_encountered_abnormal(_n_encountered_abnormal) {}

/// \brief Update `event_state` for event `event_index` in the current state
/// and return the event rate; if the event is no longer allowed and not
/// included in any event group, free the event.
///
/// Calculate event rates:
/// - If no event group manager:
///   - set rate in event selector
/// - If event group manager:
///   - Set event data in event group manager
///   - If event is ungrouped (group=0):
///     - set rate in event selector
///   - If event is grouped (group!=0):
///     - set 0.0 rate in event selector
///
/// \param event_index Linear index of event in event_list.allowed_event_map
///
/// \return rate, The rate of the event in the current state
template <bool DebugMode>
double AllowedEventCalculator<DebugMode>::calculate_rate(Index event_index) {
  AllowedEventData const &allowed_event_data =
      event_list.allowed_event_map.events()[event_index];

  // auto &log = CASM::log();
  // log.increase_indent();
  if (!allowed_event_data.is_assigned) {
    // log.indent() << "- calculate " << event_index << " / {"
    //              << allowed_event_data.event_id.unitcell_index << ","
    //              << allowed_event_data.event_id.prim_event_index
    //              << "}: (not assigned) " << std::endl;

    event_state.is_allowed = false;
    event_state.rate = 0.0;

    // this value goes into the ungrouped event selector:
    // log.decrease_indent();
    return event_state.rate;
  } else {
    this->calculate_rate(allowed_event_data.event_id);

    if (!event_state.is_allowed) {
      // log.indent() << "- calculate " << event_index << " / {"
      //              << allowed_event_data.event_id.unitcell_index << ","
      //              << allowed_event_data.event_id.prim_event_index
      //              << "}: (group 0, not allowed) rate=" << event_state.rate
      //              << std::endl;
      event_list.allowed_event_map.free(allowed_event_data.event_id);
    } else {
      // log.indent() << "- calculate " << event_index << " / {"
      //              << allowed_event_data.event_id.unitcell_index << ","
      //              << allowed_event_data.event_id.prim_event_index
      //              << "}: (group 0, is allowed) rate=" << event_state.rate
      //              << std::endl;
    }

    // This value goes into the ungrouped event selector:
    // log.decrease_indent();
    return event_state.rate;
  }
  // TODO: add event group handling
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

  /// Construct the event group manager if state graph construction is requested
  if (_options.state_graph_options.has_value()) {
    // TODO
    if constexpr (DebugMode) {
      Log &log = CASM::log();
      log.custom("Construct state graph");
      auto const &state_graph_options = *this->options.state_graph_options;
      log.indent() << "- n_recent_events = "
                   << state_graph_options.n_recent_events << std::endl;
      log.indent() << "- n_state = " << state_graph_options.n_states
                   << std::endl
                   << std::endl;
    }
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
      prim_event_list, prim_event_calculators, *event_list, event_group_manager,
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
  auto apply_selected_event_f = [&](state_type &state,
                                    monte::OccLocation &occ_location,
                                    SelectedEvent &selected_event) {
    // Log &log = CASM::log();
    if (this->event_group_manager == nullptr) {
      // log.indent() << "- use_event_groups=false" << std::endl;

      begin_section<DebugMode>("Evaluate selected event functions");
      if (collector.has_value()) {
        debug_collect<DebugMode>(collector.value());
        collector->collect();
      }
      end_section<DebugMode>();

      begin_section<DebugMode>("Apply selected event");
      occ_location.apply(selected_event.event_data->event,
                         get_occupation(state));
      end_section<DebugMode>();

      begin_section<DebugMode>("Update impacted events");
      this->event_selector->set_impacted_events(selected_event.event_index);
      this->make_event_selector_if_resized();
      this->event_selector->update_impacted_event_rates();
      this->event_selector->clear_impacted_events();
      end_section<DebugMode>();

    } else {
      // log.indent() << "- use_event_groups=true" << std::endl;

      auto &m = *this->event_group_manager;

      static std::set<Index> new_and_modified_groups;

      begin_section<DebugMode>("Restore selected event group exit state");
      //{
      //  auto &g = m.group(selected_event.group);
      //  g.restore_state(selected_event.group_state);
      //  g.last_time = selected_event.time;
      //}
      end_section<DebugMode>();

      // This section does the following:
      // - Set the impacted events & groups
      // - Evolve impacted groups to the event time
      // - Regroup impacted events as necessary
      //   - If a new group, add initial state
      //   - Add current site and event information for existing states
      // - Clear the impacted events
      begin_section<DebugMode>("Regroup impacted events");
      this->event_selector->set_impacted_events(selected_event.event_index);
      // this->regroup_impacted_events(new_and_modified_groups,
      //                               selected_event.group);
      end_section<DebugMode>();

      begin_section<DebugMode>("Evaluate selected event functions");
      if (collector.has_value()) {
        debug_collect<DebugMode>(collector.value());
        collector->collect();
      }
      end_section<DebugMode>();

      begin_section<DebugMode>("Apply selected event");
      occ_location.apply(selected_event.event_data->event,
                         get_occupation(state));
      end_section<DebugMode>();

      // This section does the following:
      // - Copies the current state to the next available state for each group
      // - Updates impacted event rates
      // - Saves states
      // - Select the next event and time for each group
      begin_section<DebugMode>("Update impacted event groups");
      // this->update_impacted_event_groups(new_and_modified_groups);
      end_section<DebugMode>();
    }
  };

  // Run Kinetic Monte Carlo at a single condition
  this->kmc_data = _kmc_data;
  kinetic_monte_carlo_v2<DebugMode>(
      state, occ_location, *_kmc_data, selected_event, set_selected_event_f,
      apply_selected_event_f, collector, run_manager, event_system);
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
  if (this->options.state_graph_options.has_value() &&
      this->event_group_manager == nullptr) {
    // TODO: construct EventGroupManager
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

/// \brief Update for given state, conditions, occupants, event filters
///
/// - If there are no event groups, just select from the event selector
/// - If there are event groups, determine which event group moves next
template <typename EventSelectorType, bool DebugMode>
void AllowedKineticEventData<EventSelectorType, DebugMode>::select_event(
    SelectedEvent &selected_event, bool requires_event_state) {
  Index selected_event_index;

  if (this->event_group_manager == nullptr) {
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
    auto &m = *this->event_group_manager;
    m.set_next_event_group();

    // Set selected event
    selected_event.group = m.next_event_group();
    log.indent() << "- next_event_group=" << m.next_event_group() << std::endl;
    auto const &g = m.group(m.next_event_group());
    log.indent() << "- g.next_global_event_index=" << g.next_global_event_index
                 << std::endl;
    log.indent() << "- g.next_time=" << g.next_time << std::endl;
    log.indent() << "- g.next_time_increment=" << g.next_time_increment
                 << std::endl;
    log.indent() << "- g.next_state=" << g.next_state << std::endl;
    selected_event_index = g.next_global_event_index;
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

    log.indent() << "- selected_group=" << m.next_event_group() << std::endl;
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

  //{
  //  Log &log = CASM::log();
  //  log.indent() << "- selected={" << prim_event_data.event_type_name << ","
  //               << prim_event_data.equivalent_index << ","
  //               << prim_event_data.is_forward << "," <<
  //               event_id.unitcell_index
  //               << "} / {" << event_id.unitcell_index << ","
  //               << event_id.prim_event_index << "}" << std::endl
  //               << std::endl;
  //  log.end_section();
  //}

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
