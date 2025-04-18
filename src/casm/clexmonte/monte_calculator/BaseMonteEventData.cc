#include "casm/clexmonte/monte_calculator/BaseMonteEventData.hh"

// #include "casm/casm_io/SafeOfstream.hh"
#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/events/event_methods.hh"
#include "casm/clexmonte/events/io/json/event_data_json_io.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/system/System.hh"
// #include "casm/clexmonte/system/io/json/system_data_json_io.hh"
#include "casm/configuration/Configuration.hh"
#include "casm/configuration/io/json/Configuration_json_io.hh"
// #include "casm/configuration/occ_events/io/json/OccEvent_json_io.hh"
#include "casm/monte/run_management/State.hh"

namespace CASM {
namespace clexmonte {

/// \brief Constructor
EventStateCalculator::EventStateCalculator(std::shared_ptr<system_type> _system,
                                           std::string _event_type_name)
    : m_system(_system),
      m_event_type_name(_event_type_name),
      m_custom_event_state_calculation(false),
      m_custom_event_state_calculation_f(nullptr),
      m_calc_method(CALCMETHOD::INVALID) {}

/// \brief Reset pointer to state currently being calculated
void EventStateCalculator::set(state_type const *state) {
  // supercell-specific
  m_state = state;
  if (m_state == nullptr) {
    throw std::runtime_error(
        "Error setting EventStateCalculator state: state is empty");
  }
  m_temperature = &m_state->conditions.scalar_values.at("temperature");
  m_formation_energy_clex = get_clex(*m_system, *m_state, "formation_energy");

  // Option: CALCMETHOD::CLEX_ENTROPY if a "entropy" cluster expansion is found
  auto const &supercell_clex = get_supercell_data(*m_system, *m_state).clex;
  auto it = supercell_clex.find("entropy");
  if (it == supercell_clex.end()) {
    m_calc_method = CALCMETHOD::CONST_ENTROPY;
  } else {
    m_calc_method = CALCMETHOD::CLEX_ENTROPY;
    m_entropy_clex = it->second;
    CASM::clexmonte::set(*m_entropy_clex, *m_state);
  }

  // set and validate event clex
  LocalMultiClexData event_local_multiclex_data =
      get_local_multiclex_data(*m_system, m_event_type_name);
  m_event_clex = get_local_multiclex(*m_system, *m_state, m_event_type_name);
  m_event_values.resize(m_event_clex->coefficients().size());
  std::map<std::string, Index> _glossary =
      event_local_multiclex_data.coefficients_glossary;

  auto _check_coeffs = [&](Index &coeff_index, std::string key) {
    if (!_glossary.count(key)) {
      std::stringstream ss;
      ss << "Error constructing " << m_event_type_name
         << " EventStateCalculator: No " << key << " cluster expansion";
      throw std::runtime_error(ss.str());
    }
    coeff_index = _glossary.at(key);
    if (coeff_index < 0 || coeff_index >= m_event_clex->coefficients().size()) {
      std::stringstream ss;
      ss << "Error constructing " << m_event_type_name
         << " EventStateCalculator: " << key << " index out of range";
      throw std::runtime_error(ss.str());
    }
  };

  if (m_calc_method == CALCMETHOD::CONST_ENTROPY) {
    _check_coeffs(m_kra_index, "kra");
    _check_coeffs(m_freq_index, "freq");
  } else if (m_calc_method == CALCMETHOD::CLEX_ENTROPY) {
    _check_coeffs(m_Ef_kra_index, "Ef_kra");
    _check_coeffs(m_S_kra_index, "S_kra");
  }
}

/// \brief Set custom event state calculation function
void EventStateCalculator::set_custom_event_state_calculation(
    CustomEventStateCalculationFunction f) {
  m_custom_event_state_calculation = true;
  m_custom_event_state_calculation_f = f;
}

/// \brief Clear custom event state calculation function
void EventStateCalculator::clear_custom_event_state_calculation() {
  m_custom_event_state_calculation = false;
  m_custom_event_state_calculation_f = nullptr;
}

/// \brief Calculate the state of an event
///
/// If a custom event state calculation function is set, it is called; otherwise
/// the default event state calculation is used.
///
/// \param state
/// \param unitcell_index
/// \param linear_site_index
/// \param prim_event_data
///
void EventStateCalculator::calculate_event_state(
    EventState &state, Index unitcell_index,
    std::vector<Index> const &linear_site_index,
    PrimEventData const &prim_event_data) const {
  // Initialize event state
  state.formation_energy_delta_corr = nullptr;
  state.local_corr = nullptr;

  // Check if event is allowed based on current occupation
  clexulator::ConfigDoFValues const *dof_values =
      m_formation_energy_clex->get();
  state.is_allowed =
      event_is_allowed(linear_site_index, *dof_values, prim_event_data);
  if (!state.is_allowed) {
    state.rate = 0.0;
    return;
  }

  // Calculate event state
  if (this->m_custom_event_state_calculation) {
    // Set current event details for access by custom event state calculation:
    m_unitcell_index = unitcell_index;
    m_linear_site_index = &linear_site_index;
    m_prim_event_data = &prim_event_data;

    // Call custom event state calculation function
    this->m_custom_event_state_calculation_f(std::ref(state), *this);
    return;
  } else {
    this->_default_event_state_calculation(state, unitcell_index,
                                           linear_site_index, prim_event_data);
  }
}

/// \brief Calculate the state of an event
void EventStateCalculator::_default_event_state_calculation(
    EventState &state, Index unitcell_index,
    std::vector<Index> const &linear_site_index,
    PrimEventData const &prim_event_data) const {
  // calculate change in energy to final state
  //  state.dE_final = m_formation_energy_clex->occ_delta_value(
  //      event_data.event.linear_site_index, prim_event_data.occ_final);

  if (m_calc_method == CALCMETHOD::CONST_ENTROPY) {
    // calculate change in energy to final state
    // - and save pointer to delta correlations
    state.formation_energy_delta_corr =
        &m_formation_energy_clex->correlations().occ_delta(
            linear_site_index, prim_event_data.occ_final);
    state.dEf_final = m_formation_energy_clex->coefficients() *
                      (*state.formation_energy_delta_corr);

    // calculate KRA and attempt frequency
    // - add save pointer to local correlations
    state.local_corr = &m_event_clex->correlations().local(
        unitcell_index, prim_event_data.equivalent_index);
    for (int i = 0; i < m_event_clex->coefficients().size(); ++i) {
      m_event_values(i) = m_event_clex->coefficients()[i] * (*state.local_corr);
    }
    state.Ef_kra = m_event_values[m_kra_index];
    state.freq = m_event_values[m_freq_index];
    state.reverse_freq = state.freq;

    // calculate energy in activated state, check if "normal", calculate rate
    state.dEf_activated = state.dEf_final * 0.5 + state.Ef_kra;

    //
    state.dE_activated = state.dEf_activated;
    state.dE_final = state.dEf_final;

    state.is_normal =
        (state.dE_activated > 0.0) && (state.dE_activated > state.dE_final);
    if (state.dE_activated < state.dE_final)
      state.dE_activated = state.dE_final;
    if (state.dE_activated < 0.0) state.dE_activated = 0.0;

    // rate = freq * exp(-beta * dEf_activated)
    //      = freq * exp(-beta * dE_activated)
    // reverse_rate = reverse_freq * exp(-beta * (dEf_activated - dEf_final))
    //              = reverse_rate * exp(beta * (dE_activated - dE_final))

    // calculate rate
    state.rate = state.freq * exp(-this->beta() * state.dE_activated);

    // calculate reverse rate
    state.reverse_rate =
        state.reverse_freq *
        exp(-this->beta() * (state.dE_activated - state.dE_final));

  } else if (m_calc_method == CALCMETHOD::CLEX_ENTROPY) {
    // calculate change in energy and entropy to final state
    // - and save pointer to delta correlations
    state.formation_energy_delta_corr =
        &m_formation_energy_clex->correlations().occ_delta(
            linear_site_index, prim_event_data.occ_final);
    state.dEf_final = m_formation_energy_clex->coefficients() *
                      (*state.formation_energy_delta_corr);
    state.entropy_delta_corr = &m_entropy_clex->correlations().occ_delta(
        linear_site_index, prim_event_data.occ_final);
    state.dS_final =
        m_entropy_clex->coefficients() * (*state.entropy_delta_corr);

    // calculate Ef_kra and S_kra
    // - add save pointer to local correlations
    state.local_corr = &m_event_clex->correlations().local(
        unitcell_index, prim_event_data.equivalent_index);
    for (int i = 0; i < m_event_clex->coefficients().size(); ++i) {
      m_event_values(i) = m_event_clex->coefficients()[i] * (*state.local_corr);
    }
    state.Ef_kra = m_event_values[m_Ef_kra_index];
    state.S_kra = m_event_values[m_S_kra_index];

    // calculate energy in activated state, check if "normal", calculate rate
    state.dEf_activated = state.dEf_final * 0.5 + state.Ef_kra;
    state.dS_activated = state.dS_final * 0.5 + state.S_kra;

    //
    state.dE_activated =
        state.dEf_activated - this->temperature() * state.dS_activated;
    state.dE_final = state.dEf_final - this->temperature() * state.dS_final;

    state.is_normal =
        (state.dE_activated > 0.0) && (state.dE_activated > state.dE_final);
    if (state.dE_activated < state.dE_final) {
      state.dE_activated = state.dE_final;
    }
    if (state.dE_activated < 0.0) {
      state.dE_activated = 0.0;
    }

    // rate = (kT/h) * exp(beta * T * dS_activated) * exp(-beta * dEf_activated)
    //      = (kT/h) * exp(-beta * (dEf_activated - T * dS_activated))
    //      = (kT/h) * exp(-beta * dE_activated)
    // reverse_rate = (kT/h) * exp(beta * T * (dS_activated-dS_final))
    //                * exp(-beta * (dEf_activated - dEf_final))
    // reverse_rate = (kT/h) * exp(
    //   -beta * (
    //     (dEf_activated - T * dS_activated) - (dEf_final - T * dS_final)
    //   )
    // )
    //              = (kT/h) * exp(-beta * (dE_activated - dE_final))

    // calculate attempt frequency
    state.freq =
        (this->kT() / CASM::PLANCK) * exp(state.dS_activated / CASM::KB);
    state.reverse_freq = (this->kT() / CASM::PLANCK) *
                         exp((state.dS_activated - state.dS_final) / CASM::KB);

    // calculate rate
    double prefactor = this->kT() / CASM::PLANCK;
    state.rate = prefactor * exp(-this->beta() * state.dE_activated);

    state.reverse_rate =
        prefactor * exp(-this->beta() * (state.dE_activated - state.dE_final));

  } else {
    throw std::runtime_error(
        "Error in EventStateCalculator::_default_event_state_calculation: "
        "invalid CALCMETHOD");
  }
}

namespace {

/// \brief Create JSON object for local configuration and event data
///
/// Format:
/// \code
/// "local_configuration": libcasm.local_configuration.LocalConfiguration
/// \endcode
///
jsonParser local_configuration_to_json(state_type const &state,
                                       EventState const &event_state,
                                       EventData const &event_data,
                                       PrimEventData const &prim_event_data) {
  jsonParser local_config_json;
  bool write_prim_basis = false;
  to_json(state.configuration, local_config_json["configuration"],
          write_prim_basis);
  local_config_json["pos"] = std::vector<Index>(
      {event_data.unitcell_index, prim_event_data.equivalent_index});
  return local_config_json;
}

/// \brief Create JSON object for local configuration and event data
///
/// Format:
///
/// \code
/// {
///   "event_state": libcasm.clexmonte.EventState
///   "unitcell_index": int
///   "linear_site_index": list[int]
///   "prim_event_data": libcasm.clexmonte.PrimEventData
/// }
/// \endcode
jsonParser event_to_json(state_type const &state, EventState const &event_state,
                         EventData const &event_data,
                         PrimEventData const &prim_event_data) {
  jsonParser event_json;
  to_json(event_state, event_json["event_state"]);
  event_json["unitcell_index"] = event_data.unitcell_index;
  event_json["linear_site_index"] = event_data.event.linear_site_index;
  to_json(prim_event_data, event_json["prim_event_data"]);
  return event_json;
}

}  // namespace

/// \brief Constructor
///
/// \param _event_kind One of "encountered" or "selected"
/// \param _do_throw If true, throw an exception when an abnormal event
///     is encountered or selected
/// \param _do_warn If true, warn when an abnormal event is encountered
///     or selected
/// \param _disallow If true, disallow the event when an abnormal event
///     is encountered (not valid for `m_event_kind`=="selected")
/// \param _n_write If >0, write the event to file when an event without
///     barrier is encountered or selected, up to `m_n_write` times, only
///     including the local configurations for which the local_corr are unique
/// \param _output_dir Output directory
/// \param _tol Tolerance for local_corr comparison
///
BasicAbnormalEventHandler::BasicAbnormalEventHandler(
    std::string _event_kind, bool _do_throw, bool _do_warn, bool _disallow,
    Index _n_write, std::optional<fs::path> _output_dir, double _tol)
    : m_event_kind(_event_kind),
      m_do_throw(_do_throw),
      m_do_warn(_do_warn),
      m_disallow(_disallow),
      m_n_write(_n_write),
      m_output_dir(_output_dir.value_or("output")),
      m_tol(_tol) {
  if (m_event_kind != "encountered" && m_event_kind != "selected") {
    throw std::runtime_error(
        "Error in BasicAbnormalEventHandler: event_kind must be "
        "'encountered' or 'selected'");
  }
  if (m_disallow && m_event_kind == "selected") {
    throw std::runtime_error(
        "Error in BasicAbnormalEventHandler: for `event_kind`=="
        "\"selected\", `disallow` cannot be true");
  }

  m_local_configurations_file =
      m_output_dir / (m_event_kind + "_abnormal_events.jsonl");

  if (m_n_write > 0) {
    _read_local_corr();
  }

  m_event_log.log = CASM::err_log();
}

/// \brief Handle an abnormal event
///
/// \param n_abnormal_events Number of events without barrier of the
///     current event type
/// \param event_state Event state
/// \param event_data Event data
/// \param prim_event_data Prim event data
/// \param state Current state
///
void BasicAbnormalEventHandler::operator()(
    Index n_abnormal_events, std::reference_wrapper<EventState> event_state,
    std::reference_wrapper<EventData const> event_data,
    std::reference_wrapper<PrimEventData const> prim_event_data,
    std::reference_wrapper<state_type const> state) {
  // Warning message ----------------------------------------------------
  if (m_do_warn && n_abnormal_events == 1) {
    Log &log = m_event_log.log;
    log << "\n";
    if (m_event_kind == "selected") {
      log << "## WARNING: SELECTED ABNORMAL EVENT #################\n"
             "#                                                   #\n"
             "# - The event was selected.                         #\n"
             "#                                                   #\n";
    } else {
      log << "## WARNING: ENCOUNTERED ABNORMAL EVENT ##############\n"
             "#                                                   #\n"
             "# - The event was encountered when calculating      #\n"
             "#   event rates.                                    #\n"
             "# - The event might not be selected.                #\n"
             "#                                                   #\n";
    }

    log << "# This warning is only printed once per event type. #\n"
           "#                                                   #\n"
           "# Event info:                                       #\n";

    jsonParser event_json =
        event_to_json(state.get(), event_state.get(), event_data.get(),
                      prim_event_data.get());

    log << event_json << "\n"
        << "#                                                   #\n"
           "#####################################################\n"
        << std::endl;
  }

  // Write local configuration --------------------------------------------

  // Write up to `n_write` selected events of each type
  // to "selected_abnormal_events.jsonl" file
  std::string const &event_type_name = prim_event_data.get().event_type_name;

  auto local_corr_it = m_local_corr.find(event_type_name);
  if (local_corr_it == m_local_corr.end()) {
    monte::FloatLexicographicalCompare compare(m_tol);
    local_corr_it = m_local_corr.emplace(event_type_name, compare).first;
  }

  if (local_corr_it->second.size() < m_n_write) {
    Eigen::VectorXd const *local_corr = event_state.get().local_corr;
    if (local_corr == nullptr) {
      throw std::runtime_error(
          "Error in BasicAbnormalEventHandler: local_corr==nullptr");
    }

    auto result = local_corr_it->second.insert(*local_corr);
    if (result.second) {
      jsonParser json;
      json["local_configuration"] =
          local_configuration_to_json(state.get(), event_state.get(),
                                      event_data.get(), prim_event_data.get());
      json["event"] = event_to_json(state.get(), event_state.get(),
                                    event_data.get(), prim_event_data.get());

      // Append `json` to file on a new line:
      fs::create_directories(m_output_dir);
      std::ofstream file(m_local_configurations_file, std::ios::app);
      json.print(file, -1 /*indent=-1 -> no new lines*/);
      file << "\n";
    }
  }

  // Throw --------------------------------------------
  if (m_do_throw) {
    throw std::runtime_error("Error: " + m_event_kind +
                             " abnormal event, which is not allowed.");
  }

  // Disallow --------------------------------------------
  if (m_disallow) {
    event_state.get().rate = 0.0;
  }
}

/// \brief Read local_corr from an existing file
void BasicAbnormalEventHandler::_read_local_corr() {
  // Read existing file, if it exists:
  if (fs::exists(m_local_configurations_file)) {
    std::ifstream ifs(m_local_configurations_file);
    // Read each line as a separate JSON object:
    std::string line;
    Index i_line = 0;
    while (std::getline(ifs, line)) {
      jsonParser json = jsonParser::parse(line);
      ParentInputParser parser(json);

      std::string event_type_name;
      parser.require(event_type_name,
                     fs::path("event") / "prim_event_data" / "event_type_name");

      Eigen::VectorXd local_corr;
      parser.require(local_corr,
                     fs::path("event") / "event_state" / "local_corr");

      auto local_corr_it = m_local_corr.find(event_type_name);
      if (local_corr_it == m_local_corr.end()) {
        monte::FloatLexicographicalCompare compare(m_tol);
        local_corr_it = m_local_corr.emplace(event_type_name, compare).first;
      }
      local_corr_it->second.insert(local_corr);

      i_line += 1;
    }
  }
}

}  // namespace clexmonte
}  // namespace CASM
