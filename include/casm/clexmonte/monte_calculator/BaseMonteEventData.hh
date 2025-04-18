#ifndef CASM_clexmonte_BaseMonteEventData
#define CASM_clexmonte_BaseMonteEventData

#include <random>

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/events/event_data.hh"
#include "casm/clexmonte/events/state_graph.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/monte/methods/kinetic_monte_carlo.hh"
#include "casm/monte/misc/LexicographicalCompare.hh"
#include "casm/monte/sampling/SelectedEventFunctions.hh"

namespace CASM {
namespace clexmonte {

struct EventFilterGroup;
struct StateData;

class EventStateCalculator;

typedef std::function<void(std::reference_wrapper<EventState> state,
                           EventStateCalculator const &calculator)>
    CustomEventStateCalculationFunction;

enum class CALCMETHOD { CONST_ENTROPY, CLEX_ENTROPY, INVALID };

/// \brief Event rate calculation for a particular KMC event
///
/// EventStateCalculator is used to separate the event calculation from the
/// event definition data in PrimEventData. All symmetrically equivalent
/// events can use the same EventStateCalculator, but a simple approach
/// is to create one for each distinct event associated with the primitive
/// cell.
///
/// Two calculation methods are supported:
/// - `CALCMETHOD::CONST_ENTROPY`: uses a constant vibrational entropy
///   assumption and a local cluster expansion (usually just a constant term)
///   sets the attempt frequency.
///   - This is the default, used if no "entropy" cluster expansion is found.
///   - This requires "kra" and "freq" local cluster expansions for each event.
/// - `CALCMETHOD::CLEX_ENTROPY`: uses a cluster expansion for the vibrational
///   entropy and a local cluster expansion for the kinetically-resolved
///   activation entropy, S_kra.
///   - This is used if an "entropy" cluster expansion is found.
///   - This requires "Ef_kra" and "S_kra" local cluster expansions for each
///   event.
class EventStateCalculator {
 public:
  /// \brief Constructor
  EventStateCalculator(std::shared_ptr<system_type> _system,
                       std::string _event_type_name);

  /// \brief Reset pointer to state currently being calculated
  void set(state_type const *state);

  /// \brief Set custom event state calculation function
  void set_custom_event_state_calculation(
      CustomEventStateCalculationFunction f);

  /// \brief Clear custom event state calculation function
  void clear_custom_event_state_calculation();

  /// \brief Pointer to current state
  state_type const *state() const { return m_state; }

  /// \brief Calculate the state of an event
  void calculate_event_state(EventState &state, Index unitcell_index,
                             std::vector<Index> const &linear_site_index,
                             PrimEventData const &prim_event_data) const;

  /// \brief Return the event type name
  std::string const &event_type_name() const { return m_event_type_name; }

  /// \brief Current state's temperature, :math:`T`
  double temperature() const { return *m_temperature; }

  /// \brief The value of :math:`k_{B}T` for the current state
  double kT() const { return CASM::KB * *this->m_temperature; }

  /// \brief Current state's reciprocal temperature, :math:`\beta = 1/(k_{B}T)`
  double beta() const { return 1.0 / (CASM::KB * *this->m_temperature); }

  /// \brief Get the unitcell index for the event currently being calculated
  /// (valid for custom event state calculations only)
  Index curr_unitcell_index() const { return m_unitcell_index; }

  /// \brief Get the linear site indices for the event currently being
  /// calculated (valid for custom event state calculations only)
  std::vector<Index> const &curr_linear_site_index() const {
    return *m_linear_site_index;
  }

  /// \brief Get the prim event data for the event currently being calculated
  /// (valid for custom event state calculations only)
  PrimEventData const &curr_prim_event_data() const {
    return *m_prim_event_data;
  }

  /// Set default event state (short cut for custom event state calculations)
  void set_default_event_state(EventState &state) const {
    _default_event_state_calculation(state, this->curr_unitcell_index(),
                                     this->curr_linear_site_index(),
                                     this->curr_prim_event_data());
  }

  /// Get the formation energy cluster expansion
  std::shared_ptr<clexulator::ClusterExpansion> formation_energy_clex() const {
    return m_formation_energy_clex;
  }

  /// Get the formation energy coefficients
  clexulator::SparseCoefficients const &formation_energy_coefficients() const {
    if (m_formation_energy_clex == nullptr) {
      throw std::runtime_error(
          "EventStateCalculator::formation_energy_coefficients: "
          "m_formation_energy_clex == nullptr");
    }
    return m_formation_energy_clex->coefficients();
  }

  /// Get the entropy cluster expansion
  std::shared_ptr<clexulator::ClusterExpansion> entropy_clex() const {
    return m_entropy_clex;
  }

  /// Get the entropy coefficients
  clexulator::SparseCoefficients const &entropy_coefficients() const {
    if (m_entropy_clex == nullptr) {
      throw std::runtime_error(
          "EventStateCalculator::entropy_coefficients: "
          "m_entropy_clex == nullptr");
    }
    return m_entropy_clex->coefficients();
  }

  /// Get the event multi-local cluster expansion
  std::shared_ptr<clexulator::MultiLocalClusterExpansion> event_clex() const {
    return m_event_clex;
  }

  /// The method used to calculate the event state
  CALCMETHOD calc_method() const { return m_calc_method; }

  /// The index of the event multi-local cluster expansion output that
  /// corresponds to the KRA value
  Index kra_index() const { return m_kra_index; }

  /// The index of the event multi-local cluster expansion output that
  /// corresponds to the attempt frequency value
  Index freq_index() const { return m_freq_index; }

  /// The index of the event multi-local cluster expansion output that
  /// corresponds to the Ef_kra value
  Index Ef_kra_index() const { return m_Ef_kra_index; }

  /// The index of the event multi-local cluster expansion output that
  /// corresponds to the S_kra value
  Index S_kra_index() const { return m_S_kra_index; }

  /// Get the attempt frequency coefficients for a specific event
  clexulator::SparseCoefficients const &freq_coefficients() const {
    if (m_event_clex == nullptr) {
      throw std::runtime_error(
          "EventStateCalculator::freq_coefficients: m_event_clex == nullptr");
    }
    return m_event_clex->coefficients()[m_freq_index];
  }

  /// Get the KRA coefficients for a specific event
  clexulator::SparseCoefficients const &kra_coefficients() const {
    if (m_event_clex == nullptr) {
      throw std::runtime_error(
          "EventStateCalculator::kra_coefficients: m_event_clex == nullptr");
    }
    return m_event_clex->coefficients()[m_kra_index];
  }

  /// Get the Ef_kra coefficients for a specific event
  clexulator::SparseCoefficients const &Ef_kra_coefficients() const {
    if (m_event_clex == nullptr) {
      throw std::runtime_error(
          "EventStateCalculator::Ef_kra_coefficients: m_event_clex == nullptr");
    }
    return m_event_clex->coefficients()[m_Ef_kra_index];
  }

  /// Get the S_kra coefficients for a specific event
  clexulator::SparseCoefficients const &S_kra_coefficients() const {
    if (m_event_clex == nullptr) {
      throw std::runtime_error(
          "EventStateCalculator::S_kra_coefficients: m_event_clex == nullptr");
    }
    return m_event_clex->coefficients()[m_S_kra_index];
  }

 private:
  /// \brief Calculate the state of an event
  void _default_event_state_calculation(
      EventState &state, Index unitcell_index,
      std::vector<Index> const &linear_site_index,
      PrimEventData const &prim_event_data) const;

  /// System pointer
  std::shared_ptr<system_type> m_system;

  /// Event type name
  std::string m_event_type_name;

  /// State to use
  state_type const *m_state;

  /// Current state's temperature
  double const *m_temperature;

  std::shared_ptr<clexulator::ClusterExpansion> m_formation_energy_clex;
  std::shared_ptr<clexulator::ClusterExpansion> m_entropy_clex;
  std::shared_ptr<clexulator::MultiLocalClusterExpansion> m_event_clex;
  mutable Eigen::VectorXd m_event_values;

  CALCMETHOD m_calc_method;

  // -- CALCMETHOD::CONST_ENTROPY: --
  Index m_kra_index;
  Index m_freq_index;

  // -- CALCMETHOD::CLEX_ENTROPY: --
  Index m_Ef_kra_index;
  Index m_S_kra_index;

  /// If true, use custom event state calculation function
  bool m_custom_event_state_calculation;

  /// Custom event state calculation function
  CustomEventStateCalculationFunction m_custom_event_state_calculation_f;

  /// Current event's unitcell index
  mutable Index m_unitcell_index;

  /// Current event's linear site indices
  mutable std::vector<Index> const *m_linear_site_index;

  /// Current event's prim event data
  mutable PrimEventData const *m_prim_event_data;
};

/// \brief A function to control what happens when an abnormal event is
///     encountered or selected
typedef std::function<void(
    Index n_abnormal_events, /*of the current type*/
    std::reference_wrapper<EventState> event_state,
    std::reference_wrapper<EventData const> event_data,
    std::reference_wrapper<PrimEventData const> prim_event_data,
    std::reference_wrapper<state_type const> state)>
    AbnormalEventHandlingFunction;

/// \brief A default AbnormalEventHandlingFunction
///
/// Options:
/// - Print a warning message the first time an abnormal event of a
///   particular type is encountered or selected
/// - Throw an exception when an abnormal event is encountered or
///   selected
/// - Disallow the event when an abnormal event is encountered
/// - Write the local configuration and event state to file when an
///   abnormal event is encountered or selected. The file is named
///   "selected_abnormal_events.jsonl" or "encountered_abnormal_events.jsonl".
///   The format is a "JSON Lines" format, where each line is a JSON object.
///   This is slightly less convenient for reading, but can be written to by
///   just appending.
struct BasicAbnormalEventHandler {
  BasicAbnormalEventHandler(std::string _event_kind, bool _do_throw,
                            bool _do_warn, bool _disallow, Index _n_write,
                            std::optional<fs::path> _output_dir, double _tol);

  /// \brief Check if abnormal events are handled or not based on the
  /// constructor parameters - if not, might be able to avoid re-calculating
  /// event rates when events are selected
  bool handling_on() {
    return (m_do_warn || m_do_throw || m_n_write > 0 || m_disallow);
  }

  /// \brief Handle an abnormal event
  void operator()(Index n_abnormal_events,
                  std::reference_wrapper<EventState> event_state,
                  std::reference_wrapper<EventData const> event_data,
                  std::reference_wrapper<PrimEventData const> prim_event_data,
                  std::reference_wrapper<state_type const> state);

 private:
  /// \brief One of "encountered" or "selected"
  std::string m_event_kind;

  /// \brief If true, throw an exception when an abnormal event is
  ///     encountered or selected
  bool m_do_throw;

  /// \brief If true, warn when an abnormal event is encountered or
  ///     selected
  bool m_do_warn;

  /// \brief If true, disallow the event when an abnormal event is
  ///     encountered (not valid for `m_event_kind`=="selected")
  bool m_disallow;

  /// \brief If >0, write the event to file when an abnormal event is
  ///     encountered or selected, up to `m_n_write` times, only including the
  ///     local configurations for which the local_corr are unique
  int m_n_write;

  /// \brief Output directory
  fs::path m_output_dir;

  /// \brief Log for warning messages (CASM::err_log())
  monte::MethodLog m_event_log;

  /// \brief Tolerance for local_corr comparison (CASM::TOL)
  double m_tol;

  /// \brief Local correlations, by event type name, of local configurations
  ///     that have already been written to file
  std::map<std::string,
           std::set<Eigen::VectorXd, monte::FloatLexicographicalCompare>>
      m_local_corr;

  /// \brief File to write local configurations to
  fs::path m_local_configurations_file;

  /// \brief Read local_corr from an existing file
  void _read_local_corr();
};

struct EventDataOptions {
  /// \brief Output directory (default="output")
  fs::path output_dir;

  /// \brief Tolerance for local_corr comparison (default=CASM::TOL)
  double local_corr_compare_tol;

  /// \brief If true (default), print a warning message when abnormal events
  ///     are encountered.
  bool warn_if_encountered_event_is_abnormal = true;

  /// \brief If true (default), throw an exception when abnormal events
  ///     are encountered.
  bool throw_if_encountered_event_is_abnormal = true;

  /// \brief If true (default), set EventState.rate=0.0 for abnormal events.
  bool disallow_if_encountered_event_is_abnormal = false;

  /// \brief If true (default), write local configurations to file when abnormal
  ///     events are encountered, up to this number.
  Index n_write_if_encountered_event_is_abnormal = 100;

  /// \brief If true (default), print a warning message when abnormal events
  ///     are selected.
  bool warn_if_selected_event_is_abnormal = true;

  /// \brief If true (default), throw an exception when abnormal events
  ///     are selected.
  bool throw_if_selected_event_is_abnormal = true;

  /// \brief If true (default), write local configurations to file when
  ///     abnormal events are selected, up to this number.
  Index n_write_if_selected_event_is_abnormal = 100;

  // --- AllowedEventData options ---

  /// \brief If true, use the map index for the AllowedEventMap; If
  /// false (default), use the vector index
  ///
  /// The map index is lower memory, but may be slower; The vector index is
  /// higher memory, but may be faster.
  ///
  /// Note: current testing does not show this to be helpful - fix to false
  bool use_map_index = false;

  /// \brief If true (default), use the neighborlist impact table; else use the
  /// relative impact table
  ///
  /// Type of impact table:
  /// - Only takes effect for AllowedKineticEventData (event_data_type=default)
  /// - If true: somewhat higher memory use; somewhat faster impact list
  /// - If false: somewhat lower memory use; somewhat slower impact list
  bool use_neighborlist_impact_table = true;

  /// \brief If true (default) check if potentially impacted events are allowed
  ///     and only assign them to the event list if they are. Otherwise,
  ///     assign all potentially impacted events to the event list (whether they
  ///     are allowed will still be checked during the rate calculation).
  ///
  bool assign_allowed_events_only = true;

  // -- State saving options --

  /// \brief State graph / FPTA (First Passage Time Analysis) options
  std::optional<state_graph::Options> state_graph_options = std::nullopt;
};

/// \brief Base class to provide access to event data for a Monte Carlo
/// simulation
class BaseMonteEventData {
 public:
  typedef default_engine_type engine_type;
  typedef monte::KMCData<config_type, statistics_type, engine_type>
      kmc_data_type;
  typedef clexmonte::run_manager_type<engine_type> run_manager_type;

  BaseMonteEventData() = default;
  virtual ~BaseMonteEventData() = default;

  /// The system
  std::shared_ptr<system_type> system;

  /// The `prim events`, one translationally distinct instance
  /// of each event, associated with origin primitive cell
  std::vector<clexmonte::PrimEventData> prim_event_list;

  /// Information about what sites may impact each prim event
  std::vector<clexmonte::EventImpactInfo> prim_impact_info_list;

  /// Custom event state calculation functions
  std::map<std::string, CustomEventStateCalculationFunction>
      custom_event_state_calculation_f;

  /// Function to handle encountered abnormal events
  AbnormalEventHandlingFunction encountered_abnormal_event_handling_f;

  /// \brief If true, the selected abnormal event handling function
  /// is activated
  bool encountered_abnormal_event_handling_on = false;

  /// Function to handle selected abnormal events
  AbnormalEventHandlingFunction selected_abnormal_event_handling_f;

  /// \brief If true, the selected abnormal event handling function
  /// is activated
  bool selected_abnormal_event_handling_on = false;

  /// \brief Count not-normal events (key == event_type_name; value == count)
  ///
  /// An "encounter" with a "not normal" event occurs when (i) the event state
  /// calculation is performed and (ii) the configuration allows the event to
  /// be possible and (iii) there is no barrier between the initial and final
  /// states.
  ///
  /// Note that the KMC implementation may require re-calculating the
  /// rates when the local environment impacts whether the event is possible
  /// the event rate changes, or for other implementation-specific reasons
  /// such as the event list size changing, or for performing an exploration of
  /// the local energy landscape.
  std::map<std::string, Index> n_encountered_abnormal;

  /// \brief Count not-normal events (key == event_type_name; value == count)
  ///
  /// Counts how many "selected" events (used to update the configuration)
  /// are "not normal" (meaning there is no barrier between the initial and
  /// final states).
  std::map<std::string, Index> n_selected_abnormal;

  // -- Options --

  /// \brief Various options, that control event handling (if applicable)
  virtual EventDataOptions const &event_data_options() const = 0;

  // -- Customize event state calculation and handling functions --

  /// \brief Set a custom event state calculation function
  void set_custom_event_state_calculation(
      std::string const &event_type_name,
      CustomEventStateCalculationFunction f) {
    custom_event_state_calculation_f[event_type_name] = f;
  }

  /// \brief Erase a custom event state calculation function
  void set_custom_event_state_calculation_off(
      std::string const &event_type_name) {
    custom_event_state_calculation_f.erase(event_type_name);
  }

  /// \brief Set the encountered abnormal event handling function
  void set_encountered_abnormal_event_handling(
      AbnormalEventHandlingFunction handling_f) {
    encountered_abnormal_event_handling_f = handling_f;
    encountered_abnormal_event_handling_on = true;
  }

  /// \brief Turn off encountered abnormal event handling
  void set_encountered_abnormal_event_handling_off() {
    encountered_abnormal_event_handling_f = nullptr;
    encountered_abnormal_event_handling_on = false;
  }

  /// \brief Set the selected abnormal event handling function
  void set_selected_abnormal_event_handling(
      AbnormalEventHandlingFunction handling_f) {
    selected_abnormal_event_handling_f = handling_f;
    selected_abnormal_event_handling_on = true;
  }

  /// \brief Turn off selected abnormal event handling
  void set_selected_abnormal_event_handling_off() {
    selected_abnormal_event_handling_f = nullptr;
    selected_abnormal_event_handling_on = false;
  }

  // -- Update and run --

  virtual void update(
      std::shared_ptr<StateData> _state_data,
      std::optional<std::vector<EventFilterGroup>> _event_filters,
      std::shared_ptr<engine_type> engine) = 0;

  virtual void run(state_type &state, monte::OccLocation &occ_location,
                   SelectedEvent &selected_event,
                   std::optional<monte::SelectedEventDataCollector> &collector,
                   run_manager_type &run_manager,
                   std::shared_ptr<kmc_data_type> _kmc_data,
                   std::shared_ptr<occ_events::OccSystem> event_system) = 0;

  // -- Select Event --

  /// Select an event to apply
  virtual void select_event(SelectedEvent &selected_event,
                            bool requires_event_state) = 0;

  // -- Event list summary info --

  /// The size of the event list
  virtual Index n_events() const = 0;

  /// Return the current total event rate
  virtual double total_rate() const = 0;

  // -- Event list iteration --

  /// Construct new internal iterator and return its index
  virtual Index new_iterator(bool is_end) = 0;

  /// Copy internal iterator and return the new iterator index
  virtual Index copy_iterator(Index i) = 0;

  /// Erase internal iterator
  virtual void erase_iterator(Index i) = 0;

  /// Check if two internal iterators are equal
  virtual bool equal_iterator(Index i, Index j) = 0;

  /// Advance internal iterator by one event
  virtual void advance_iterator(Index i) = 0;

  /// The event ID for the current state of the internal iterator
  virtual EventID const &event_id(Index i) const = 0;

  // -- Event info (accessed by EventID) --

  /// The monte::OccEvent that can apply the specified event. Reference is
  /// valid until the next call to this method.
  virtual monte::OccEvent const &event_to_apply(EventID const &id) const = 0;

  /// Return the current rate for a specific event
  virtual double event_rate(EventID const &id) const = 0;

  /// Calculate event state data. Reference is valid until the next call to this
  /// method.
  virtual EventState const &event_state(EventID const &id) const = 0;

  /// The events that must be updated if the specified event occurs. Reference
  /// is valid until the next call to this method.
  virtual std::vector<EventID> const &impact(EventID const &id) const = 0;
};

}  // namespace clexmonte
}  // namespace CASM

#endif
