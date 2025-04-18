#ifndef CASM_clexmonte_MonteEventData
#define CASM_clexmonte_MonteEventData

#include "casm/casm_io/Log.hh"
#include "casm/clexmonte/monte_calculator/BaseMonteEventData.hh"
#include "casm/clexmonte/monte_calculator/StateData.hh"
#include "casm/monte/sampling/SelectedEventFunctions.hh"

namespace CASM {

class Log;

size_t memory_used(bool resident = false);
std::string convert_size(size_t size_bytes);

namespace clexmonte {

class MonteEventListIterator {
 public:
  using iterator_category = std::input_iterator_tag;
  using difference_type = std::ptrdiff_t;
  using value_type = EventID;
  using pointer = EventID const *;    // or also value_type*
  using reference = EventID const &;  // or also value_type&

  MonteEventListIterator() : m_data(nullptr), m_lib(nullptr), m_index(-1) {}

  MonteEventListIterator(std::shared_ptr<BaseMonteEventData> _data,
                         std::shared_ptr<RuntimeLibrary> _lib, bool _is_end)
      : m_data(_data), m_lib(_lib), m_index(m_data->new_iterator(_is_end)) {}

  ~MonteEventListIterator() {
    // ensure BaseMonteEventData is deleted before library
    m_data.reset();
  }

  MonteEventListIterator &operator=(MonteEventListIterator const &other) {
    m_data = other.m_data;
    m_lib = other.m_lib;
    m_index = m_data->copy_iterator(other.m_index);
    return *this;
  }

  reference operator*() const { return m_data->event_id(m_index); }
  pointer operator->() { return &m_data->event_id(m_index); }

  // Prefix increment (++it)
  MonteEventListIterator &operator++() {
    m_data->advance_iterator(m_index);
    return *this;
  }

  // Postfix increment (it++)
  MonteEventListIterator operator++(int) {
    MonteEventListIterator tmp = *this;
    tmp.m_index = m_data->copy_iterator(m_index);
    m_data->advance_iterator(m_index);
    return tmp;
  }

  bool operator==(MonteEventListIterator const &other) const {
    if (m_data == nullptr) {
      return false;
    }
    return m_data == other.m_data &&
           m_data->equal_iterator(m_index, other.m_index);
  }

  bool operator!=(MonteEventListIterator const &other) const {
    return !(*this == other);
  }

 private:
  std::shared_ptr<BaseMonteEventData> m_data;
  std::shared_ptr<RuntimeLibrary> m_lib;
  Index m_index;
};

/// Allows iterating over EventID
class MonteEventList {
 public:
  MonteEventList(std::shared_ptr<BaseMonteEventData> _data,
                 std::shared_ptr<RuntimeLibrary> _lib)
      : m_data(_data), m_lib(_lib) {}

  ~MonteEventList() {
    // ensure BaseMonteEventData is deleted before library
    m_data.reset();
  }

  MonteEventListIterator begin() const {
    return MonteEventListIterator(m_data, m_lib, false);
  }

  MonteEventListIterator end() const {
    return MonteEventListIterator(m_data, m_lib, true);
  }

  /// The number of events
  Index size() const {
    if (!m_data) {
      throw std::runtime_error("Error in MonteEventList::size: No data");
    }
    return m_data->n_events();
  }

  /// The current total event rate
  double total_rate() const { return m_data->total_rate(); }

 private:
  std::shared_ptr<BaseMonteEventData> m_data;
  std::shared_ptr<RuntimeLibrary> m_lib;
};

class MonteEventData {
 public:
  MonteEventData(std::shared_ptr<BaseMonteEventData> _data,
                 std::shared_ptr<RuntimeLibrary> _lib)
      : m_data(_data), m_lib(_lib), m_event_list(_data, _lib) {}

  ~MonteEventData() {
    // ensure BaseMonteEventData is deleted before library
    m_data.reset();
  }

  /// The system
  std::shared_ptr<system_type> system() const { return m_data->system; }

  /// Output directory
  fs::path const &output_dir() const {
    return m_data->event_data_options().output_dir;
  }

  /// The `prim events`, one translationally distinct instance
  /// of each event, associated with origin primitive cell
  std::vector<clexmonte::PrimEventData> const &prim_event_list() const {
    return m_data->prim_event_list;
  }

  /// Information about what sites may impact each prim event
  std::vector<clexmonte::EventImpactInfo> const &prim_impact_info_list() const {
    return m_data->prim_impact_info_list;
  }

  // -- Custom event state calculation and handling functions --

  /// Set custom event state calculation functions
  void set_custom_event_state_calculation(
      std::string const &event_type_name,
      CustomEventStateCalculationFunction f) {
    m_data->set_custom_event_state_calculation(event_type_name, f);
  }

  /// Set custom event state calculation functions off
  void set_custom_event_state_calculation_off(
      std::string const &event_type_name) {
    m_data->set_custom_event_state_calculation_off(event_type_name);
  }

  /// \brief Set the encountered abnormal event handling function
  void set_encountered_abnormal_event_handling(
      AbnormalEventHandlingFunction handling_f) {
    m_data->set_encountered_abnormal_event_handling(handling_f);
  }

  /// \brief Turn off the encountered abnormal event handling function
  void set_encountered_abnormal_event_handling_off() {
    m_data->set_encountered_abnormal_event_handling_off();
  }

  /// \brief Set the selected abnormal event handling function
  void set_selected_abnormal_event_handling(
      AbnormalEventHandlingFunction handling_f) {
    m_data->set_selected_abnormal_event_handling(handling_f);
  }

  /// \brief Turn off the selected abnormal event handling function
  void set_selected_abnormal_event_handling_off() {
    m_data->set_selected_abnormal_event_handling_off();
  }

  /// \brief Turn off both the encountered and selected abnormal event
  /// handling functions
  void set_abnormal_event_handling_off() {
    m_data->set_encountered_abnormal_event_handling_off();
    m_data->set_selected_abnormal_event_handling_off();
  }

  /// Return number of encountered abnormal events calculated, by type
  std::map<std::string, Index> const &n_encountered_abnormal() {
    return m_data->n_encountered_abnormal;
  }

  /// Return number of selected abnormal events, by type
  std::map<std::string, Index> const &n_selected_abnormal() const {
    return m_data->n_selected_abnormal;
  }

  // -- Event list --

  /// Get EventID
  MonteEventList const &event_list() const { return m_event_list; }

  // -- Event info (accessed by EventID) --

  /// The monte::OccEvent that can apply the event for the current state of the
  /// internal iterator
  monte::OccEvent const &event_to_apply(EventID const &id) const {
    return m_data->event_to_apply(id);
  }

  /// Return the current rate for a specific event
  double event_rate(EventID const &id) const { return m_data->event_rate(id); }

  /// Calculate event state data
  EventState const &event_state(EventID const &id) const {
    return m_data->event_state(id);
  }

  /// The events that must be updated if the specified event occurs
  std::vector<EventID> const &event_impact(EventID const &id) const {
    return m_data->impact(id);
  }

 private:
  std::shared_ptr<BaseMonteEventData> m_data;
  std::shared_ptr<RuntimeLibrary> m_lib;
  MonteEventList m_event_list;
};

struct EventTypeStats {
  typedef monte::CountType CountType;

  EventTypeStats(
      std::vector<std::string> const &_partion_names_by_type,
      std::vector<std::string> const &_partion_names_by_equivalent_index,
      double _initial_begin, double _bin_width, bool _is_log,
      Index _max_size = 10000);

  CountType n_total;
  double min;
  double max;
  double sum;
  double mean;
  monte::PartitionedHistogram1D hist_by_type;
  monte::PartitionedHistogram1D hist_by_equivalent_index;

  void insert(int partition_by_type, int partition_by_equivalent_index,
              double value);
};

struct EventDataSummary {
  // -- Primary data --

  std::shared_ptr<StateData> state_data;
  MonteEventData event_data;
  std::vector<PrimEventData> const &prim_event_list;

  typedef monte::CountType CountType;
  typedef std::string TypeKey;
  typedef std::pair<std::string, Index> EquivKey;

  TypeKey type_key(Index prim_event_index) const {
    return prim_event_list.at(prim_event_index).event_type_name;
  }

  TypeKey type_key(EventID const &id) const {
    return type_key(id.prim_event_index);
  }

  EquivKey equiv_key(Index prim_event_index) const {
    return std::make_pair(
        prim_event_list.at(prim_event_index).event_type_name,
        prim_event_list.at(prim_event_index).equivalent_index);
  }

  EquivKey equiv_key(EventID const &id) const {
    return equiv_key(id.prim_event_index);
  }

  EventDataSummary(std::shared_ptr<StateData> const &_state_data,
                   MonteEventData const &_event_data,
                   double energy_bin_width = 0.1, double freq_bin_width = 0.1,
                   double rate_bin_width = 0.1);

  /// The event type index for each prim event index
  std::vector<Index> to_event_type;
  std::vector<std::string> event_type_names;

  /// The {event type + equivalent index} index for each prim event index
  std::vector<Index> to_equivalent_index;
  std::vector<std::string> equivalent_index_names;

  std::set<TypeKey> all_types;
  std::map<TypeKey, std::set<EquivKey>> equiv_keys_by_type;
  std::set<EquivKey> all_equiv_keys;

  struct IntCountByType {
    std::map<TypeKey, Index> by_type;
    std::map<EquivKey, Index> by_equivalent_index;
  };

  struct FloatCountByType {
    std::map<TypeKey, double> by_type;
    std::map<EquivKey, double> by_equivalent_index;
  };

  struct ImpactTable {
    // value: number impacted events of type (first key) when type (second key)
    // occurs
    std::map<TypeKey, std::map<TypeKey, double>> by_type;
    std::map<EquivKey, std::map<EquivKey, double>> by_equivalent_index;
  };

  // -- Count & total rate info --

  Index n_events_allowed;
  Index n_events_possible;
  Index n_abnormal_total;
  Index event_list_size;
  double total_rate;
  double mean_time_increment;

  IntCountByType n_possible;
  IntCountByType n_allowed;
  IntCountByType n_abnormal;

  FloatCountByType rate;

  // -- Memory usage --

  size_t resident_bytes_used;
  double resident_MiB_used;

  // -- Impact info --

  //  FloatCountByType n_impact;

  /// The number of impacted events of type (second key) when type (first key)
  /// occurs. Values are only added for allowed events, so n_allowed can be
  /// used to normalize.
  ImpactTable impact_table;

  std::map<TypeKey, CountType> neighborhood_size_total;
  std::map<TypeKey, CountType> neighborhood_size_formation_energy;
  std::map<TypeKey, CountType> neighborhood_size_kra;
  std::map<TypeKey, CountType> neighborhood_size_freq;

  // -- State info --

  std::vector<std::string> stats_labels;
  std::vector<EventTypeStats> stats;

 private:
  void _add_count(EventID const &id, EventState const &state);
  void _add_impact(EventID const &id, EventState const &state);
  void _add_stats(EventID const &id, EventState const &state);
};

template <int VerbosityLevel = Log::standard>
void print(Log &log, EventDataSummary const &event_data_summary);

}  // namespace clexmonte

jsonParser &to_json(clexmonte::EventTypeStats const &stats, jsonParser &json);

jsonParser &to_json(clexmonte::EventDataSummary::IntCountByType const &count,
                    jsonParser &json);

jsonParser &to_json(clexmonte::EventDataSummary::FloatCountByType const &count,
                    jsonParser &json);

jsonParser &to_json(clexmonte::EventDataSummary const &event_data_summary,
                    jsonParser &json);

}  // namespace CASM

// --- Implementation ---

namespace CASM {
namespace clexmonte {

template <int VerbosityLevel>
void print(Log &log, EventDataSummary const &event_data_summary) {
  EventDataSummary const &x = event_data_summary;
  log.begin_section<VerbosityLevel>();
  log.indent() << "Event data summary:\n";
  log.indent() << "- Number of unitcells = " << x.state_data->n_unitcells
               << std::endl;
  // Number of events
  log.indent() << "- Number of events (total) = " << x.n_events_allowed
               << std::endl;
  log.indent() << "- Number of events (by type):" << std::endl;
  for (auto const &pair : x.equiv_keys_by_type) {
    std::string event_type_name = pair.first;
    log.indent() << "  - " << event_type_name << " = "
                 << x.n_allowed.by_type.at(event_type_name);
    if (x.n_abnormal.by_type.at(event_type_name) != 0.0) {
      log.indent() << " (** abnormal = "
                   << x.n_abnormal.by_type.at(event_type_name) << " **)";
    }
    log.indent() << std::endl;

    for (auto const &equiv_key : pair.second) {
      Index equivalent_index = equiv_key.second;
      log.indent() << "    - " << event_type_name << "." << equivalent_index
                   << " = " << x.n_allowed.by_equivalent_index.at(equiv_key);
      if (x.n_abnormal.by_equivalent_index.at(equiv_key) != 0.0) {
        log.indent() << " (** abnormal = "
                     << x.n_abnormal.by_equivalent_index.at(equiv_key)
                     << " **)";
      }
      log.indent() << std::endl;
    }
  }
  // Rate info:
  log.indent() << "- Event rate (total) (1/s) = " << x.total_rate << std::endl;
  log.indent() << "- Event rate (by type) (1/s):" << std::endl;
  for (auto const &pair : x.equiv_keys_by_type) {
    std::string event_type_name = pair.first;
    log.indent() << "  - " << event_type_name << " = "
                 << x.rate.by_type.at(event_type_name) << std::endl;

    for (auto const &equiv_key : pair.second) {
      Index equivalent_index = equiv_key.second;
      log.indent() << "    - " << event_type_name << "." << equivalent_index
                   << " = " << x.rate.by_equivalent_index.at(equiv_key)
                   << std::endl;
    }
  }
  log.indent() << "- Mean time increment (total) (s) = "
               << x.mean_time_increment << std::endl;
  log.indent() << "- Mean time increment (by type) (s):" << std::endl;
  for (auto const &pair : x.equiv_keys_by_type) {
    std::string event_type_name = pair.first;
    log.indent() << "  - " << event_type_name << " = "
                 << 1.0 / x.rate.by_type.at(event_type_name) << std::endl;

    for (auto const &equiv_key : pair.second) {
      Index equivalent_index = equiv_key.second;
      log.indent() << "    - " << event_type_name << "." << equivalent_index
                   << " = " << 1.0 / x.rate.by_equivalent_index.at(equiv_key)
                   << std::endl;
    }
  }
  // Memory usage:
  log.indent() << "- Memory used (RAM) = "
               << convert_size(x.resident_bytes_used) << std::endl;
  log.indent() << "- Event list size = " << x.event_list_size << std::endl;
  // Impact neighborhood sizes:
  log.indent() << "- Impact neighborhood sizes (#sites): total (Ef / Ef_kra / "
                  "freq)"
               << std::endl;
  for (auto const &pair : x.neighborhood_size_total) {
    std::string event_type_name = pair.first;
    log.indent() << "  - " << event_type_name << " = " << pair.second << " ("
                 << x.neighborhood_size_formation_energy.at(event_type_name)
                 << " / " << x.neighborhood_size_kra.at(event_type_name)
                 << " / " << x.neighborhood_size_freq.at(event_type_name) << ")"
                 << std::endl;
  }
  // Impact table:
  log.indent() << "- Impact number (occurring type -> impacted type) "
                  "(Avg. # over allowed events):"
               << std::endl;
  for (auto const &pair : x.equiv_keys_by_type) {
    std::string event_type_name = pair.first;
    log.indent() << "  - " << event_type_name << " = ";
    for (auto const &pair2 : x.impact_table.by_type.at(event_type_name)) {
      log.indent() << pair2.second / x.n_allowed.by_type.at(event_type_name)
                   << " ";
    }
    log.indent() << std::endl;

    for (auto const &equiv_key : pair.second) {
      Index equivalent_index = equiv_key.second;
      log.indent() << "    - " << event_type_name << "." << equivalent_index
                   << " = ";
      for (auto const &pair2 :
           x.impact_table.by_equivalent_index.at(equiv_key)) {
        log.indent() << pair2.second /
                            x.n_allowed.by_equivalent_index.at(equiv_key)
                     << " ";
      }
      log.indent() << std::endl;
    }
  }
  log.indent() << std::endl;
  log.end_section();
}

}  // namespace clexmonte

}  // namespace CASM

#endif
