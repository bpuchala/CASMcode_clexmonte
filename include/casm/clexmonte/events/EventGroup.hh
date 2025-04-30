#ifndef CASM_clexmonte_events_EventGroup
#define CASM_clexmonte_events_EventGroup

#include <set>
#include <vector>

#include "casm/casm_io/Log.hh"
#include "casm/clexmonte/definitions.hh"
// #include "casm/clexmonte/events/event_data.hh"
#include "casm/clexmonte/events/lotto.hh"
#include "casm/clexmonte/misc/MapLike.hh"
// #include "casm/clexmonte/state/Configuration.hh"
// #include "casm/monte/events/OccLocation.hh"

namespace CASM {
namespace clexmonte {
namespace event_group {

// EventGroup overview:
// - Stores data needed for saving and restoring states and events for
//   KMC acceleration using passage time analysis
//
// Include an event in the group if:
// - it is a transition between states
// - it is an allowed event that is impacted by a transition between states
//
// Include a site in the group if:
// - it is modified by a transition between states
//
// If removing a state:
// - Remove known transitions to that state
// - Remove events that are no longer impacted by any transition
// - Remove sites that are no longer modified by any transition between states
//
// EventGroupManager should be the main (only?) interface users use to interact
// with the EventGroup, states, etc.
//
// Conventions:
// - Delete should be "fast", just indicating that an item in a container is no
//   longer in use as quickly as possible (may be `clear`, may be setting a
//   value to indicate the item is no longer assigned and ready to be assigned
//   when a new item is needed).
// - When adding an item, should check if container needs expanding and set
//   all item data to valid.

typedef MapLikeTraits::size_type size_type;
typedef MapLikeTraits::key_type key_type;

/// \brief Event info
template <typename SavedEventDataType>
struct EventHolder {
  /// \brief Index in AllowedEventMap
  Index global_event_index;

  /// \brief Event data in each state, by state key
  std::vector<SavedEventDataType> data;

  /// \brief Keys of group transitions this site is modified by
  std::set<key_type> transitions;

  // -- Methods --

  /// \brief Default constructor
  EventHolder() : global_event_index(-1) {}

  /// \brief Reset the event, clearing members
  void reset() {
    global_event_index = -1;
    data.clear();
    transitions.clear();
  }
};

/// \brief Site info
template <typename SavedSiteDataType>
struct SiteHolder {
  /// \brief Linear site index of site
  Index linear_site_index;

  /// \brief Site data, by state key
  std::vector<SavedSiteDataType> data;

  /// \brief Keys of group transitions this site is modified by
  std::set<key_type> transitions;

  // -- Methods --

  /// \brief Default constructor, unassigned
  SiteHolder() : linear_site_index(-1) {}

  /// \brief Reset the site so it is unassigned
  void reset() {
    linear_site_index = -1;
    data.clear();
    transitions.clear();
  }
};

/// \brief Some state info
///
/// - State event info is stored in Event
/// - State site info is stored in Site
///
template <typename SavedStateDataType>
struct StateHolder {
  /// \brief State data, excluding state event or site data
  SavedStateDataType data;

  /// \brief Markov chain the state is included in (as a transient state)
  key_type chain_key;

  /// \brief Index of this state in the Markov chain
  Index chain_state_index;

  /// \brief Keys of group transitions this state is connected by
  std::set<key_type> transitions;

  // -- Methods --

  /// \brief Constructor
  StateHolder() : chain_key(0), chain_state_index(-1) {}

  /// \brief Reset the state so it is unassigned
  void reset() {
    data.reset();
    chain_key = -1;
    chain_state_index = -1;
    transitions.clear();
  }
};

/// \brief Transitions between known states
///
/// - The event, initial state, and final state are saved to allow building
///   an absorbing Markov chain transition matrix
/// - The impacted events and sites are saved so when a transition is removed,
///   the impacted events and sites can be removed from the event group if
///   they are no longer impacted by any other transitions
struct Transition {
  key_type event_key;
  key_type initial_state;
  key_type final_state;
  std::vector<Index> impacted_events;
  std::vector<Index> impacted_sites;

  Transition() : event_key(0), initial_state(0), final_state(0) {}

  void reset() {
    event_key = 0;
    initial_state = 0;
    final_state = 0;
    impacted_events.clear();
    impacted_sites.clear();
  }
};

// struct TransientState {
//   /// \brief Group state index of this state
//   Index group_state_index;
//
//   /// \brief Group event index of events connecting this state to other
//   ///     transient states
//   std::vector<Index> transient_event;
//
//   /// \brief Group state index of connected transient states
//   std::vector<Index> transient_state;
//
//   /// \brief Group event index of events connecting this transient state to
//   ///     known absorbing states
//   std::vector<Index> known_absorbing_event;
//
//   /// \brief Group state index of connected absorbing states
//   std::vector<Index> known_absorbing_state;
//
//   /// \brief Group event index of events connecting with unknown absorbing
//   ///     states
//   std::vector<Index> unknown_absorbing_event;
//
//   // -- Methods --
//
//   /// \brief Default constructor
//   TransientState();
//
//   /// \brief Constructor
//   explicit TransientState(Index _group_state_index);
//
//   void reset();
// };

// struct AbsorbingMarkovChain {
//   /// \brief True if this object is currently assigned a chain
//   bool is_assigned;
//
//   /// \brief Transient states
//   std::vector<TransientState> transient_state;
//
//   /// \brief Rate matrix, M
//   ///
//   /// - Matrix of size (n+1, n+1), where n = n_transient, using one state
//   ///   (index n) to handle all absorbing states
//   /// - M(i,j) = -rate j to i, if i != j; sum_k rate i to k if i = j;
//   ///   using rate n to k = 0 for all k (absorbing state condition)
//   Eigen::MatrixXd rate_matrix;
//
//   /// \brief Rate matrix eigenvalues, \lambda
//   Eigen::VectorXd eigenvalues;
//
//   /// \brief Rate matrix eigenvectors, V
//   Eigen::MatrixXd eigenvectors;
//
//   /// \brief Rate matrix inverse eigenvectors, V^{-1}
//   Eigen::MatrixXd inv_eigenvectors;
//
//   /// \brief Mean residence time in transient states, \tau^{1}
//   ///
//   /// - Vector of size (n,)
//   /// - The mean time spent per visit, for each transient state
//   Eigen::VectorXd mean_residence_time;
//
//   /// \brief Transition probability matrix, T
//   ///
//   /// - Matrix of size (n, n)
//   /// - T(i,j) = (rate j to i ) * (mean residence time in j)
//   Eigen::MatrixXd transition_probability_matrix;
//
//   /// \brief Occupation probability sum matrix, (I - T)^{-1}
//   ///
//   /// - Each column gives the summed occupation probability given that system
//   ///   started in the corresponding state.
//   Eigen::MatrixXd occupation_probability_sum_matrix;
//
//   // -- Methods --
//
//   /// \brief Default constructor
//   AbsorbingMarkovChain() : is_assigned(false) {}
//
//   /// \brief Reset the chain so it is unassigned
//   void reset() {
//     is_assigned = false;
//     transient_state.clear();
//   }
// };

template <typename _SavedEventDataType, typename _SavedSiteDataType,
          typename _SavedStateDataType, typename _StateDataType, bool DebugMode>
struct EventGroup {
  typedef default_engine_type engine_type;

  using SavedEventDataType = _SavedEventDataType;
  using SavedSiteDataType = _SavedSiteDataType;
  using SavedStateDataType = _SavedStateDataType;
  using StateDataType = _StateDataType;

  // -- Group state --

  /// \brief Key for this group
  Index key;

  /// \brief The time at which the event group entered its current state
  monte::TimeType last_time;

  /// \brief The last / latest state of the event group
  key_type last_state;

  /// \brief The previous state of the event group
  ///
  /// - used to set the initial state when adding a transition
  key_type prev_state;

  // -- Data structures for event group --

  /// \brief Events in this group
  MapLike<EventHolder<SavedEventDataType>> event;

  /// \brief Sites with changing occupation between saved states
  MapLike<SiteHolder<SavedSiteDataType>> site;

  /// \brief Saved states
  MapLike<StateHolder<SavedStateDataType>> state;

  /// \brief Transitions between states
  MapLike<Transition> transition;

  // -- Selected event data --

  /// \brief Index in AllowedEventMap of the selected event for this group
  Index next_global_event_index;

  /// \brief The time between the last event and next event for this group
  double next_time_increment;

  /// \brief Index of group state from which the next
  ///     selected event occurs, if applicable
  key_type next_state;

  /// \brief The time the next event occurs for this group
  double next_time;

  // -- Constructor  --

  /// \brief Default constructor
  EventGroup()
      : key(0),
        last_time(0.0),
        last_state(0),
        next_global_event_index(-1),
        next_time_increment(std::numeric_limits<double>::max()),
        next_time(std::numeric_limits<double>::max()),
        next_state(0) {}

  void reset() {
    // -- Group state --
    key = 0;
    last_time = 0.0;
    last_state = 0;

    // -- Data structures for event group --
    event.clear();
    site.clear();
    state.clear();

    // -- Selected event data --
    next_global_event_index = -1;
    next_time_increment = 0.0;
    next_time = 0.0;
    next_state = 0;
  }

  // -- Methods --

  /// \brief Find a state in the event group
  typename MapLike<StateHolder<SavedStateDataType>>::const_iterator find_state(
      StateDataType const &state_data) const;

  /// \brief Evolve the group's state to the given time, under the condition
  ///     that the group has remained in the transient states in the current
  ///     chain
  void resolve_state(monte::TimeType time,
                     lotto::RandomGeneratorT<engine_type> &random_generator);

  /// \brief Select the next event for this group
  void select_next_event(
      lotto::RandomGeneratorT<engine_type> &random_generator);
};

/// \brief Global event info
///
/// - Used to lookup global_event_index -> group_key, group_event_key
/// - Used to hold calculated event data for each allowed event
template <typename EventDataType>
struct GlobalEventHolder {
  /// \brief Index of the event group the event belongs to
  key_type group_key = 0;

  /// \brief Index of the event in the group
  ///
  /// - Note that events are not explicitly added to group 0, so this value is
  ///   arbitrary for group 0
  key_type group_event_key = 0;

  /// \brief Event rate, dE_activated, etc.
  EventDataType data;

  GlobalEventHolder() : group_key(0), group_event_key(0) {}

  void reset() {
    group_key = 0;
    group_event_key = 0;
    data.reset();
  }
};

/// \brief Global site info
///
/// - Used to lookup linear_site_index -> group_key, group_site_key
/// - Does not store site data, which must be stored in the state data
///
struct GlobalSiteHolder {
  key_type group_key = 0;
  key_type group_site_key = 0;

  GlobalSiteHolder() : group_key(0), group_site_key(0) {}

  void reset() {
    group_key = 0;
    group_site_key = 0;
  }
};

/// \brief Interface class for grouping events, saving states, and KMC
///     acceleration
///
/// Usage (high level concept):
/// - Create an EventGroupManager object
/// - By default, all events are ungrouped and stored in group 0
/// - Save calculated event data (rate, dE_activated) for each event either
///   currently allowed or allowed in a saved state
/// - When the system is trapped in a set of states:
///   - Create a new event group
///   - Add to the group all the events that are either (i) transitions between
///     those frequently occurring states or (ii) impacted by those transitions
///   - Also save to the event group all the sites that are modified by those
///     transitions
///   - For each state, save the event data (rate, dE_activated) for the events
///     in the group and site data (occupation, trajectory information) for the
///     occupants on the sites in the group
///   - States may be saved in a single "chain" (absorbing Markov chain) or
///     multiple chains, which may be merged later if it is determined
///     necessary.
/// - Each group will be evolved independently, storing the "last state" the
///   the group entered and the "last time", which is the time at which the
///   group entered that state.
/// - For each group, the next event to occur, excluding saved transitions
///   within the current chain, is determined.
///
/// The main loop:
/// - The next event to occur overall is determined.
/// - For time-based sampling:
///   - While the next sample time is before the next event time:
///     - The state of each group is evolved to the next sample time
///     - Sampling is performed
///     - The next event is selected for each group
///     - The next event to occur overall is determined
/// - Before applying an event, the impact of the selected event is determined:
///   - The selected group is evolved to the state it will be in at the time of
///     the selected event.
///   - Other impacted groups are evolved to the state they will be in at the
///     time of the event.
///   - Event groups are updated (add a new group, add ungrouped events to an
///     existing group, merge groups, etc.) as necessary.
/// - After applying the event:
///   - The transition and new state may be saved for the selected group, or the
///     group may be deleted.
///   - The next event is selected for the selected, new, and impacted groups
/// - For count-based sampling:
///   - The state of each group is evolved to the current time
///   - Sampling is performed
///   - The next event is selected for each group
///   - The next event to occur overall is determined

///
/// \tparam EventDataType Data structure expected to include
///     `EventDataType::rate` and `EventDataType::dE_activated`
/// \tparam SiteDataType Data structure expected to include data necessary to
///     save, compare, and restore states
/// \tparam StateDataType Data structure expected to provide the KMC state
/// \tparam DebugMode If true, print debug information
///
/// - Required functions:
///   - `void save(Index linear_site_index, SavedSiteDataType &site_data,
///     StateDataType const &state_data)` to save the site data
///   - `void restore(Index linear_site_index,
///     SavedSiteDataType const &site_data, StateDataType &state_data)` to
///     restore the site data
///   - `bool is_equal(Index linear_site_index,
///     SavedSiteDataType const &site_data,
///     StateDataType const &state_data)` to compare site data to the current
///     KMC state for purposes of checking if the KMC state is the same as a
///     saved state
template <typename EventGroupType, bool DebugMode>
class EventGroupManager {
 public:
  typedef typename MapLike<EventGroupType>::const_iterator const_iterator;
  using SavedEventData = typename EventGroupType::SavedEventDataType;
  using SavedSiteData = typename EventGroupType::SavedSiteDataType;
  using SavedStateData = typename EventGroupType::SavedStateDataType;
  using StateDataType = typename EventGroupType::StateDataType;

  EventGroupManager();

  // -- SavedEventData methods --

  /// \brief Add ungrouped event with global_event_index to the event data
  ///
  /// - If the specified event is already added, this resets it to be ungrouped
  ///   (in group 0)
  ///
  /// \param global_event_index The global event index
  ///
  /// \raises std::runtime_error if the event is assigned to a group
  void add_global_event_data(Index global_event_index) {
    m_event_data.insert(global_event_index);
    m_event_data[global_event_index].reset();
  }

  /// \brief Delete ungrouped event with global_event_index from the event data
  ///
  /// - This should only be called if the event is ungrouped (in group 0)
  ///
  /// \param global_event_index The global event index
  ///
  /// \raises std::runtime_error if the event is grouped (in a group other than
  ///     0)
  void delete_global_event_data(Index global_event_index) {
    auto const &x = m_event_data[global_event_index];
    if (x.group_key != 0) {
      std::stringstream ss;
      ss << "Error in EventGroupManager::delete_global_event_data: "
         << "Event (=" << global_event_index << ") is assigned to group "
         << x.group_key;
      throw std::runtime_error(ss.str());
    }
    m_event_data.erase(global_event_index);
  }

  /// \brief Clear all global event data
  void clear_global_event_data() { m_event_data.clear(); }

  /// \brief Get the global event data for the specified global event index
  SavedEventData const &global_event_data(Index global_event_index) const {
    return m_event_data[global_event_index].data;
  }

  /// \brief Set the global event data for the specified global event index
  void set_global_event_data(Index global_event_index,
                             SavedEventData const &event_data) {
    m_event_data[global_event_index].data = event_data;
  }

  // -- Site methods --

  void resize_global_site_data(size_type n_sites) {
    m_site_data.resize(n_sites);
  }

  void clear_global_site_data() { m_site_data.clear(); }

  // -- EventGroup methods --

  std::size_t size() const { return m_event_group.size(); }

  const_iterator begin() const { return m_event_group.begin(); }

  const_iterator end() const { return m_event_group.end(); }

  /// \brief The indices of current groups
  std::set<key_type> const &current_groups() const { return m_current_groups; }

  EventGroupType const &operator[](key_type group_key) const {
    return m_event_group[group_key];
  }

  size_type count(key_type group_key) const {
    return m_event_group.count(group_key);
  }

  EventGroupType const &group(key_type group_key) const {
    return m_event_group[group_key];
  }

  /// \brief Create the default group
  void add_default_group(monte::TimeType last_time);

  /// \brief Create a new event group
  key_type add_group(std::vector<Index> const &global_event_indices,
                     monte::TimeType last_time);

  /// \brief Delete an event group and ungroup all its events and sites
  void delete_group(key_type group_key);

  /// \brief Add ungrouped events to the specified event group
  void add_ungrouped_events_to_existing_group(
      std::vector<Index> const &global_event_indices, key_type group_key);

  // -- Add sites, states, and transitions --

  /// \brief Add ungrouped sites to the specified event group
  void add_ungrouped_sites_to_existing_group(
      std::vector<Index> const &linear_site_indices,
      StateDataType const &state_data, key_type group_key);

  /// \brief Check if the current KMC state is an existing saved state
  bool in_existing_state(StateDataType const &state_data, key_type group_key);

  /// \brief Add a state to the event group
  void save_state(StateDataType const &state_data, key_type group_key);

  /// \brief Restore the current KMC state from a saved state
  void restore_state(StateDataType &state_data, key_type state_key,
                     key_type group_key);

  /// \brief Add a transition
  void add_transition(Index transition,
                      std::vector<Index> const &impacted_events,
                      std::vector<Index> const &impacted_sites,
                      key_type group_key);

  // -- Evolution and event selection methods --

  /// \brief This is the next event group selected to occur
  key_type next_event_group() const { return m_next_event_group; }

  /// \brief Set `next_event_group` by finding which group moves next
  void set_next_event_group();

 private:
  EventGroupType &_group(key_type group_key) {
    return m_event_group[group_key];
  }

  /// \brief Event groups for state saving and first passage time analysis
  MapLike<EventGroupType> m_event_group;

  /// \brief Lookup for global event ID to group index and group event index
  MapLike<GlobalEventHolder<SavedEventData>> m_event_data;

  /// \brief Lookup for global site index to group index and group site index
  std::vector<GlobalSiteHolder> m_site_data;

  /// \brief The current event groups
  std::set<key_type> m_current_groups;

  /// \brief If `use_event_groups` is true, this is the next event group
  ///     selected to occur
  key_type m_next_event_group;
};

// -- Implementations --

// -- EventGroup --

/// \brief Find a state in the event group
template <typename _SavedEventDataType, typename _SavedSiteDataType,
          typename _SavedStateDataType, typename _StateDataType,
          bool _DebugMode>
typename MapLike<StateHolder<_SavedStateDataType>>::const_iterator
EventGroup<_SavedEventDataType, _SavedSiteDataType, _SavedStateDataType,
           _StateDataType, _DebugMode>::find_state(StateDataType const
                                                       &state_data) const {
  auto _is_equal_state = [&](key_type state_key,
                             StateDataType const &state_data) {
    for (auto const &s : this->site) {
      if (!is_equal(s.linear_site_index, s.data[state_key], state_data)) {
        return false;
      }
    }
    return true;
  };

  auto it = this->state.begin();
  auto end = this->state.end();
  for (; it != end; ++it) {
    if (_is_equal_state(it.key(), state_data)) {
      return it;
    }
  }
  return end;
}

/// \brief Evolve the group's state to the given time, under the condition
///     that the group has remained in the transient states in the current
///     chain
template <typename _SavedEventDataType, typename _SavedSiteDataType,
          typename _SavedStateDataType, typename _StateDataType,
          bool _DebugMode>
void EventGroup<_SavedEventDataType, _SavedSiteDataType, _SavedStateDataType,
                _StateDataType,
                _DebugMode>::resolve_state(monte::TimeType time,
                                           lotto::RandomGeneratorT<engine_type>
                                               &random_generator) {
  if (true) {
    // if constexpr (DebugMode) {
    Log &log = CASM::log();
    log.custom("Resolve state");
    log.indent() << "- Event group: " << this->key << std::endl;
    log.indent() << "- Resolve state at time: " << time << std::endl;
    log.indent() << "- Current time: " << this->last_time << std::endl;
  }

  // TODO

  if (true) {
    // if constexpr (DebugMode) {
    Log &log = CASM::log();
    log.indent() << "- Resolve state DONE" << std::endl << std::endl;
    log.end_section();
  }
}

// \brief Select the next event for this group
template <typename _SavedEventDataType, typename _SavedSiteDataType,
          typename _SavedStateDataType, typename _StateDataType,
          bool _DebugMode>
void EventGroup<_SavedEventDataType, _SavedSiteDataType, _SavedStateDataType,
                _StateDataType, _DebugMode>::
    select_next_event(lotto::RandomGeneratorT<engine_type> &random_generator) {
  size_type tsize = this->event.size();
  size_type i;
  static std::vector<double> m_tsum;

  // Sum the rates of all assigned events
  m_tsum.resize(tsize + 1);
  m_tsum[0] = 0.;
  i = 0;
  for (auto const &e : this->event) {
    m_tsum[i + 1] = m_tsum[i] + e.data[this->last_state].rate;
    ++i;
  }

  if (m_tsum.back() == 0.0) {
    this->next_state = -1;
    this->next_global_event_index = -1;
    this->next_time = std::numeric_limits<double>::max();
    return;
  }

  // (0, m_tsum.back()]
  double rand = random_generator.sample_unit_interval() * m_tsum.back();

  // Select event
  this->next_global_event_index = -1;
  i = 0;
  for (auto const &e : this->event) {
    if (rand < m_tsum[i + 1]) {
      this->next_global_event_index = this->event[i].global_event_index;
      break;
    }
    ++i;
  }

  this->next_state = -1;
  this->next_time_increment =
      -std::log(random_generator.sample_unit_interval()) / m_tsum.back();
  this->next_time = this->last_time + this->next_time_increment;
}

// -- EventGroupManager --

template <typename EventGroupType, bool DebugMode>
EventGroupManager<EventGroupType, DebugMode>::EventGroupManager()
    : m_next_event_group(0) {}

/// \brief Create the default group
template <typename EventGroupType, bool DebugMode>
void EventGroupManager<EventGroupType, DebugMode>::add_default_group(
    monte::TimeType last_time) {
  m_event_group.insert(0);
  auto &g = _group(0);
  g.reset();
  g.key = 0;
  g.last_time = last_time;
}

/// \brief Create a new event group
///
/// \param allowed_event_map The AllowedEventMap, which tracks all events
///     allowed in any group, which group they belong to, and their index in
///     the group
/// \param global_event_indices The indices of events in the AllowedEventMap
/// \param last_time The time for the new group
/// \return group_key, The key for the new group
template <typename EventGroupType, bool DebugMode>
key_type EventGroupManager<EventGroupType, DebugMode>::add_group(
    std::vector<Index> const &global_event_indices, monte::TimeType last_time) {
  key_type group_key = m_event_group.insert();
  m_current_groups.insert(group_key);
  auto &g = _group(group_key);
  g.reset();
  g.key = group_key;
  g.last_time = last_time;

  add_ungrouped_events_to_existing_group(global_event_indices, group_key);

  return group_key;
}

/// \brief Delete an event group and ungroup all its events and sites
///
/// - This also moves all events in the group to group 0
/// - This also moves all sites in the group to group 0
///
/// \tparam DebugMode
/// \param group_key, The key of the group to delete
template <typename EventGroupType, bool DebugMode>
void EventGroupManager<EventGroupType, DebugMode>::delete_group(
    key_type group_key) {
  if (group_key == 0) {
    throw std::runtime_error(
        "Error in EventGroupManager::delete_group: "
        "Cannot delete group 0");
  }
  if (!m_event_group.count(group_key)) {
    throw std::runtime_error(
        "Error in EventGroupManager::delete_group: "
        "group (=" +
        std::to_string(group_key) + ") does not exist");
  }

  auto &g = _group(group_key);
  for (auto &e : g.event) {
    auto &x = m_event_data[e.global_event_index];
    x.group_key = 0;
    x.group_event_key = 0;
  }
  for (auto &s : g.site) {
    auto &x = m_site_data[s.linear_site_index];
    x.group_key = 0;
    x.group_site_key = 0;
  }

  m_event_group.erase(group_key);
  m_current_groups.erase(group_key);
}

/// \brief Check if the current KMC state is an existing saved state
///
/// - If true, this will also update the group's "last_state" and "prev_state"
///
/// \param state_data
/// \param group_key
/// \return true if the current KMC state is an existing saved state, false
/// otherwise
template <typename EventGroupType, bool DebugMode>
bool EventGroupManager<EventGroupType, DebugMode>::in_existing_state(
    StateDataType const &state_data, key_type group_key) {
  auto &g = group(group_key);
  auto it = g.find_state(state_data);
  if (it != g.state.end()) {
    // update "last_state" and "prev_state"
    g.prev_state = g.last_state;
    g.last_state = it.key();
    return true;
  }
  return false;
}

/// \brief Add a state to the event group
///
/// - This will also update the group's "last_state" and "prev_state"
///
/// \param state_data The current KMC state
/// \param group_key The key of the group the state should be added to
template <typename EventGroupType, bool DebugMode>
void EventGroupManager<EventGroupType, DebugMode>::save_state(
    StateDataType const &state_data, key_type group_key) {
  auto &g = _group(group_key);

  key_type state_key = g.state.insert();

  // set state data
  auto &s = g.state[state_key];
  s.reset();
  s.data.reset();  // TODO

  // save the state's event data
  for (auto &event : g.event) {
    auto &x = m_event_data[event.global_event_index];
    if (event.data.size() <= state_key) {
      event.data.resize(state_key + 1);
    }
    event.data[state_key] = x.data;
  }

  // save the state's site data
  for (auto &site : g.site) {
    if (site.data.size() <= state_key) {
      site.data.resize(state_key + 1);
    }
    save(site.linear_site_index, site.data[state_key], state_data);
  }

  // update "last_state" and "prev_state"
  g.prev_state = g.last_state;
  g.last_state = state_key;
}

/// \brief Restore the current KMC state from a saved state
///
/// - This will also modify the global event data to set it from the saved state
/// event data
/// - This will also update the group's "last_state" and "prev_state"
///
/// \param state_data The current KMC state to be modified
/// \param state_key The key of the saved state to restore
/// \param group_key The key of the group to restore the state from
template <typename EventGroupType, bool DebugMode>
void EventGroupManager<EventGroupType, DebugMode>::restore_state(
    StateDataType &state_data, key_type state_key, key_type group_key) {
  auto &g = _group(group_key);
  auto const &s = g.state[state_key];

  // restore the site data
  for (auto const &site : g.site) {
    restore(site.linear_site_index, s.data[state_key], state_data);
  }

  // restore the event data
  for (auto const &event : g.event) {
    auto const &e = event.data[state_key];
    auto &x = m_event_data[event.global_event_index];
    x.data = e.data[state_key];
  }

  // set the "last state" of the event group
  g.prev_state = g.last_state;
  g.last_state = state_key;
}

/// \brief Add events to the specified event group
///
/// - This adds the event's current rate and dE_activated to all existing saved
/// states
///   in the group.
///
/// \param global_event_indices The indices of events to be added
/// \param group_key The key of the group the events should be added to
///
/// \raises std::runtime_error if the event is already assigned to a group that
///     is not group 0
///
template <typename EventGroupType, bool DebugMode>
void EventGroupManager<EventGroupType, DebugMode>::
    add_ungrouped_events_to_existing_group(
        std::vector<Index> const &global_event_indices, key_type group_key) {
  auto &g = _group(group_key);

  size_type state_size_allocated = g.state.size_allocated();

  for (Index global_event_index : global_event_indices) {
    auto &x = m_event_data[global_event_index];
    if (x.group_key != 0) {
      std::stringstream ss;
      ss << "Error in "
            "EventGroupManager::add_ungrouped_events_to_existing_group: "
         << "Event (=" << global_event_index
         << ") is already assigned to group " << x.group_key;
      throw std::runtime_error(ss.str());
    }

    // assign an event object in the group
    key_type group_event_key = g.event.insert();
    auto &new_group_event = g.event[group_event_key];
    new_group_event.global_event_index = global_event_index;

    // update the global event's group indices
    x.group_key = group_key;
    x.group_event_key = group_event_key;

    // set event data for the new event, for all states, using global event data
    new_group_event.data.resize(g.state.size_allocated());
    auto state_it = g.state.begin();
    auto state_end = g.state.end();
    for (; state_it != state_end; ++state_it) {
      new_group_event.data[state_it.key()] = x.data;
    }
  }
}

/// \brief Add ungrouped sites to the specified event group
///
/// \param linear_site_indices Sites to add
/// \param state_data The current KMC state
/// \param group_key The group to add the sites to
///
/// \raises std::runtime_error if the site is already assigned to a group that
///     is not group 0
///
template <typename EventGroupType, bool DebugMode>
void EventGroupManager<EventGroupType, DebugMode>::
    add_ungrouped_sites_to_existing_group(
        std::vector<Index> const &linear_site_indices,
        StateDataType const &state_data, key_type group_key) {
  auto &g = _group(group_key);
  for (Index linear_site_index : linear_site_indices) {
    auto &x = m_site_data[linear_site_index];
    if (x.group_key != 0) {
      std::stringstream ss;
      ss << "Error in "
            "EventGroupManager::add_ungrouped_sites_to_existing_group: "
         << "Site (=" << linear_site_index << ") is already assigned to group "
         << x.group_key;
      throw std::runtime_error(ss.str());
    }

    // assign a site object in the group
    key_type group_site_key = g.site.insert();
    auto &new_group_site = g.site[group_site_key];
    new_group_site.linear_site_index = linear_site_index;

    // update the site's group indices
    x.group_key = group_key;
    x.group_site_key = group_site_key;

    // set site data for the new site, for all states, using global site data
    new_group_site.data.resize(g.site.size_allocated());
    auto state_it = g.state.begin();
    auto state_end = g.state.end();
    for (; state_it != state_end; ++state_it) {
      save(linear_site_index, new_group_site.data[state_it.key()], state_data);
    }
  }
}

/// \brief Add a transition
///
/// - This method will raise if any of the transition, impacted events, or
///   impacted sites are not already in the specified group
///
/// \param transition Global event index of the transition event
/// \param impacted_events The global event indices of events impacted by the
///     transition event
/// \param impacted_sites The linear site indices of sites modified by the
///     transition
/// \param group_key The group that the transition event was selected for
template <typename EventGroupType, bool DebugMode>
void EventGroupManager<EventGroupType, DebugMode>::add_transition(
    Index transition, std::vector<Index> const &impacted_events,
    std::vector<Index> const &impacted_sites, key_type group_key) {
  auto &g = _group(group_key);
  auto &x_global = m_event_data[transition];
  if (x_global.group_key != group_key) {
    throw std::runtime_error(
        "Error in EventGroupManager::add_transition: "
        "transition event (=" +
        std::to_string(transition) + ") is not in group " +
        std::to_string(group_key));
  }

  // add transition
  key_type group_transition_key = g.transition.insert();
  auto &t = g.transition[group_transition_key];
  t.initial_state = g.prev_state;
  t.final_state = g.last_state;
  t.event_key = x_global.group_event_key;
  t.impacted_events = impacted_events;
  t.impacted_sites = impacted_sites;

  // link transition and initial / final states
  g.state[t.initial_state].transitions.insert(group_transition_key);
  g.state[t.final_state].transitions.insert(group_transition_key);

  // link transition and impacted events
  for (Index global_event_index : impacted_events) {
    auto &x = m_event_data[global_event_index];
    if (x.group_key != group_key) {
      throw std::runtime_error(
          "Error in EventGroupManager::add_transition: "
          "Event (=" +
          std::to_string(global_event_index) + ") is not in group " +
          std::to_string(group_key));
    }

    auto &e = g.event[x.group_event_key];
    e.transitions.insert(group_transition_key);
  }

  // link transition and impacted sites
  for (Index linear_site_index : impacted_sites) {
    auto &x = m_site_data[linear_site_index];
    if (x.group_key != group_key) {
      throw std::runtime_error(
          "Error in EventGroupManager::add_transition: "
          "Site (=" +
          std::to_string(linear_site_index) + ") is not in group " +
          std::to_string(group_key));
    }
    auto &s = g.site[x.group_site_key];
    s.transitions.insert(group_transition_key);
  }
}

template <typename EventGroupType, bool DebugMode>
void EventGroupManager<EventGroupType, DebugMode>::set_next_event_group() {
  double min_time;
  auto it = m_event_group.begin();
  auto end = m_event_group.end();
  m_next_event_group = it.key();
  min_time = it->next_time;

  ++it;
  for (; it != end; ++it) {
    if (it->next_time < min_time) {
      min_time = it->next_time;
      m_next_event_group = it.key();
    }
  }
}

}  // namespace event_group
}  // namespace clexmonte
}  // namespace CASM

#endif  // CASM_clexmonte_events_EventGroup
