#ifndef CASM_clexmonte_events_EventGroup
#define CASM_clexmonte_events_EventGroup

#include <vector>

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/events/event_data.hh"
#include "casm/clexmonte/events/lotto.hh"

namespace CASM {
namespace clexmonte {

struct StateData;

namespace event_group_v2 {

// Include an event in the group if:
// - it is a transition between states
// - it is an allowed event that is impacted by a transition between states

// Include a site in the group if:
// - it is modified by a transition between states

// If adding a transition:
// -

// If removing a state:
// - Remove known transitions to that state
// - Remove events that are no longer impacted by any transition
// - Remove sites that are no longer modified by any transition between states

// EventGroup:
// -> GroupEvent[]
// -> GroupSite[]
// -> GroupState[]
// -> GroupTransition[]
// GroupEvent:
// -> global_event_index
// -> GroupState[] is_allowed
// -> GroupSite[] modified_sites
// GroupState:
// -> {GroupEvent -> rate, dE_activated)} (allowed only)
// -> {GroupSite -> occ, traj}
// GroupSite:
// -> linear_site_index
// GroupTransition:
// -> GroupState initial_state;
// -> GroupState final_state;
// -> GroupEvent forward_event;
// -> GroupEvent reverse_event;
// -> GroupEvent[] impacted;
// -> GroupSite[] modified;
// TransientState:
// -> Index state;
// -> GroupEvent[] transient_event;
// -> GroupState[] transient_state_final;
// -> GroupEvent[] known_absorbing_event;
// -> GroupState[] known_absorbing_state_final;
// -> GroupEvent[] unknown_absorbing_event;
// AbsorbingMarkovChain:
// -> TransientState[] transient_state;

// KnownTransitions:
// -> GroupState initial
// -> GroupState final
// -> GroupEvent forward
// -> GroupEvent reverse
// -> GroupEvent[] impacted
// -> GroupSite[] modified

// event.assign(global_event_index)

}  // namespace event_group_v2

namespace event_group {

/// \brief Event info
struct Event {
  /// \brief True if this object is currently assigned an event
  bool is_assigned;

  /// \brief Index in AllowedEventMap
  Index global_event_index;

  /// \brief Rate of the event in each state
  ///
  /// rate[i_state] -> rate in that state
  std::vector<double> rate;

  /// \brief Energy of the event barrier in each state
  ///
  /// state_dE_activated[i_state] -> dE_activated in that state
  std::vector<double> dE_activated;

  /// \brief Index of transitions this event's state is modified by
  std::set<Index> transitions;

  /// \brief Default constructor, not assigned
  Event() : is_assigned(false), global_event_index(-1) {}

  /// \brief Constructor, assigned but no current rate or dE_activated known
  explicit Event(Index _global_event_index)
      : is_assigned(true), global_event_index(_global_event_index) {}

  void reset() {
    is_assigned = false;
    global_event_index = -1;
    rate.clear();
    dE_activated.clear();
    transitions.clear();
  }

  void set(Index state_index, double rate, double dE_activated) {
    if (state_index >= rate.size()) {
      rate.resize(state_index + 1);
      dE_activated.resize(state_index + 1);
    }
    rate[state_index] = rate;
    dE_activated[state_index] = dE_activated;
  }
};

/// \brief Site info
struct Site {
  /// \brief True if this object is currently assigned a site
  bool is_assigned;

  /// \brief Linear site index of site
  Index linear_site_index;

  /// \brief Occupation of the site in each state
  std::vector<int> occ;

  /// \brief Index of transitions this site is modified by
  std::set<Index> transitions;

  /// \brief Default constructor
  Site() : is_assigned(false), linear_site_index(-1) {}

  /// \brief Constructor
  explicit Site(Index _linear_site_index)
      : is_assigned(true), linear_site_index(_linear_site_index) {}

  /// \brief Constructor
  Site(Index _linear_site_index)
      : is_assigned(true), linear_site_index(_linear_site_index) {}

  void reset() {
    linear_site_index = -1;
    occ.clear();
    transitions.clear();
  }

  void set(Index state_index, int value) {
    if (state_index >= occ.size()) {
      occ.resize(state_index + 1);
    }
    occ[state_index] = value;
  }
};

/// \brief Some state info
///
/// - State event info is stored in Event
/// - State site info is stored in Site
///
struct State {
  /// \brief State energy / generalized enthalpy (relative to state 0)
  double dE;

  /// \brief State sum of all event rates
  ///
  /// - rate_sum -> rate sum over all EventGroup events
  /// - Must be updated when additional events are added to the group
  double rate_sum;

  /// \brief Markov chain the state is included in (as a transient state)
  Index chain;

  /// \brief Index of this state in the Markov chain
  Index chain_index;

  /// \brief Constructor
  State() : dE(0.0), rate_sum(0.0), chain(-1), chain_index(-1) {}

  void reset() {
    dE = 0.0;
    //    occ.clear();
    //    rate.clear();
    //    dE_activated.clear();
    rate_sum = 0.0;
    chain = -1;
    chain_index = -1;
  }
};

/// \brief Transitions between known states
struct Transition {
  /// \brief True if this object is currently assigned a transition
  bool is_assigned;

  /// \brief Index of initial state
  Index initial_state;

  /// \brief Index of final state
  Index final_state;

  /// \brief Index of event that transforms the initial state to the final state
  Index forward_event;

  /// \brief Index of event that transforms the final state to the initial state
  Index reverse_event;

  /// \brief The other events that are impacted by this transition
  ///
  /// - Should only be events are allowed in either the initial or final state
  std::vector<Index> impacted_events;

  /// \brief The sites that are modified by this transition
  std::vector<Index> impacted_sites;
};

struct TransientState {
  /// \brief Group state index of this state
  Index group_state_index;

  /// \brief Group event index of events connecting this state to other
  ///     transient states
  std::vector<Index> transient_event;

  /// \brief Group state index of connected transient states
  std::vector<Index> transient_state;

  /// \brief Group event index of events connecting this transient state to
  ///     known absorbing states
  std::vector<Index> known_absorbing_event;

  /// \brief Group state index of connected absorbing states
  std::vector<Index> known_absorbing_state;

  /// \brief Group event index of events connecting with unknown absorbing
  ///     states
  std::vector<Index> unknown_absorbing_event;

  /// \brief Default constructor
  TransientState() : group_state_index(-1) {}
};

struct AbsorbingMarkovChain {
  /// \brief Transient states
  std::vector<TransientState> transient_state;

  /// \brief Rate matrix, M
  ///
  /// - Matrix of size (n+1, n+1), where n = n_transient, using one state
  ///   (index n) to handle all absorbing states
  /// - M(i,j) = -rate j to i, if i != j; sum_k rate i to k if i = j;
  ///   using rate n to k = 0 for all k (absorbing state condition)
  Eigen::MatrixXd rate_matrix;

  /// \brief Rate matrix eigenvalues, \lambda
  Eigen::VectorXd eigenvalues;

  /// \brief Rate matrix eigenvectors, V
  Eigen::MatrixXd eigenvectors;

  /// \brief Rate matrix inverse eigenvectors, V^{-1}
  Eigen::MatrixXd inv_eigenvectors;

  /// \brief Mean residence time in transient states, \tau^{1}
  ///
  /// - Vector of size (n,)
  /// - The mean time spent per visit, for each transient state
  Eigen::VectorXd mean_residence_time;

  /// \brief Transition probability matrix, T
  ///
  /// - Matrix of size (n, n)
  /// - T(i,j) = (rate j to i ) * (mean residence time in j)
  Eigen::MatrixXd transition_probability_matrix;

  /// \brief Occupation probability sum matrix, (I - T)^{-1}
  ///
  /// - Each column gives the summed occupation probability given that system
  ///   started in the corresponding state.
  Eigen::MatrixXd occupation_probability_sum_matrix;

  /// \brief Default constructor
  AbsorbingMarkovChain() {}
};

template <bool DebugMode>
struct EventGroup {
  typedef default_engine_type engine_type;

  /// \brief Index of this group
  Index group;

  /// \brief Current state data
  std::shared_ptr<StateData> state_data;

  /// \brief The time at which the event group entered its current state
  monte::TimeType current_time;

  /// \brief Index in `event` where event rates calculated for the current
  ///     state should be stored
  Index current_state;

  // -- Data structures for event group --

  /// \brief Events in this group
  std::vector<Event> event;

  /// \brief Index of available objects in `event`
  std::vector<Index> available_event;

  /// \brief Sites with changing occupation between saved states
  std::vector<Site> site;

  /// \brief Index of available objects in `site`
  std::vector<Index> available_site;

  /// \brief Saved states
  std::vector<State> state;

  /// \brief Index of available objects in `state`
  std::vector<Index> available_state;

  /// \brief Transitions between states
  std::vector<Transition> transition;

  /// \brief Index of available objects in `transition`
  std::vector<Index> available_transition;

  // -- Chains --

  /// \brief Absorbing Markov chains
  std::vector<AbsorbingMarkovChain> chain;

  // -- Selected event data --

  /// \brief Index of group state from which the next
  ///     selected event occurs, if applicable
  Index next_state;

  /// \brief Index in AllowedEventMap of the selected event for this group
  Index next_global_event_index;

  /// \brief The time between the last event and next event for this group
  double next_time_increment;

  /// \brief The time the next event occurs for this group
  double next_time;

  // -- Constructor  --

  /// \brief Constructor
  explicit EventGroup(std::shared_ptr<StateData> _state_data,
                      monte::TimeType _current_time, Index _group);

  // -- Methods --

  /// \brief Save the current state
  void save_state(StateData const& state_data);

  /// \brief Restore a saved state
  void restore_state(Index state_index, StateData& state_data) const;

  /// \brief Evolve the group's state to the given time, under the condition
  ///     that the group has remained in the transient states in the current
  ///     chain
  void resolve_state(monte::TimeType time,
                     lotto::RandomGeneratorT<engine_type>& random_generator);

  /// \brief Select the next event for this group
  void select_next_event(
      lotto::RandomGeneratorT<engine_type>& random_generator);
};

}  // namespace event_group
}  // namespace clexmonte
}  // namespace CASM

#endif  // CASM_clexmonte_events_EventGroup
