#ifndef CASM_clexmonte_events_state_graph
#define CASM_clexmonte_events_state_graph

#include <map>
#include <optional>
#include <queue>
#include <vector>

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/events/event_data.hh"

namespace CASM {
namespace clexmonte {

class AllowedEventMap;

namespace state_graph {

enum class StateSelectionMethod {
  BASIN_JUMP_FIRST,
  N_JUMP_FIRST,
  BASIN_LOOK_FIRST,
  N_LOOK_FIRST,
};

enum class EventSelectionMethod {
  FIRST_PASSAGE_TIME_ANALYSIS,
  MEAN_RATE_METHOD,
  MEMORY_ONLY,
};

struct Options {
  Options();

  /// \brief Choice of method for saving states in the state graph
  /// (default=N_JUMP_FIRST)
  StateSelectionMethod state_selection_method;

  /// \brief Choice of method for selecting events (default=MEMORY_ONLY)
  EventSelectionMethod event_selection_method;

  /// \brief Number of recent events to track (default=100)
  ///
  /// If the number of unique events in the recent events queue is less than
  /// `start_saving_unique_recent_events_frac` of the total number of recent
  /// events, then states will start being saved.
  Index n_recent_events;

  /// \brief Number of recent states to save
  ///
  /// Applicable for N_JUMP_FIRST and N_LOOK_FIRST methods
  Index n_states;

  //  /// \brief Start saving states when the fraction of recent events which
  //  are
  //  /// unique is less than this threshold (default=0.1)
  //  std::optional<double> start_saving_unique_recent_events_frac;
  //
  //  /// \brief Start saving states when the energy (per unitcell) is less than
  //  /// this threshold
  //  std::optional<double> start_saving_energy_per_unitcell;
  //
  //  /// \brief Stop saving threshold (default=100.0)
  //  ///
  //  /// If the sum of the occupation probabilities in the transient states
  //  over
  //  /// all possible numbers of jumps is less than this, then stop saving
  //  states. std::optional<double> stop_saving_occ_prob_sum;
  //
  //  /// \brief Stop saving states when the energy (per unitcell) is greater
  //  than
  //  /// this threshold
  //  std::optional<double> stop_saving_energy_per_unitcell;
  //
  //  /// \brief Equilibrating basin threshold (default = 1e7)
  //  ///
  //  /// The equilibrating basin approximation is only used in conjunction with
  //  /// FPTA (1) after there is an ill-conditioning error, and (2) when we
  //  find a
  //  /// set of states between which transition rates are this many times
  //  faster
  //  /// than any transition rate to a state outside of the set.
  //  std::optional<double> equilibrating_basin_threshold;
};

/// \brief A StateGraph reference state
///
/// This structure holds the reference state for the state graph. It holds the
/// energy of the first saved state, and the occupation and event rates for
/// those sites and events which change when the system transition between
/// saved states.
struct ReferenceState {
  ReferenceState();

  /// \brief The energy (per supercell) of the reference state.
  double energy_per_supercell;

  /// \brief Indices of sites that change as the system transitions between
  /// saved states
  std::vector<Index> linear_site_index;

  /// \brief Occupation of the sites in `linear_site_index` at the time this
  /// reference state was saved.
  std::vector<int> occ;

  std::vector<EventID> event_id;
  std::vector<double> rate;

  void reset();
};

/// \brief A StateGraph state
///
/// Represents a state in the state graph, which includes the energy,
/// occupation, and rates of events in the configuration at the time the state
/// was saved.
struct State {
  State();

  /// \brief The energy (per supercell) in this state, relative to the
  /// reference state
  double dE;

  /// \brief The occupation of the configuration in this state for the
  /// sites stored in `ReferenceState.linear_site_index` at the time this
  /// state is saved.
  std::vector<int> occ;

  /// \brief The rate of events in this configuration for the events stored
  /// in `ReferenceState.event_id` at the time this state is saved.
  std::vector<double> rate;

  /// \brief The total rate of events in this state.
  ///
  /// Note that this includes both the event rates explicitly stored in `rate`
  /// *and* those not explicitly stored in `rate`
  double total_rate;
};

/// \brief Return true if current configuration `config` and `state` have the
/// same occupation; false otherwise
bool is_equal(config_type const &config, State const &state,
              ReferenceState const &reference);

/// \brief A StateGraph edge (transition) between two states
struct Edge {
  Edge();

  int state_init;
  int state_final;
  double rate_init_to_final;
  double rate_final_to_init;
};

// struct MultiStateGraph {
//   std::vector<StateGraph> state_graphs;
//   std::map<EventID, Index> event_id_to_state_graph_index;
//
//   void clear() {
//     this->state_graphs.clear();
//     this->event_id_to_state_graph_index.clear();
//   }
// };

struct StateGraph {
  // -- Data --

  /// \brief Options
  Options opt;

  /// \brief The reference state
  ReferenceState reference;

  /// \brief The saved states
  std::vector<State> state;

  /// \brief The known transitions between states
  std::vector<Edge> edge;

  /// \brief The current state index
  Index current_state;

  /// \brief Flag to indicate whether states should be saved or not
  bool do_save_states;

  /// \brief Queue to keep track of recent events
  /// (up to size opt.n_recent_events)
  std::queue<EventID> recent_events;

  /// \brief A map to keep track of the count of each EventID in the recent
  /// events queue
  std::map<EventID, Index> recent_events_count;

  // -- Methods --

  StateGraph(Options const &_opt);

  // -- Recent event tracking methods --

  /// \brief Store the most recent event ID in the recent events queue
  void push_recent_event(EventID const &event_id);

  /// \brief Calculate the fraction of recent events which are unique
  double unique_recent_events_frac() const;

  /// \brief Clear recent events queue and count
  void clear_recent_events();

  // -- State saving methods --

  /// \brief Clear saved states, including the reference state
  void clear_states();

  /// \brief Save state
  template <typename EventSelectorType>
  void save_current_state(config_type const &config,
                          SelectedEvent const &selected_event,
                          EventSelectorType const &event_selector,
                          std::optional<Index> previous_state);

  /// \brief Save state
  template <typename EventSelectorType>
  void save_current_state(config_type const &config,
                          SelectedEvent const &selected_event,
                          EventSelectorType const &event_selector,
                          AllowedEventMap const &allowed_event_map,
                          std::optional<Index> previous_state);

  /// \brief Update the reference state *before* applying `selected_event`
  void update_reference_occ(config_type const &config,
                            SelectedEvent const &selected_event);

  /// \brief Update the reference state *before* applying `selected_event`
  template <typename EventSelectorType>
  void update_reference_rate(config_type const &config,
                             SelectedEvent const &selected_event,
                             EventSelectorType const &event_selector);

  /// \brief Update the reference state *before* applying `selected_event`
  template <typename EventSelectorType>
  void update_reference_rate(config_type const &config,
                             SelectedEvent const &selected_event,
                             EventSelectorType const &event_selector,
                             AllowedEventMap const &allowed_event_map);

  /// \brief Find a state equal to `config`
  std::vector<State>::const_iterator find_state(
      config_type const &config) const;

  /// \brief Set the current state
  void set_current_state(Index state_index) {
    this->current_state = state_index;
  }

  /// \brief Get the current state
  Index get_current_state() const { return this->current_state; }
};

}  // namespace state_graph
}  // namespace clexmonte
}  // namespace CASM

#endif  // CASM_clexmonte_events_state_graph
