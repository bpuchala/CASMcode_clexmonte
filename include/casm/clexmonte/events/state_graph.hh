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

// TODO:
//
// Possible options:
// - Track recent events and start saving states when the fraction of recent
//   events which are unique is less than some threshold (default=0.1)
// -

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

}  // namespace state_graph
}  // namespace clexmonte
}  // namespace CASM

#endif  // CASM_clexmonte_events_state_graph
