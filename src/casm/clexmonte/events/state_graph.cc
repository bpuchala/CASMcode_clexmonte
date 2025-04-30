#include "casm/clexmonte/events/state_graph.hh"

namespace CASM {
namespace clexmonte {
namespace state_graph {

Options::Options()
    : state_selection_method(StateSelectionMethod::N_JUMP_FIRST),
      event_selection_method(EventSelectionMethod::MEMORY_ONLY),
      n_recent_events(100),
      n_states(20) {}

}  // namespace state_graph
}  // namespace clexmonte
}  // namespace CASM
