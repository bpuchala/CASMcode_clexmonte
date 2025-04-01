
#include "casm/clexmonte/events/state_graph.hh"

#include "casm/clexmonte/state/Configuration.hh"

namespace CASM {
namespace clexmonte {
namespace state_graph {

Options::Options()
    : n_recent_events(100),
      start_saving_unique_recent_events_frac(0.1),
      start_saving_energy_per_unitcell(std::nullopt),
      stop_saving_occ_prob_sum(100.0),
      stop_saving_energy_per_unitcell(std::nullopt) {}

ReferenceState::ReferenceState() : energy_per_supercell(0.0) {}

void ReferenceState::reset() {
  energy_per_supercell = 0.0;
  linear_site_index.clear();
  occ.clear();
  event_id.clear();
  rate.clear();
}

State::State() : dE(0.0), total_rate(0.0) {}

/// \brief Return true if current configuration `config` and `state` have the
/// same occupation; false otherwise
bool is_equal(config_type const &config, State const &state,
              ReferenceState const &reference) {
  Index l;
  Index i = 0;
  Eigen::VectorXi const &occupation = config.dof_values.occupation;
  for (; i < state.occ.size(); ++i) {
    l = reference.linear_site_index[i];
    if (occupation(l) != state.occ[i]) {
      return false;
    }
  }
  for (; i < reference.occ.size(); ++i) {
    l = reference.linear_site_index[i];
    if (occupation(l) != reference.occ[i]) {
      return false;
    }
  }
  return true;
}

Edge::Edge()
    : state_init(-1),
      state_final(-1),
      rate_init_to_final(0.0),
      rate_final_to_init(0.0) {}

StateGraph::StateGraph(Options const &_opt)
    : opt(_opt), do_save_states(false) {}

/// \brief Store the most recent event ID in the recent events queue
///
/// This also updates the count of each event ID in the recent events count map,
/// and updates whether states should be saved based on the fraction of recent
/// events that are unique
void StateGraph::push_recent_event(EventID const &event_id) {
  // Add event to recent events queue
  recent_events.push(event_id);

  // Keep track of the count of each event ID in the recent events, use
  // iterators:
  auto it = recent_events_count.find(event_id);
  if (it != recent_events_count.end()) {
    // If the event ID already exists in the count map, increment its count
    it->second++;
  } else {
    // If the event ID does not exist in the count map,
    // add it with a count of 1
    recent_events_count.emplace(event_id, 1);
  }

  // If the queue exceeds the specified size, remove the oldest event
  if (recent_events.size() > this->opt.n_recent_events) {
    // Decrease the count for the oldest event
    auto it_oldest = recent_events_count.find(recent_events.front());
    if (it_oldest != recent_events_count.end()) {
      it_oldest->second--;
      // If the count reaches zero, remove it from the map
      if (it_oldest->second == 0) {
        recent_events_count.erase(it_oldest);
      }
    } else {
      throw std::runtime_error(
          "Error in StateGraph::push: event missing form "
          "recent_events_count");
    }
    recent_events.pop();  // Remove it from the queue
  }

  // Update whether states should be saved based on the fraction of recent
  // events that are unique
  if (do_save_states == false &&
      recent_events.size() == this->opt.n_recent_events) {
    double unique_frac = this->unique_recent_events_frac();
    if (unique_frac < this->opt.start_saving_unique_recent_events_frac) {
      // If the unique fraction is below the threshold, do not save states
      do_save_states = true;
    }
  }
}

/// \brief Calculate the fraction of recent events which are unique
double StateGraph::unique_recent_events_frac() const {
  // Calculate the fraction of unique events in the recent events
  if (this->recent_events.size() == 0) {
    return 0.0;  // Avoid division by zero
  }
  return static_cast<double>(this->recent_events_count.size()) /
         this->recent_events.size();
}

/// \brief Clear recent events queue and count
void StateGraph::clear_recent_events() {
  do_save_states = false;
  recent_events = std::queue<EventID>();  // Clear the queue
  recent_events_count.clear();
}

// -- State saving methods --

/// \brief Clear saved states, including the reference state
void StateGraph::clear_states() {
  this->reference.reset();
  this->state.clear();
  this->edge.clear();
}

/// \brief Update the reference state *before* applying `selected_event`
void StateGraph::update_reference_occ(config_type const &config,
                                      SelectedEvent const &selected_event) {
  if (selected_event.event_data == nullptr) {
    throw std::runtime_error("Error in StateGraph::add: event_data is null");
  }
  monte::OccEvent const &event = selected_event.event_data->event;
  Eigen::VectorXi const &occupation = config.dof_values.occupation;

  // If any new sites are involved in the next selected event,
  // save the site index and current occupant
  for (Index l : event.linear_site_index) {
    bool already_existing = false;
    for (Index l_existing : reference.linear_site_index) {
      if (l == l_existing) {
        already_existing = true;
        break;
      }
    }
    if (!already_existing) {
      reference.linear_site_index.push_back(l);
      reference.occ.push_back(occupation(l));
    }
  }
};

/// \brief Find a state equal to `config`
std::vector<State>::const_iterator StateGraph::find_state(
    config_type const &config) const {
  auto it = this->state.begin();
  auto end = this->state.end();
  for (; it != end; ++it) {
    if (is_equal(config, *it, this->reference)) {
      return it;
    }
  }
  return this->state.end();
}

}  // namespace state_graph
}  // namespace clexmonte
}  // namespace CASM
