#ifndef CASM_clexmonte_events_state_graph_impl
#define CASM_clexmonte_events_state_graph_impl

#include "casm/clexmonte/events/AllowedEventList.hh"
#include "casm/clexmonte/events/event_data.hh"
#include "casm/clexmonte/events/lotto.hh"
#include "casm/clexmonte/events/state_graph.hh"
#include "casm/clexmonte/state/Configuration.hh"

namespace CASM {
namespace clexmonte {
namespace state_graph {

/// \brief Save state
///
/// \param config The configuration of the state being added
/// \param selected_event The SelectedEvent from the previous configuration
/// to `config`
template <typename EventSelectorType>
void StateGraph::save_current_state(config_type const &config,
                                    SelectedEvent const &selected_event,
                                    EventSelectorType const &event_selector,
                                    std::optional<Index> previous_state) {
  Index new_state_index = this->state.size();
  this->state.emplace_back();
  State &new_state = this->state.back();

  // -- Save energy --
  if (previous_state.has_value()) {
    if (selected_event.event_state == nullptr) {
      throw std::runtime_error(
          "Error in StateGraph::save_state: "
          "previous_state has value and selected_event.event_state is null");
    }

    EventState const &event_state = *selected_event.event_state;
    new_state.dE = this->state[*previous_state].dE + event_state.dE_final;
  } else {
    new_state.dE = 0.0;
  }

  // -- Save occ --
  // store the occ value on the existing mutated sites from `config`
  Eigen::VectorXi const &occupation = config.dof_values.occupation;
  for (Index l : this->reference.linear_site_index) {
    new_state.occ.push_back(occupation(l));
  }

  // -- Save rates --
  // store the rate for the existing affected events from `config`
  for (EventID const &event_id : this->reference.event_id) {
    new_state.rate.push_back(event_selector->get_rate(event_id));
  }

  // -- Save total rate --
  new_state.total_rate = event_selector.total_rate();

  // -- Add edge --
  if (previous_state.has_value()) {
    EventState const &event_state = *selected_event.event_state;
    this->edge.emplace_back();
    Edge &new_edge = this->edge.back();
    new_edge.state_init = *previous_state;
    new_edge.state_final = new_state_index;
    new_edge.rate_init_to_final = event_state.rate;
    new_edge.rate_final_to_init = event_state.reverse_rate;
  }
}

/// \brief Save state
///
/// \param config The configuration of the state being added
/// \param selected_event The SelectedEvent from the previous configuration
/// to `config`
template <typename EventSelectorType>
void StateGraph::save_current_state(config_type const &config,
                                    SelectedEvent const &selected_event,
                                    EventSelectorType const &event_selector,
                                    AllowedEventMap const &allowed_event_map,
                                    std::optional<Index> previous_state) {
  if (selected_event.event_data == nullptr) {
    throw std::runtime_error(
        "Error in StateGraph::save_state: "
        "selected_event.event_data is null");
  }

  monte::OccEvent const &event = selected_event.event_data->event;

  Index new_state_index = this->state.size();
  this->state.emplace_back();
  State &new_state = this->state.back();

  // -- Save energy --
  if (previous_state.has_value()) {
    EventState const &event_state = *selected_event.event_state;
    new_state.dE = this->state[*previous_state].dE + event_state.dE_final;
  } else {
    new_state.dE = 0.0;
  }

  // -- Save occ --
  // store the occ value on the existing mutated sites from `config`
  Eigen::VectorXi const &occupation = config.dof_values.occupation;
  for (Index l : this->reference.linear_site_index) {
    new_state.occ.push_back(occupation(l));
  }

  // -- Save rates --
  // store the rate for the existing affected events from `config`
  for (EventID const &event_id : this->reference.event_id) {
    auto it = allowed_event_map.find(event_id);
    if (it == allowed_event_map.events().end()) {
      new_state.rate.push_back(0.0);
    } else {
      new_state.rate.push_back(event_selector->get_rate(*it));
    }
  }

  // -- Save total rate --
  new_state.total_rate = event_selector.total_rate();

  // -- Add edge --
  if (previous_state.has_value()) {
    EventState const &event_state = *selected_event.event_state;
    this->edge.emplace_back();
    Edge &new_edge = this->edge.back();
    new_edge.state_init = *previous_state;
    new_edge.state_final = new_state_index;
    new_edge.rate_init_to_final = event_state.rate;
    new_edge.rate_final_to_init = event_state.reverse_rate;
  }
}

/// \brief Update the reference state *before* applying `selected_event`
template <typename EventSelectorType>
void StateGraph::update_reference_rate(
    config_type const &config, SelectedEvent const &selected_event,
    EventSelectorType const &event_selector) {
  std::vector<EventID> &impacted_events = event_selector.get_impacted_events();
  for (EventID const &id : impacted_events) {
    bool already_existing = false;
    for (EventID const &id_existing : this->reference.event_id) {
      if (id == id_existing) {
        already_existing = true;
        break;
      }
    }
    if (!already_existing) {
      this->reference.event_id.push_back(id);
      this->reference.rate.push_back(event_selector.get_rate(id));
    }
  }
}

template <typename EventSelectorType>
void StateGraph::update_reference_rate(
    config_type const &config, SelectedEvent const &selected_event,
    EventSelectorType const &event_selector,
    AllowedEventMap const &allowed_event_map) {
  std::vector<Index> &impacted_events = event_selector.get_impacted_events();
  for (Index event_index : impacted_events) {
    EventID const &id = allowed_event_map.event_id(event_index);
    bool already_existing = false;
    for (EventID const &id_existing : this->reference.event_id) {
      if (id == id_existing) {
        already_existing = true;
        break;
      }
    }
    if (!already_existing) {
      this->reference.event_id.push_back(id);
      this->reference.rate.push_back(event_selector.get_rate(event_index));
    }
  }
}

}  // namespace state_graph
}  // namespace clexmonte
}  // namespace CASM

#endif  // CASM_clexmonte_events_state_graph_impl