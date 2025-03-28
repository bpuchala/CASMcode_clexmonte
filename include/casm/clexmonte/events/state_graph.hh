#ifndef CASM_clexmonte_events_state_graph
#define CASM_clexmonte_events_state_graph

namespace CASM {
namespace clexmonte {
namespace state_graph {

struct ReferenceState {
  double energy;

  std::vector<Index> linear_site_index;
  std::vector<int> occ;

  std::vector<EventID> event_id;
  std::vector<double> rate;
};

struct State {
  double energy;
  std::vector<int> occ;
  std::vector<double> rate;
};

struct Edge {
  int state_init;
  int state_final;
  double dE_activated;
  double rate_init_to_final;
  double rate_final_to_init;
}

struct StateGraph {
  /// \brief The reference state
  ReferenceState reference;

  /// \brief The saved states
  std::vector<State> state;

  /// \brief The known transitions between states
  std::vector<Edge> edge;

  StateGraph() : n_states(0) {}

  Index n_mutated_sites() { return this->reference.linear_site_index.size(); }

  bool is_mutated_site(Index l) {
    for (Index l_check : this->reference.linear_site_index) {
      if (l_check == l) {
        return true;
      }
    }
    return false;
  }

  Index get_mutated_site_index(Index l) {
    Index mutated_site_index = 0;
    for (Index l_check : this->reference.linear_site_index) {
      if (l_check == l) {
        return mutated_site_index;
      }
      ++mutated_site_index;
    }
    return mutated_site_index;
  }

  /// \brief Add state
  ///
  /// \param config The configuration of the state being added
  /// \param selected_event The SelectedEvent from the previous configuration
  /// to `config`
  template <typename EventSelectorType>
  void add_state(std::optional<Index> state_init;
                 config_type const &config, SelectedEvent const &selected_event,
                 AllowedEventMap const &allowed_event_map,
                 EventSelectorType const &event_selector) {
    EventState const &event_state = selected_event.event_state;
    monte::OccEvent const &event = selected_event.event_data->event;

    Index new_state_index = this->state.size();
    this->state.emplace_back();
    new_state = this->state.back();

    // -- energy --
    new_state.energy = 0.0;
    if (state_init.has_value()) {
      new_state.energy =
          this->state[state_init.value()].energy + event_state.dE_final;
    }

    // -- occ --
    // store the occ value on the existing mutated sites from `config`
    Eigen::VectorXi const &occupation = configuration.dof_values.occupation;
    for (Index l : this->reference.linear_site_index) {
      state.occ.push_back(occupation(l));
    }

    // -- rates --
    // store the occ value on the existing mutated sites from `config`
    for (EventID const &event_id : this->reference.event_id) {
      auto it = allowed_event_map.find(event_id);
      if (it == allowed_event_map.end()) {
        state.rate.push_back(0.0);
      } else {
        state.rate.push_back(event_selector->get_rate(*it));
      }
    }

    if (state_init.has_value()) {
      this->edge.emplace_back();
      Edge &new_edge = this->edge.back();
      new_edge.state_init = state_init.value();
      new_edge.state_final = new_state_index;
      new_edge.dE_activated =
          this->state[state_init.value()].energy + event_state.dE_activated;
      new_edge.rate_init_to_final = event_state.rate;
      new_edge.rate_final_to_init;
    }
  }

  void add_edge(Index state_init, Index state_final, config_type const &config,
                SelectedEvent const &selected_event) {}

  bool is_state(Index s, config_type const &config,
                SelectedEvent const &selected_event) {
    Index l;
    Index l_occ;
    for (Occupation const &occupation : this->occ) {
      l = occupation.linear_site_index;
      l_occ = occupation.reference;
      if (s < occupation.state.size()) {
        l_occ = occupation.state[s];
      }
      if (config.occupation[l] != l_occ) {
        return false;
      }
    }
    return true;
  }

  std::pair<bool, Index> find_state(config_type const &config,
                                    SelectedEvent const &selected_event) {
    monte::OccEvent const &event = selected_event.event_data->event;
    Index i_state = 0;

    Index _n_states = this->n_states();
    for (Index i = 0; i < _n_states; ++i) {
      for (Index i_occ = 0; i_occ < this->occ.size(); ++i_occ) {
      }
    }

    return {false, -1};
  }

  void update(std::optional<Index> state_init, config_type const &config,
              SelectedEvent const &selected_event) {
    if (selected_event.event_data == nullptr) {
      throw std::runtime_error("Error in StateGraph::add: event_data is null");
    }
    monte::OccEvent const &event = selected_event.event_data->event;

    if (!state_init.has_value()) {
      this->add_state(config, selected_event);
    } else {
      // check if in an existing state
      bool in_existing_state;
      Index state_index;
      std::tie(in_existing_state, state_index) =
          find_state(config, selected_event);

      if (in_existing_state) {
        // add an edge
        this->add_edge(state_index, state_index, config, selected_event);
      } else {
        // add a new state

        // add an edge
      }
    }
  }
};

void add(StateGraph &state_graph, Index state_init, Index state_final,
         double dE_activated, double rate_init_to_final,
         double rate_final_to_init) {
  Edge edge = {state_init, state_final, dE_activated, rate_init_to_final,
               rate_final_to_init};
  state_graph.edge.push_back(edge);
}

}  // namespace state_graph
}  // namespace clexmonte
}  // namespace CASM
