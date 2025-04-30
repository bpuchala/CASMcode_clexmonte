#include "casm/clexmonte/events/EventGroup.hh"

#include "casm/clexmonte/monte_calculator/StateData.hh"

namespace CASM {
namespace clexmonte {
namespace event_group {

//// -- State --
//
///// \brief Constructor
// State::State()
//     : is_assigned(false), dE(0.0), rate_sum(0.0), chain(-1), chain_index(-1)
//     {}
//
///// \brief Assign the state to a chain
// void State::assign(Index _chain, Index _chain_index) {
//   is_assigned = true;
//   chain = _chain;
//   chain_index = _chain_index;
// }
//
///// \brief Reset the state so it is unassigned
// void State::reset() {
//   is_assigned = false;
//   dE = 0.0;
//   rate_sum = 0.0;
//   chain = -1;
//   chain_index = -1;
// }
//
//// -- Transition --
//
///// \brief Default constructor, unassigned
// Transition::Transition()
//     : is_assigned(false),
//       initial_state(-1),
//       final_state(-1),
//       forward_event(-1),
//       reverse_event(-1) {}
//
///// \brief Constructor, assigned
// Transition::Transition(Index _initial_state, Index _final_state,
//                        Index _forward_event, Index _reverse_event)
//     : is_assigned(true),
//       initial_state(_initial_state),
//       final_state(_final_state),
//       forward_event(_forward_event),
//       reverse_event(_reverse_event) {}
//
///// \brief Assign a transition
// void Transition::assign(Index _initial_state, Index _final_state,
//                         Index _forward_event, Index _reverse_event) {
//   is_assigned = true;
//   initial_state = _initial_state;
//   final_state = _final_state;
//   forward_event = _forward_event;
//   reverse_event = _reverse_event;
// }
//
///// \brief Reset the transition so it is unassigned
// void Transition::reset() {
//   is_assigned = false;
//   initial_state = -1;
//   final_state = -1;
//   forward_event = -1;
//   reverse_event = -1;
//   impacted_events.clear();
//   impacted_sites.clear();
// }
//
//// -- TransientState --
//
///// \brief Default constructor
// TransientState::TransientState() : group_state_index(-1) {}
//
///// \brief Constructor
// TransientState::TransientState(Index _group_state_index)
//     : group_state_index(_group_state_index) {}
//
// void TransientState::reset() {
//   group_state_index = -1;
//   transient_event.clear();
//   transient_state.clear();
//   known_absorbing_event.clear();
//   known_absorbing_state.clear();
//   unknown_absorbing_event.clear();
// }

// -- EventGroup --

///// \brief Default constructor
// template <bool DebugMode>
// EventGroup<DebugMode>::EventGroup()
//     : key(0),
//       last_time(0.0),
//       last_state(-1),
//       next_state(-1),
//       next_global_event_index(-1),
//       next_time_increment(std::numeric_limits<double>::max()),
//       next_time(std::numeric_limits<double>::max()) {}
//
///// \brief Reset the event group
// template <bool DebugMode>
// void EventGroup<DebugMode>::reset() {
//   // -- Group state --
//   key = 0;
//   last_time = 0.0;
//   last_state = -1;
//
//   // -- Data structures for event group --
//   event.clear();
//   site.clear();
//   state.clear();
//   available_state.clear();
//   transition.clear();
//   available_transition.clear();
//
//   // -- Chains --
//   chain.clear();
//   available_chain.clear();
//
//   // -- Selected event data --
//   next_state = -1;
//   next_global_event_index = -1;
//   next_time_increment = 0.0;
//   next_time = 0.0;
// }

// template <bool DebugMode>
// void EventGroup<DebugMode>::throw_if_chain_is_out_of_range(Index chain_index)
// {
//   if (chain_index < 0 || chain_index >= this->chain.size()) {
//     std::stringstream ss;
//     ss << "Error in EventGroup: "
//        << "chain_index (=" << chain_index << ") is out of range [0,"
//        << this->chain.size() << ")";
//     throw std::runtime_error(ss.str());
//   }
// }
//
// template <bool DebugMode>
// void EventGroup<DebugMode>::throw_if_chain_is_not_assigned(Index chain_index)
// {
//   if (!this->chain[chain_index].is_assigned) {
//     std::stringstream ss;
//     ss << "Error in EventGroup: "
//        << "chain_index (=" << chain_index << ") is not assigned";
//     throw std::runtime_error(ss.str());
//   }
// }
//
///// \brief Delete the specified chain and its states
// template <bool DebugMode>
// void EventGroup<DebugMode>::delete_chain(Index chain_index) {
//   if (true) {
//     // if constexpr (DebugMode) {
//     Log &log = CASM::log();
//     log.custom("Delete chain");
//     log.indent() << "- Event group: " << this->key << std::endl;
//     log.indent() << "- Deleting chain: " << chain_index << std::endl;
//   }
//
//   // validate
//   this->throw_if_chain_is_out_of_range(chain_index);
//   this->throw_if_chain_is_not_assigned(chain_index);
//   auto &_chain = this->chain[chain_index];
//
//   // delete chain states
//   if (true) {
//     // if constexpr (DebugMode) {
//     Log &log = CASM::log();
//     log.indent() << "- Deleting chain states... " << std::endl;
//     log.increase_indent();
//   }
//   for (TransientState const &transient_state : _chain.transient_state) {
//     this->throw_if_state_is_out_of_range(transient_state.group_state_index);
//     this->throw_if_state_is_not_assigned(transient_state.group_state_index);
//     this->state[transient_state.group_state_index].reset();
//     this->available_state.push_back(transient_state.group_state_index);
//   }
//   if (true) {
//     // if constexpr (DebugMode) {
//     Log &log = CASM::log();
//     log.indent() << "- Deleting chain states... DONE" << std::endl;
//     log.decrease_indent();
//   }
//
//   // delete the chain
//   this->chain[chain_index].reset();
//   this->available_chain.push_back(chain_index);
//
//   if (true) {
//     // if constexpr (DebugMode) {
//     Log &log = CASM::log();
//     log.indent() << "- Delete chain DONE" << std::endl << std::endl;
//     log.end_section();
//   }
// }

// template <bool DebugMode>
// void EventGroup<DebugMode>::throw_if_state_is_out_of_range(Index state_index)
// {
//   if (state_index < 0 || state_index >= this->state.size()) {
//     std::stringstream ss;
//     ss << "Error in EventGroup: "
//        << "state_index (=" << state_index << ") is out of range [0,"
//        << this->state.size() << ")";
//     throw std::runtime_error(ss.str());
//   }
// }
//
// template <bool DebugMode>
// void EventGroup<DebugMode>::throw_if_state_is_not_assigned(Index state_index)
// {
//   if (!this->state[state_index].is_assigned) {
//     std::stringstream ss;
//     ss << "Error in EventGroup: "
//        << "state_index (=" << state_index << ") is not assigned";
//     throw std::runtime_error(ss.str());
//   }
// }
//
// template <bool DebugMode>
// void EventGroup<DebugMode>::throw_if_group_state_is_out_of_range() {
//   if (this->last_state < 0 || this->last_state >= this->state.size()) {
//     std::stringstream ss;
//     ss << "Error in EventGroup: "
//        << "last_state (=" << this->last_state << ") is out of range [0,"
//        << this->state.size() << ")";
//     throw std::runtime_error(ss.str());
//   }
// }
//
// template <bool DebugMode>
// void EventGroup<DebugMode>::throw_if_group_state_is_not_assigned() {
//   if (!this->state[this->last_state].is_assigned) {
//     std::stringstream ss;
//     ss << "Error in EventGroup: "
//        << "last_state (=" << this->last_state << ") is not assigned";
//     throw std::runtime_error(ss.str());
//   }
// }
//
///// \brief Set the current state to the next available
// template <bool DebugMode>
// void EventGroup<DebugMode>::set_group_state_to_next_available() {
//   if (true) {
//     // if constexpr (DebugMode) {
//     Log &log = CASM::log();
//     log.custom("Set current state to next available");
//     log.indent() << "- Event group: " << this->key << std::endl;
//     log.indent() << "- Original current state: " << this->last_state
//                  << std::endl;
//   }
//
//   if (this->available_state.empty()) {
//     auto &ref = this->state.emplace_back();
//     this->last_state = this->state.size() - 1;
//   } else {
//     this->last_state = this->available_state.back();
//     this->available_state.pop_back();
//   }
//
//   if (true) {
//     // if constexpr (DebugMode) {
//     Log &log = CASM::log();
//     log.indent() << "- New current state: " << this->last_state << std::endl;
//     log.indent() << "- Set current state... DONE" << std::endl << std::endl;
//     log.end_section();
//   }
// }
//
///// \brief Set the current state to the next available and copy the event info
/////
///// - This does copy Event data
///// - This does not copy Site data
// template <bool DebugMode>
// void EventGroup<DebugMode>::copy_group_state_to_next_available() {
//   if (true) {
//     // if constexpr (DebugMode) {
//     Log &log = CASM::log();
//     log.custom("Copy current state to next available");
//     log.indent() << "- Event group: " << this->key << std::endl;
//     log.indent() << "- Original current state: " << this->last_state
//                  << std::endl;
//   }
//
//   // validate
//   this->throw_if_group_state_is_out_of_range();
//   this->throw_if_group_state_is_not_assigned();
//
//   // copy event info from last state to new current state
//   Index last_state = this->last_state;
//   this->set_group_state_to_next_available();
//   if (last_state != -1) {
//     for (Event &e : this->event) {
//       e.copy_state(this->last_state, last_state);
//     }
//   }
//
//   if (true) {
//     // if constexpr (DebugMode) {
//     Log &log = CASM::log();
//     log.indent() << "- New current state: " << this->last_state << std::endl;
//     log.indent() << "- Copy current state... DONE" << std::endl << std::endl;
//     log.end_section();
//   }
// }
//
///// \brief Save the current state
/////
///// - Sets the occupation and atom info for each `Site` in `EventGroup::site`
/////   from the values in `EventGroup::state` and `EventGroup::occ_location`.
// template <bool DebugMode>
// void EventGroup<DebugMode>::save_state(GlobalState const &global_state) {
//   if (true) {
//     // if constexpr (DebugMode) {
//     Log &log = CASM::log();
//     log.custom("Save state");
//     log.indent() << "- Event group: " << this->key << std::endl;
//     log.indent() << "- Saving current state: " << this->last_state <<
//     std::endl;
//   }
//
//   // validate
//   this->throw_if_group_state_is_out_of_range();
//   this->throw_if_group_state_is_not_assigned();
//
//   // save the current state's site info
//   for (Site &_site : this->site) {
//     _site.save(this->last_state, global_state);
//   }
//
//   // save state as new transient state in current chain
//   auto &_state = this->state[this->last_state];
//   Index chain_index = _state.chain_index;
//   this->throw_if_chain_is_out_of_range(chain_index);
//   this->throw_if_chain_is_not_assigned(chain_index);
//   auto &_chain = this->chain[chain_index];
//   _state.assign(chain_index, _chain.transient_state.size());
//
//   // set the transient state attributes...
//   auto &_transient_state = _chain.transient_state.emplace_back();
//   _transient_state.group_state_index = this->last_state;
//
//   if (true) {
//     // if constexpr (DebugMode) {
//     Log &log = CASM::log();
//     log.indent() << "- Save state DONE" << std::endl << std::endl;
//     log.end_section();
//   }
// }
//
///// \brief Delete the specified state
// template <bool DebugMode>
// void EventGroup<DebugMode>::delete_state(Index state_index) {
//   if (true) {
//     // if constexpr (DebugMode) {
//     Log &log = CASM::log();
//     log.custom("Delete state");
//     log.indent() << "- Event group: " << this->key << std::endl;
//     log.indent() << "- Deleting state: " << state_index << std::endl;
//   }
//
//   // validate
//   this->throw_if_state_is_out_of_range(state_index);
//   this->throw_if_state_is_not_assigned(state_index);
//   // delete state
//   this->state[state_index].reset();
//   this->available_state.push_back(state_index);
//
//   if (true) {
//     // if constexpr (DebugMode) {
//     Log &log = CASM::log();
//     log.indent() << "- Delete state DONE" << std::endl << std::endl;
//     log.end_section();
//   }
// }
//
///// \brief Restore a saved state
/////
///// - Sets the occupation and atom info for each `Site` in `EventGroup::site`
/////   to `EventGroup::state` and `EventGroup::occ_location`.
///// - Sets `EventGroup::last_state` to `state_index`.
// template <bool DebugMode>
// void EventGroup<DebugMode>::restore_state(Index state_index) {
//   if (true) {
//     // if constexpr (DebugMode) {
//     Log &log = CASM::log();
//     log.custom("Restore state");
//     log.indent() << "- Event group: " << this->key << std::endl;
//     log.indent() << "- Restore state: " << state_index << std::endl;
//   }
//
//   // Restore the state
//   this->last_state = state_index;
//
//   // Set the state occupation
//   for (Site &_site : this->site) {
//     _site.restore(state_index, this->_occupation(), this->_occ_location());
//   }
//
//   // Set the event rates
//   // TODO
//
//   if (true) {
//     // if constexpr (DebugMode) {
//     Log &log = CASM::log();
//     log.indent() << "- Restore state DONE" << std::endl << std::endl;
//     log.end_section();
//   }
// }

///// \brief Evolve the group's state to the given time, under the condition
/////     that the group has remained in the transient states in the current
/////     chain
// template <typename EventDataType, typename SiteDataType, bool DebugMode>
// void EventGroup<EventDataType, SiteDataType, DebugMode>::resolve_state(
//     monte::TimeType time,
//     lotto::RandomGeneratorT<engine_type> &random_generator) {
//   if (true) {
//     // if constexpr (DebugMode) {
//     Log &log = CASM::log();
//     log.custom("Resolve state");
//     log.indent() << "- Event group: " << this->key << std::endl;
//     log.indent() << "- Resolve state at time: " << time << std::endl;
//     log.indent() << "- Current time: " << this->last_time << std::endl;
//   }
//
//   // TODO
//
//   if (true) {
//     // if constexpr (DebugMode) {
//     Log &log = CASM::log();
//     log.indent() << "- Resolve state DONE" << std::endl << std::endl;
//     log.end_section();
//   }
// }
//
//// \brief Select the next event for this group
// template <typename EventDataType, typename SiteDataType, bool DebugMode>
// void EventGroup<EventDataType, SiteDataType, DebugMode>::select_next_event(
//     lotto::RandomGeneratorT<engine_type> &random_generator) {
//   size_type tsize = this->event.size();
//   size_type i;
//   static std::vector<double> m_tsum;
//
//   // Sum the rates of all assigned events
//   m_tsum.resize(tsize + 1);
//   m_tsum[0] = 0.;
//   i = 0;
//   for (Event const &e : this->event) {
//     m_tsum[i + 1] = m_tsum[i] + e.data[this->last_state].rate;
//     ++i;
//   }
//
//   if (m_tsum.back() == 0.0) {
//     this->next_state = -1;
//     this->next_global_event_index = -1;
//     this->next_time = std::numeric_limits<double>::max();
//     return;
//   }
//
//   // (0, m_tsum.back()]
//   double rand = random_generator.sample_unit_interval() * m_tsum.back();
//
//   // Select event
//   this->next_global_event_index = -1;
//   i = 0;
//   for (Event const &e : this->event) {
//     if (rand < m_tsum[i + 1]) {
//       this->next_global_event_index = this->event[i].global_event_index;
//       break;
//     }
//     ++i;
//   }
//
//   this->next_state = -1;
//   this->next_time_increment =
//       -std::log(random_generator.sample_unit_interval()) / m_tsum.back();
//   this->next_time = this->last_time + this->next_time_increment;
// }
//
// template <typename EventDataType, typename SiteDataType, typename
// StateDataType,
//           bool DebugMode>
// EventGroupManager<EventDataType, SiteDataType, StateDataType,
//                   DebugMode>::EventGroupManager()
//     : m_next_event_group(0) {}
//
///// \brief Create a new event group
/////
///// \param allowed_event_map The AllowedEventMap, which tracks all events
/////     allowed in any group, which group they belong to, and their index in
/////     the group
///// \param global_event_indices The indices of events in the AllowedEventMap
///// \param last_time The time for the new group
///// \return group_key, The key for the new group
// template <typename EventDataType, typename SiteDataType, typename
// StateDataType,
//           bool DebugMode>
// key_type EventGroupManager<EventDataType, SiteDataType, StateDataType,
//                            DebugMode>::add_group(std::vector<Index> const
//                                                      &global_event_indices,
//                                                  monte::TimeType last_time) {
//   key_type group_key = m_event_group.insert();
//   m_current_groups.insert(group_key);
//   auto &g = _group(group_key);
//   g.key = group_key;
//   g.last_time = last_time;
//
//   add_ungrouped_events_to_existing_group(global_event_indices, group_key);
//
//   return group_key;
// }
//
///// \brief Delete an event group
/////
///// - This also moves all events in the group to group 0
/////
///// \tparam DebugMode
///// \param group_key, The key of the group to delete
// template <typename EventDataType, typename SiteDataType, typename
// StateDataType,
//           bool DebugMode>
// void EventGroupManager<EventDataType, SiteDataType, StateDataType,
//                        DebugMode>::delete_group(key_type group_key) {
//   if (group_key == 0) {
//     throw std::runtime_error(
//         "Error in EventGroupManager::delete_group: "
//         "Cannot delete group 0");
//   }
//   if (!m_event_group.count(group_key)) {
//     throw std::runtime_error(
//         "Error in EventGroupManager::delete_group: "
//         "group (=" +
//         std::to_string(group_key) + ") does not exist");
//   }
//
//   for (auto &e : group(group_key).event) {
//     auto &x = m_event_data[e.global_event_index];
//     x.group_key = 0;
//     x.group_event_key = 0;
//   }
//
//   m_event_group.erase(group_key);
//   m_current_groups.erase(group_key);
// }
//
///// \brief Add events to the specified event group
/////
///// - This adds the event's current rate and dE_activated to all existing
/// saved
///// states
/////   in the group.
/////
///// \param global_event_indices The indices of events to be added
///// \param group_key The key of the group the events should be added to
/////
///// \raises std::runtime_error if the event is already assigned to a group
/// that
/////     is not group 0
/////
// template <typename EventDataType, typename SiteDataType, typename
// StateDataType,
//           bool DebugMode>
// void EventGroupManager<EventDataType, SiteDataType, StateDataType,
// DebugMode>::
//     add_ungrouped_events_to_existing_group(
//         std::vector<Index> const &global_event_indices, key_type group_key) {
//   auto &g = _group(group_key);
//   for (Index global_event_index : global_event_indices) {
//     auto &x = m_event_data[global_event_index];
//     if (x.group_key != 0) {
//       std::stringstream ss;
//       ss << "Error in "
//             "EventGroupManager::add_ungrouped_events_to_existing_group: "
//          << "Event (=" << global_event_index
//          << ") is already assigned to group " << x.group_key;
//       throw std::runtime_error(ss.str());
//     }
//
//     // assign an event object in the group
//     key_type group_event_key = g.event.insert();
//     g.event[group_event_key].global_event_index = global_event_index;
//
//     // update the event's group indices
//     x.group_key = group_key;
//     x.group_event_key = group_event_key;
//
//     // update the event data for all states
//     // TODO
//   }
// }
//
// template <typename EventDataType, typename SiteDataType, typename
// StateDataType,
//           bool DebugMode>
// void EventGroupManager<EventDataType, SiteDataType, StateDataType,
//                        DebugMode>::set_next_event_group() {
//   double min_time;
//   auto it = m_event_group.begin();
//   auto end = m_event_group.end();
//   m_next_event_group = it.key();
//   min_time = it->next_time;
//
//   ++it;
//   for (; it != end; ++it) {
//     if (it->next_time < min_time) {
//       min_time = it->next_time;
//       m_next_event_group = it.key();
//     }
//   }
// }

//// DebugMode=false
// template class EventGroup<false>;
// template class EventGroupManager<false>;
//
//// DebugMode=true
// template class EventGroup<true>;
// template class EventGroupManager<true>;

}  // namespace event_group
}  // namespace clexmonte
}  // namespace CASM
