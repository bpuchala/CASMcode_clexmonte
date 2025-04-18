#include "casm/clexmonte/events/EventGroup.hh"

#include "casm/clexmonte/monte_calculator/StateData.hh"

namespace CASM {
namespace clexmonte {
namespace event_group {

/// \brief Constructor
template <bool DebugMode>
EventGroup<DebugMode>::EventGroup(std::shared_ptr<StateData> _state_data,
                                  monte::TimeType _current_time, Index _group)
    : state_data(_state_data),
      current_time(_current_time),
      group(_group),
      current_state(-1),
      next_state(-1),
      next_event_index(-1),
      next_time_increment(std::numeric_limits<double>::max()),
      next_time(std::numeric_limits<double>::max()) {}

/// \brief Save the current state
template <bool DebugMode>
void EventGroup<DebugMode>::save_state() {
  if constexpr (DebugMode) {
    Log &log = CASM::log();
    log.custom("Save state");
    log.indent() << "- Event group: " << this->group << std::endl;
  }

  // Save the state:
  auto const &occupation = this->state_data->state.occupation;
  for (Site &_site : this->site) {
    if (_site.is_assigned) {
      _site.set(this->current_state, occupation(_site.linear_site_index));
    }
  }

  if constexpr (DebugMode) {
    Log &log = CASM::log();
    log.indent() << "- Save state DONE" << std::endl << std::endl;
    log.end_section();
  }
}

/// \brief Restore a saved state
template <bool DebugMode>
void EventGroup<DebugMode>::restore_state(Index state_index) const {
  if constexpr (DebugMode) {
    Log &log = CASM::log();
    log.custom("Restore state");
    log.indent() << "- Event group: " << this->group << std::endl;
    log.indent() << "- Restore state: " << state_index << std::endl;
  }

  // Restore the state
  this->current_state = state_index;

  // Set the state occupation
  auto const &occupation = this->state_data->state.occupation;
  for (Site &_site : this->site) {
    if (_site.is_assigned) {
      occupation(_site.linear_site_index) = _site.occ[state_index];
    }
  }

  if constexpr (DebugMode) {
    Log &log = CASM::log();
    log.indent() << "- Restore state DONE" << std::endl << std::endl;
    log.end_section();
  }
}

/// \brief Evolve the group's state to the given time, under the condition
///     that the group has remained in the transient states in the current
///     chain
template <bool DebugMode>
void EventGroup<DebugMode>::resolve_state(
    monte::TimeType time,
    lotto::RandomGeneratorT<engine_type> &random_generator) {
  Index chain = this->state[this->current_state].chain;
}

// \brief Select the next event for this group
template <bool DebugMode>
void EventGroup<DebugMode>::select_next_event(
    lotto::RandomGeneratorT<engine_type> &random_generator) {
  Index tsize = this->new_state.rate.size();
  static std::vector<double> m_tsum;
  m_tsum.resize(tsize + 1);
  m_tsum[0] = 0.;
  for (Index i = 0; i < tsize; ++i) {
    m_tsum[i + 1] = m_tsum[i] + this->new_state.rate[i];
  }

  if (m_tsum.back() == 0.0) {
    this->next_state = -1;
    this->next_event_index = -1;
    this->next_time = std::numeric_limits<double>::max();
    return;
  }

  // (0, m_tsum.back()]
  double rand = random_generator.sample_unit_interval() * m_tsum.back();

  this->next_event_index = -1;
  for (Index i = 0; i < tsize; ++i) {
    if (rand < m_tsum[i + 1]) {
      this->next_event_index = this->event[i].event_index;
      break;
    }
  }

  this->next_state = -1;
  this->next_time_increment =
      -std::log(random_generator.sample_unit_interval()) / m_tsum.back();
  this->next_time = this->current_time + this->next_time_increment;
}

// DebugMode=false
template class EventGroup<false>;

// DebugMode=true
template class EventGroup<true>;

}  // namespace event_group
}  // namespace clexmonte
}  // namespace CASM
