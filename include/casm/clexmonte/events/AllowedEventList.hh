#ifndef CASM_clexmonte_events_AllowedEventList
#define CASM_clexmonte_events_AllowedEventList

#include <map>
#include <numeric>
#include <optional>
#include <vector>

#include "casm/clexmonte/events/EventGroup.hh"
#include "casm/clexmonte/events/ImpactTable.hh"
#include "casm/clexmonte/events/event_data.hh"

namespace CASM {

namespace clexulator {
struct ConfigDoFValues;
class PrimNeighborList;
class SuperNeighborList;
}  // namespace clexulator

namespace monte {
class OccLocation;
}  // namespace monte

namespace clexmonte {

struct AllowedEventData {
  AllowedEventData() : is_assigned(false), event_id(-1, -1) {}

  AllowedEventData(bool is_assigned, EventID event_id)
      : is_assigned(is_assigned), event_id(event_id) {}

  /// \brief Whether the event_id is assigned
  bool is_assigned;

  /// \brief The event ID
  ///
  /// This may be in an invalid state if `is_assigned` is false.
  EventID event_id;
};

/// \brief Data structure storing mapping between EventIDs and event index into
/// a vector containing the currently allowed events (plus some extra buffer
/// elements)
///
/// Example iteration:
/// ```
/// AllowedEventMap allowed_event_map;
/// allowed_event_map.begin_iteration();
/// while (allowed_event_map.next()) {
///   // do something with allowed_event_map.curr_event_id()
/// }
/// ```
///
class AllowedEventMap {
 public:
  AllowedEventMap(bool use_map_index = true)
      : m_use_map_index(use_map_index),
        m_n_assigned(0),
        m_has_been_resized(false) {}

  Index n_total() const { return m_events.size(); }

  Index n_assigned() const { return m_n_assigned; }

  void start_iteration() const {
    if (m_use_map_index) {
      m_map_it = m_map_index.begin();
    } else {
      m_vec_pos = 0;
      while (m_vec_pos < m_events.size() && !m_events[m_vec_pos].is_assigned) {
        m_vec_pos++;
      }
    }
  }

  bool next() const {
    if (m_use_map_index) {
      if (m_map_it == m_map_index.end()) {
        return false;
      }
      m_map_it++;
      return m_map_it != m_map_index.end();
    } else {
      while (m_vec_pos < m_events.size() && !m_events[m_vec_pos].is_assigned) {
        m_vec_pos++;
      }
      return m_vec_pos != m_events.size();
    }
  }

  EventID curr_event_id() const {
    if (m_use_map_index) {
      return m_map_it->first;
    } else {
      return m_events[m_vec_pos].event_id;
    }
  }

  Index curr_event_index() const {
    if (m_use_map_index) {
      return m_map_it->second;
    } else {
      return m_vec_pos;
    }
  }

  bool is_end() const {
    if (m_use_map_index) {
      return m_map_it == m_map_index.end();
    } else {
      return m_vec_pos == m_events.size();
    }
  }

  /// \brief Access the AllowedEventData from event index (undefined out of
  /// range)
  AllowedEventData &event_data(Index index) { return m_events[index]; }

  /// \brief Access the AllowedEventData from event index (undefined out of
  /// range)
  AllowedEventData const &event_data(Index index) const {
    return m_events[index];
  }

  /// \brief Get the EventID from event index (undefined out of range)
  EventID const &event_id(Index index) const {
    return m_events[index].event_id;
  }

  /// \brief Get the index of an assigned event (undefined if not assigned)
  Index event_index(EventID const &event_id) const {
    if (m_use_map_index) {
      return m_map_index.at(event_id);
    } else {
      return m_vec_index[event_id.unitcell_index][event_id.prim_event_index];
    }
  }

  /// \brief Return an iterator into `events()` if assigned, or
  /// `events().end()` if not assigned
  std::vector<AllowedEventData>::const_iterator find(
      EventID const &event_id) const {
    if (m_use_map_index) {
      auto it = m_map_index.find(event_id);
      if (it != m_map_index.end()) {
        return m_events.begin() + it->second;
      } else {
        return m_events.end();
      }
    } else {
      if (event_id.unitcell_index >= m_vec_index.size() ||
          event_id.prim_event_index >=
              m_vec_index[event_id.unitcell_index].size()) {
        return m_events.end();
      }
      Index index =
          m_vec_index[event_id.unitcell_index][event_id.prim_event_index];
      if (index == -1) {
        return m_events.end();
      } else {
        return m_events.begin() + index;
      }
    }
  }

  /// \brief Add unassigned events to `events` until `events.size()==n`
  void reserve(Index n) {
    while (m_events.size() < n) {
      m_available.push_back(m_events.size());
      m_events.push_back({false, EventID(-1, -1)});
      m_has_been_resized = true;
    }
  }

  /// \brief Assign an event ID to an event index
  ///
  /// If already assigned, return the event index; otherwise, assign the event
  /// and return the event index. If there are no unassigned events elements,
  /// then add a new element to the `events` list and set `has_been_resized` to
  /// true.
  Index assign(EventID const &event_id) {
    auto it = find(event_id);

    if (it != m_events.end()) {
      return std::distance(events().begin(), it);
    }

    if (m_available.empty()) {
      m_available.push_back(m_events.size());
      m_events.push_back({false, EventID()});
      m_has_been_resized = true;
    }

    Index index = m_available.back();
    if (m_use_map_index) {
      _set_map_index(event_id, index);
    } else {
      _set_vec_index(event_id, index);
    }
    AllowedEventData &event_data = m_events[index];
    event_data.is_assigned = true;
    event_data.event_id = event_id;
    m_n_assigned++;
    m_available.pop_back();
    return index;
  }

  /// \brief Assign an event ID to an event index
  ///
  /// - This overload uses an existing `find` result to avoid a second lookup
  Index assign(std::vector<AllowedEventData>::const_iterator it,
               EventID const &event_id) {
    if (it != m_events.end()) {
      return std::distance(events().begin(), it);
    }

    if (m_available.empty()) {
      m_available.push_back(m_events.size());
      m_events.push_back({false, EventID()});
      m_has_been_resized = true;
    }

    Index index = m_available.back();
    if (m_use_map_index) {
      _set_map_index(event_id, index);
    } else {
      _set_vec_index(event_id, index);
    }
    AllowedEventData &event_data = m_events[index];
    event_data.is_assigned = true;
    event_data.event_id = event_id;
    m_n_assigned++;
    m_available.pop_back();
    return index;
  }

  /// \brief If true, the `events` list has been expanded
  bool has_been_resized() const { return m_has_been_resized; }

  /// \brief Set the flag that the `events` list has been expanded to false
  void clear_has_been_resized() { m_has_been_resized = false; }

  /// \brief Free (un-assign) an element of `events` by index; do nothing if
  /// already unassigned; (undefined if index out of range)
  void free(Index index) {
    if (m_events[index].is_assigned) {
      EventID const &event_id = m_events[index].event_id;
      if (m_use_map_index) {
        m_map_index.erase(event_id);
      } else {
        m_vec_index[event_id.unitcell_index][event_id.prim_event_index] = -1;
      }
      m_events[index].is_assigned = false;
      m_n_assigned--;
      m_available.push_back(index);
    }
  }

  /// \brief Free (un-assign) an element of `events` by EventID; do nothing if
  /// already unassigned
  void free(EventID const &event_id) {
    if (m_use_map_index) {
      auto it = m_map_index.find(event_id);
      if (it != m_map_index.end()) {
        m_events[it->second].is_assigned = false;
        m_n_assigned--;
        m_available.push_back(it->second);
        m_map_index.erase(it);
      }
    } else {
      if (event_id.unitcell_index < m_vec_index.size() &&
          event_id.prim_event_index <
              m_vec_index[event_id.unitcell_index].size()) {
        Index &index =
            m_vec_index[event_id.unitcell_index][event_id.prim_event_index];
        if (index != -1) {
          m_events[index].is_assigned = false;
          m_n_assigned--;
          m_available.push_back(index);
          index = -1;
        }
      }
    }
  }

  /// \brief Returns a list {0, 1, 2, ..., events.size()-1}
  std::vector<Index> event_index_list() const {
    std::vector<Index> _list(events().size());
    std::iota(_list.begin(), _list.end(), 0);
    return _list;
  }

  /// \brief The current events list (contains assigned and unassigned elements)
  std::vector<AllowedEventData> const &events() const { return m_events; }

 private:
  void _set_map_index(EventID const &event_id, Index index) {
    m_map_index[event_id] = index;
  }

  void _set_vec_index(EventID const &event_id, Index index) {
    if (m_vec_index.size() < event_id.unitcell_index + 1) {
      m_vec_index.resize(event_id.unitcell_index + 1);
    }
    if (m_vec_index[event_id.unitcell_index].size() <
        event_id.prim_event_index + 1) {
      m_vec_index[event_id.unitcell_index].resize(event_id.prim_event_index + 1,
                                                  -1);
    }
    m_vec_index[event_id.unitcell_index][event_id.prim_event_index] = index;
  }

  /// \brief If true, use the map index; if false, use the vector index
  const bool m_use_map_index;

  // -- Map index data --
  // - Advantages: lower memory usage
  // - Disadvantages: more memory allocation and deallocation; slower lookup;
  //   slower iteration?
  std::map<EventID, Index> m_map_index;
  mutable std::map<EventID, Index>::const_iterator m_map_it;

  // -- Vector index data
  // - Disadvantages: higher memory usage
  // - Advantages: less memory allocation and deallocation; faster lookup;
  // faster
  //   iteration?
  std::vector<std::vector<Index>> m_vec_index;
  mutable Index m_vec_pos;

  // -- Used by both map and vector index methods --

  /// \brief A vector that specifies the currently allowed events
  ///
  /// This also contains a certain number of events that are not allowed to
  /// avoid having to resize the event rate tree too often.
  std::vector<AllowedEventData> m_events;

  /// \brief The available elements of `m_events` (where `is_assigned` is false)
  std::vector<Index> m_available;

  /// \brief The number of assigned events
  Index m_n_assigned;

  // \brief If true, the `events` list has been expanded
  bool m_has_been_resized;
};

/// \brief Data structure for KMC storing only the allowed events
///
/// This is designed to work with `RejectionFreeEventSelector` as follows:
///
/// - `AllowedEventList::allowed_event_map.events` stores the currently allowed
///   events, plus some extra elements to avoid resizing the event rate tree
///   too often.
/// - When `RejectionFreeEventSelector::select_event` is called, an event is
///   selected as an index into `AllowedEventList::allowed_event_map.events`.
/// - When `RejectionFreeEventSelector::set_impacted_events` is called, the
///   method `AllowedEventList::make_impact_list` is called which:
///
///   - Updates `AllowedEventList::selected_event_id` to the selected event
///     ID.
///   - Uses `AllowedEventList::relative_impact_table` to determine which
///     events are impacted and update `events` to contain all the impacted
///     events, if possible.
///   - If `events` does not have space to include all impacted events, then
///     `AllowedEventList::allowed_event_map.has_been_resized` is set to true;
///     otherwise, it is set to false.
///   - Updates `AllowedEventList::impact_list` to contain the elements of
///     `events` that must be updated before the next event is selected.
///   - If `AllowedEventList::allowed_event_map.has_been_resized` is true, then
///     the event selector must be rebuilt after applying the selected event,
///     but before selecting the next event.
///
/// - When `RejectionFreeEventSelector::select_event` is called and updates
///   rates via `rate_calculator_ptr->calculate_rate(event_id)`, any events
///   that are no longer allowed should be freed by calling
///   `AllowedEventList::allowed_event_map.free`.
struct AllowedEventList {
  AllowedEventList(
      std::vector<PrimEventData> const &prim_event_list,
      std::vector<EventImpactInfo> const &prim_impact_info_list,
      clexulator::ConfigDoFValues const &dof_values,
      monte::OccLocation const &occ_location,
      std::shared_ptr<clexulator::PrimNeighborList> prim_nlist,
      std::shared_ptr<clexulator::SuperNeighborList> supercell_nlist,
      bool use_map_index, bool use_neighborlist_impact_table,
      bool assign_allowed_events_only);

  /// \brief If true, use the neighborlist impact table; if false, use the
  /// relative impact table
  const bool use_neighborlist_impact_table;

  /// \brief A reference to the prim event list
  std::vector<PrimEventData> const &prim_event_list;

  /// \brief The relative impact table
  std::optional<RelativeEventImpactTable> relative_impact_table;

  /// \brief The neighbor list based impact table
  std::optional<NeighborlistEventImpactTable> neighborlist_impact_table;

  /// \brief The configuration degrees of freedom values
  clexulator::ConfigDoFValues const &dof_values;

  /// \brief The occupant location tracker
  monte::OccLocation const &occ_location;

  /// \brief The supercell neighbor list
  std::shared_ptr<clexulator::SuperNeighborList> supercell_nlist;

  /// \brief Maintains a mapping between EventIDs and event linear indices
  AllowedEventMap allowed_event_map;

  /// \brief The neighbor indices for the sites that change for each event in
  /// the prim_event_list
  ///
  /// `neighbor_index[prim_event_index][site_index]` is the index in the
  /// neighbor list for the `site_index`-th site in the `prim_event_index`-th
  /// event.
  std::vector<std::vector<int>> neighbor_index;

  /// \brief When an event is selected, this stores selected event ID
  EventID selected_event_id;

  /// \brief A list that gets updated based on the selected event to contain
  /// the elements of `events` that are impacted by the selected event
  std::vector<Index> impact_list;

  /// \brief If true (default) check if potentially impacted events are allowed
  ///     and only assign them to the event list if they are. Otherwise,
  ///     assign all potentially impacted events to the event list (whether they
  ///     are allowed will still be checked during the rate calculation).
  ///
  bool assign_allowed_events_only;

  /// \brief Returns a list of indices of events that are impacted by the
  /// selected event
  std::vector<Index> const &make_impact_list(Index selected_event_index);
};

struct GetImpactFromAllowedEventList {
  GetImpactFromAllowedEventList(
      std::shared_ptr<AllowedEventList> _allowed_event_list)
      : allowed_event_list(_allowed_event_list) {
    if (!allowed_event_list) {
      throw std::runtime_error(
          "GetImpactFromAllowedEventList: allowed_event_list is nullptr");
    }
  }

  std::vector<Index> const &operator()(Index selected_event_index) {
    return allowed_event_list->make_impact_list(selected_event_index);
  }

  std::shared_ptr<AllowedEventList> allowed_event_list;
};

}  // namespace clexmonte
}  // namespace CASM

#endif
