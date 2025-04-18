#include "casm/clexmonte/events/AllowedEventList.hh"

#include "casm/clexmonte/events/event_methods.hh"
#include "casm/clexulator/ConfigDoFValues.hh"
#include "casm/clexulator/NeighborList.hh"
#include "casm/monte/Conversions.hh"
#include "casm/monte/events/OccLocation.hh"

// debug
#include "casm/casm_io/Log.hh"

namespace CASM {
namespace clexmonte {

AllowedEventList::AllowedEventList(
    std::vector<PrimEventData> const &_prim_event_list,
    std::vector<EventImpactInfo> const &prim_impact_info_list,
    clexulator::ConfigDoFValues const &_dof_values,
    monte::OccLocation const &_occ_location,
    std::shared_ptr<clexulator::PrimNeighborList> prim_nlist,
    std::shared_ptr<clexulator::SuperNeighborList> _supercell_nlist,
    bool use_map_index, bool _use_neighborlist_impact_table,
    bool _assign_allowed_events_only)
    : use_neighborlist_impact_table(_use_neighborlist_impact_table),
      prim_event_list(_prim_event_list),
      dof_values(_dof_values),
      occ_location(_occ_location),
      supercell_nlist(_supercell_nlist),
      allowed_event_map(use_map_index),
      assign_allowed_events_only(_assign_allowed_events_only) {
  if (prim_event_list.size() != prim_impact_info_list.size()) {
    throw std::runtime_error(
        "Error in AllowedEventList constructor: prim_event_list and "
        "prim_impact_info_list size mismatch");
  }
  if (prim_nlist == nullptr) {
    throw std::runtime_error(
        "Error in AllowedEventList constructor: prim_nlist is nullptr");
  }
  if (supercell_nlist == nullptr) {
    throw std::runtime_error(
        "Error in AllowedEventList constructor: supercell_nlist is nullptr");
  }
  if (use_neighborlist_impact_table) {
    neighborlist_impact_table = NeighborlistEventImpactTable(
        prim_impact_info_list,
        _occ_location.convert().unitcell_index_converter(), prim_nlist,
        _occ_location.convert().transformation_matrix_to_super(),
        _supercell_nlist);
  } else {
    relative_impact_table = RelativeEventImpactTable(
        prim_impact_info_list,
        _occ_location.convert().unitcell_index_converter());
  }

  auto const &unitcell_index_converter =
      occ_location.convert().unitcell_index_converter();
  Index n_unitcells = unitcell_index_converter.total_sites();

  // Construct `neighbor_index`
  this->neighbor_index.reserve(prim_event_list.size());
  for (Index prim_event_index = 0; prim_event_index < prim_event_list.size();
       ++prim_event_index) {
    PrimEventData const &prim_event_data = prim_event_list[prim_event_index];

    std::vector<int> _neighbor_index;
    for (auto const &site : prim_event_data.sites) {
      _neighbor_index.push_back(prim_nlist->neighbor_index(site));
    }
    this->neighbor_index.push_back(_neighbor_index);
  }

  // Assign allowed events to `allowed_event_map`
  Index group = 0;
  Index index_in_group = -1;
  Index max_n_impacted = 0;
  std::vector<Index> linear_site_index;
  for (Index unitcell_index = 0; unitcell_index < n_unitcells;
       ++unitcell_index) {
    for (Index prim_event_index = 0; prim_event_index < prim_event_list.size();
         ++prim_event_index) {
      PrimEventData const &prim_event_data = prim_event_list[prim_event_index];

      // set linear_site_index
      set_event_linear_site_index(linear_site_index, unitcell_index,
                                  this->neighbor_index[prim_event_index],
                                  *supercell_nlist);

      // if event is allowed, assign
      if (event_is_allowed(linear_site_index, dof_values, prim_event_data)) {
        // set event_id
        EventID event_id(prim_event_index, unitcell_index);

        // assign event
        allowed_event_map.assign(event_id, group, index_in_group);
      }
    }
  }

  // Add elements to `events` to avoid resizing the event rate tree too often
  Index target_n_events = static_cast<Index>(
      std::ceil((allowed_event_map.events().size() + max_n_impacted) * 1.2));
  allowed_event_map.reserve(target_n_events);
}

/// \brief Returns a list of indices of events that are impacted by the
/// selected event (undefined if `selected_event_index` is out of range)
///
/// Also sets `allowed_event_map.has_been_resized()` to false if `events` does
/// not change size or true if it does
std::vector<Index> const &AllowedEventList::make_impact_list(
    Index selected_event_index) {
  // get `selected_event_id`
  EventID const &selected_event_id =
      this->allowed_event_map.event_id(selected_event_index);

  // set `impact_list` and update `events`
  impact_list.clear();
  this->allowed_event_map.clear_has_been_resized();
  std::vector<EventID> const &impacted_event_ids =
      use_neighborlist_impact_table
          ? this->neighborlist_impact_table.value()(selected_event_id)
          : this->relative_impact_table.value()(selected_event_id);

  if (assign_allowed_events_only) {
    // approach 1: include assigned events,
    // and only assign new events that are allowed
    Index group = 0;
    Index index_in_group = -1;

    static std::vector<Index> linear_site_index;
    for (auto const &event_id : impacted_event_ids) {
      // check if already assigned
      auto it = this->allowed_event_map.find(event_id);
      if (it != this->allowed_event_map.events().end()) {
        // if already assigned, add to impact list
        this->impact_list.push_back(
            std::distance(this->allowed_event_map.events().begin(), it));
      }
      // if not assigned, check if event is allowed and assign it if allowed
      else {
        // set linear_site_index
        PrimEventData const &prim_event_data =
            prim_event_list[event_id.prim_event_index];
        set_event_linear_site_index(
            linear_site_index, event_id.unitcell_index,
            this->neighbor_index[event_id.prim_event_index], *supercell_nlist);
        if (event_is_allowed(linear_site_index, dof_values, prim_event_data)) {
          // assign event
          this->impact_list.push_back(
              allowed_event_map.assign(it, event_id, group, index_in_group));
        }
      }
    }
  } else {
    Index group = 0;
    Index index_in_group = -1;

    // approach 2: include all possible impacted events
    for (auto const &event_id : impacted_event_ids) {
      this->impact_list.push_back(
          this->allowed_event_map.assign(event_id, group, index_in_group));
    }
  }

  return this->impact_list;
}

}  // namespace clexmonte
}  // namespace CASM
