#ifndef CASM_clexmonte_events_event_data
#define CASM_clexmonte_events_event_data

#include <string>
#include <vector>

#include "casm/configuration/clusterography/IntegralCluster.hh"
#include "casm/configuration/occ_events/OccEvent.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/global/definitions.hh"
#include "casm/misc/Comparisons.hh"
#include "casm/monte/events/OccEvent.hh"
#include "casm/monte/misc/LexicographicalCompare.hh"

namespace CASM {
namespace clexmonte {

/// \brief Data calculated for a single event in a single state
struct EventState {
  EventState()
      : is_allowed(false),
        is_normal(false),
        dE_final(0.0),
        dE_activated(0.0),
        freq(0.0),
        rate(0.0),
        reverse_rate(0.0),
        formation_energy_delta_corr(nullptr),
        local_corr(nullptr),
        Ef_kra(0.0),
        dEf_final(0.0),
        dEf_activated(0.0),
        entropy_delta_corr(nullptr),
        S_kra(0.0),
        dS_final(0.0),
        dS_activated(0.0) {}

  bool is_allowed;      ///< Is allowed given current configuration
  bool is_normal;       ///< Is "normal" (dEa > 0.0) && (dEa > dEf)
  double dE_final;      ///< Final state energy, relative to initial state
  double dE_activated;  ///< Activation energy, relative to initial state
  double freq;          ///< Attempt frequency
  double reverse_freq;  ///< Reverse attempt frequency
  double rate;          ///< Occurance rate
  double reverse_rate;  ///< Reverse rate

  // -- Evaluated for CALCMETHOD::CONST_ENTROPY and CALCMETHOD::CLEX_ENTROPY --
  Eigen::VectorXd const *
      formation_energy_delta_corr;  ///< Change in formation energy correlations
  Eigen::VectorXd const *local_corr;  ///< Local correlations
  double Ef_kra;                      ///< KRA energy
  double dEf_final;      ///< Change in formation energy to final state
  double dEf_activated;  ///< Change in formation energy to activated state

  // -- Evaluated for CALCMETHOD::CLEX_ENTROPY --
  Eigen::VectorXd const
      *entropy_delta_corr;  ///< Change in formation energy correlations
  double S_kra;             ///< KRA entropy
  double dS_final;          ///< Final state entropy, relative to initial state
  double dS_activated;      ///< Activation entropy, relative to initial state
};

// TODO Option:
//  Generalize the event state data structure,
//  Template the event state calculation
//
// struct GeneralizedEventState {
//   bool is_allowed;      ///< Is allowed given current configuration
//   bool is_normal;       ///< Is "normal" (dEa > 0.0) && (dEa > dEf)
//   double dE_final;      ///< Final state energy, relative to initial state
//   double dE_activated;  ///< Activation energy, relative to initial state
//   double freq;          ///< Attempt frequency
//   double rate;          ///< Occurance rate
//   double reverse_rate;  ///< Reverse rate
//
//   std::vector<double> scalar_params;                 ///< Scalar parameters
//   std::vector<Eigen::VectorXd> vector_params;        ///< Vector parameters
//   std::vector<Eigen::MatrixXd> matrix_params;        ///< Matrix parameters
//   std::vector<Eigen::VectorXd const *> corr_params;  ///< Correlation vectors
// };
//
// struct GeneralizedEventStateGlossary {
//   std::vector<std::string> scalar_param_names;  ///< Scalar parameter names
//
//   std::vector<std::string> vector_param_names;  ///< Vector parameter names
//
//   std::vector<std::string> matrix_param_names;  ///< Matrix parameter names
//
//   std::vector<std::string> corr_param_names;  ///< Correlation vector names
// };
//
// rate = freq * exp(-beta * dE_activated)
// struct ConstEntropyEventStateCalculator;
//
// rate = (kT/h) * exp(-beta * dE_activated)
// rate = (kT/h) * exp(-beta * (dEf_activated - T * dS_activated))
// struct ClexEntropyEventStateCalculator;

// TODO Option:
//  Template the event state data structure and calculation
//

/// \brief Data particular to a single translationally distinct event
struct EventData {
  /// \brief Unit cell linear index
  Index unitcell_index;

  /// \brief Used to apply event and track occupants in monte::OccLocation
  monte::OccEvent event;

  // Note: To keep all event state values, uncomment this:
  // /// \brief Holds the calculated event state
  // mutable EventState event_state;
};

/// \brief Data common to all translationally equivalent events
struct PrimEventData {
  /// \brief Event type name
  std::string event_type_name;

  /// \brief Equivalent event index
  Index equivalent_index;

  /// \brief Is forward trajectory (else reverse)
  bool is_forward;

  /// \brief Linear index for this prim event
  Index prim_event_index;

  /// \brief Defines event
  occ_events::OccEvent event;

  /// \brief Event sites, relative to origin unit cell
  std::vector<xtal::UnitCellCoord> sites;

  /// \brief Initial site occupation
  std::vector<int> occ_init;

  /// \brief Final site occupation
  std::vector<int> occ_final;
};

std::string name(PrimEventData const &prim_event_data);

/// \brief Data structure specifying information about the impact of a possible
///     event. Used to construct impact tables.
struct EventImpactInfo {
  /// \brief Sites whose DoF are modified when this event occurs
  std::vector<xtal::UnitCellCoord> phenomenal_sites;

  /// \brief The set of sites for which a change in DoF results in a
  ///     change in the event propensity
  std::set<xtal::UnitCellCoord> required_update_neighborhood;
};

struct RelativeEventID;

/// \brief Identifies an event via translation from the `origin` unit cell
///
/// This data structure is used to build a "relative impact table" listing which
/// events are impacted by the occurance of each possible event in the origin
/// unit cell.
struct RelativeEventID : public Comparisons<CRTPBase<RelativeEventID>> {
  /// \brief Index specifying a possible event in the `origin` unit cell
  Index prim_event_index;

  /// \brief Translation of the event from the origin unit cell
  xtal::UnitCell translation;

  /// \brief Less than comparison of RelativeEventID
  bool operator<(RelativeEventID const &rhs) const {
    if (this->translation < rhs.translation) {
      return true;
    }
    if (rhs.translation < this->translation) {
      return false;
    }
    return this->prim_event_index < rhs.prim_event_index;
  }

 private:
  friend struct Comparisons<CRTPBase<RelativeEventID>>;

  /// \brief Equality comparison of RelativeEventID
  bool eq_impl(RelativeEventID const &rhs) const {
    return this->translation == rhs.translation &&
           this->prim_event_index == rhs.prim_event_index;
  }
};

struct EventID;

/// \brief Identifies an event via linear unit cell index in some supercell
///
/// Thie unitcell index and prim event index can be used to lookup the correct
/// local clexulator and neighbor list information for evaluating local
/// correlations and updating global correlations.
struct EventID : public Comparisons<CRTPBase<EventID>> {
  EventID() = default;

  EventID(Index _prim_event_index, Index _unitcell_index)
      : prim_event_index(_prim_event_index), unitcell_index(_unitcell_index) {}

  /// \brief Index specifying a possible event in the `origin` unit cell
  Index prim_event_index;

  /// \brief Linear unit cell index into a supercell, as determined by
  /// xtal::UnitCellIndexConverter
  Index unitcell_index;

  /// \brief Less than comparison of EventID
  bool operator<(EventID const &rhs) const {
    if (this->unitcell_index < rhs.unitcell_index) {
      return true;
    }
    if (this->unitcell_index > rhs.unitcell_index) {
      return false;
    }
    return this->prim_event_index < rhs.prim_event_index;
  }

 private:
  friend struct Comparisons<CRTPBase<EventID>>;

  /// \brief Equality comparison of Configuration
  bool eq_impl(EventID const &rhs) const {
    return this->unitcell_index == rhs.unitcell_index &&
           this->prim_event_index == rhs.prim_event_index;
  }
};

struct EventFilterGroup {
  /// The linear unit cell index for which the group applies
  std::set<Index> unitcell_index;

  /// Whether events are included (default) or excluded as a default
  bool include_by_default = true;

  /// The prim event index of excluded/included events
  std::set<Index> prim_event_index;
};

struct SelectedEvent {
  EventID event_id;
  Index event_index;
  double total_rate;
  double time;
  double time_increment;
  PrimEventData const *prim_event_data;
  EventData const *event_data;

  /// \brief Event state data, re-calculated for the selected event if
  /// `requires_event_state` is true for one of the selected event data
  /// functions; otherwise, nullptr
  EventState const *event_state;

  /// \brief Event group index, if applicable
  Index group;

  /// \brief Event group state index, if applicable
  Index group_state;

  void reset() {
    event_id = EventID();
    event_index = -1;
    total_rate = 0.0;
    time_increment = 0.0;
    prim_event_data = nullptr;
    event_data = nullptr;
    event_state = nullptr;
    group = -1;
    group_state = -1;
  }
};

/// \brief Get selected event data needed for the collecting functions
///
/// Notes:
/// - prim_event_list should be present after calculation->reset()
struct SelectedEventInfo {
  std::vector<PrimEventData> const &prim_event_list;
  std::shared_ptr<std::vector<Index>> prim_event_index_to_index;
  std::shared_ptr<std::vector<bool>> prim_event_index_to_has_value;

  std::vector<std::string> partition_names;
  std::map<Eigen::VectorXl, std::string, monte::LexicographicalCompare>
      value_labels;

  SelectedEventInfo(std::vector<PrimEventData> const &_prim_event_list);

  /// \brief Construct `prim_event_index_to_index` so that indices differentiate
  ///     by event type
  void make_indices_by_type();

  /// \brief Construct `prim_event_index_to_index` so that indices differentiate
  ///     by event type and equivalent index
  void make_indices_by_equivalent_index();

  /// \brief Construct `prim_event_index_to_index` so that indices differentiate
  ///     by event type, equivalent index, and direction
  void make_indices_by_equivalent_index_and_direction();

  /// \brief Construct `prim_event_index_to_index` so that indices differentiate
  ///     by event type and equivalent index, but only for a single event type
  void make_indices_by_equivalent_index_per_event_type(
      std::string event_type_name);
};

}  // namespace clexmonte
}  // namespace CASM

#endif
