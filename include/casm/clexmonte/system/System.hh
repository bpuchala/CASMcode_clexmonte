#ifndef CASM_clexmonte_system_System
#define CASM_clexmonte_system_System

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/events/event_data.hh"
#include "casm/clexmonte/misc/Matrix3lCompare.hh"
#include "casm/clexmonte/state/LocalOrbitCompositionCalculator.hh"
#include "casm/clexmonte/system/system_data.hh"
#include "casm/clexulator/ClusterExpansion.hh"
#include "casm/clexulator/DoFSpace.hh"
#include "casm/clexulator/LocalClusterExpansion.hh"
#include "casm/clexulator/NeighborList.hh"
#include "casm/clexulator/OrderParameter.hh"
#include "casm/composition/CompositionCalculator.hh"
#include "casm/composition/CompositionConverter.hh"
#include "casm/configuration/Prim.hh"
#include "casm/configuration/SupercellSet.hh"
#include "casm/configuration/clusterography/IntegralCluster.hh"
#include "casm/configuration/occ_events/OccEvent.hh"
#include "casm/configuration/occ_events/OccEventRep.hh"
#include "casm/configuration/occ_events/OccSystem.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/monte/Conversions.hh"
#include "casm/monte/ValueMap.hh"
#include "casm/monte/events/OccCandidate.hh"

namespace CASM {

namespace monte {
template <typename _ConfigType>
struct State;
}

namespace clexmonte {

struct System;
struct SupercellSystemData;

/// \brief Data structure for holding Monte Carlo calculation data and methods
///     that should only exist once, and should be accessible by
///     sampling functions - occupation DoF
///
/// Notes:
/// - Use the standalone `get_supercell_data` helper methods to get
///   supercell-specific canonical Monte Carlo calculation data, constructing
///   it as necessary
/// - Use the standalone `get_clex` helper method to get
///   supercell-specific clexulator::ClusterExpansion instance for a given
///   state, constructing it as necessary
struct System {
  /// \brief Constructor
  System(std::shared_ptr<xtal::BasicStructure const> const &_shared_prim,
         composition::CompositionConverter const &_composition_converter,
         Index _n_dimensions = 3);

  // --- Crystal structure

  /// Primitive crystal structure and allowed degrees of freedom (DoF)
  std::shared_ptr<config::Prim const> prim;

  /// Number of dimensions (used for example in normalizing kinetic
  /// coefficients)
  Index n_dimensions;

  // --- Composition

  /// Composition axes and parametric composition conversions functor
  composition::CompositionConverter composition_converter;

  /// Composition calculation functor
  composition::CompositionCalculator composition_calculator;

  /// Performs index conversions in supercell
  monte::Conversions convert;

  // --- Order parameters

  /// DoFSpace that define order parameters
  std::map<std::string, std::shared_ptr<clexulator::DoFSpace const>> dof_spaces;

  /// Indices of DoFSpace basis vectors forming subspaces in which
  /// order parameter magnitudes are to be calculated
  std::map<std::string, std::vector<std::vector<Index>>> dof_subspaces;

  // --- Cluster expansions

  /// Prim neighbor list
  ///
  /// Notes:
  /// - Make sure to use this->prim_neighbor_list to construct
  ///   Clexulators stored in this->basis_sets and this->local_basis_sets
  ///   so that SupercellSystemData can properly construct
  ///   the SuperNeighborList needed to evaluate correlations
  std::shared_ptr<clexulator::PrimNeighborList> prim_neighbor_list;

  /// Cluster expansion basis sets
  ///
  /// Notes:
  /// - Maps basis set name -> Clexulator
  /// - Make sure to use this->prim_neighbor_list to construct the
  ///   Clexulator so that SupercellSystemData can propertly construct
  ///   the SuperNeighborList needed to evaluate correlations
  std::map<std::string, std::shared_ptr<clexulator::Clexulator>> basis_sets;

  /// Cluster expansion basis sets cluster info
  ///
  /// Notes:
  /// - Maps basis set name -> BasisSetClusterInfo
  std::map<std::string, std::shared_ptr<BasisSetClusterInfo const>>
      basis_set_cluster_info;

  /// Data used to construct clexulator::ClusterExpansion. Contains:
  /// - basis_set_name
  /// - clexulator::SparseCoefficients
  std::map<std::string, ClexData> clex_data;

  /// Data used to construct clexulator::MultiClusterExpansion. Contains:
  /// - basis_set_name
  /// - std::vector<clexulator::SparseCoefficients>
  std::map<std::string, MultiClexData> multiclex_data;

  // --- Local cluster expansions

  /// Local cluster expansion basis sets
  ///
  /// Notes:
  /// - Maps local basis set name -> std::vector<clexulator::Clexulator>
  /// - Make sure to use this->prim_neighbor_list to construct the
  ///   Clexulator so that SupercellSystemData can propertly construct
  ///   the SuperNeighborList needed to evaluate correlations
  std::map<std::string, std::shared_ptr<std::vector<clexulator::Clexulator>>>
      local_basis_sets;

  /// Local cluster expansion basis set equivalence info
  ///
  /// Maps local basis set name -> EquivalentsInfo
  std::map<std::string, EquivalentsInfo> equivalents_info;

  /// Local cluster expansion basis sets cluster info
  ///
  /// Notes:
  /// - Maps basis set name -> LocalBasisSetClusterInfo
  std::map<std::string, std::shared_ptr<LocalBasisSetClusterInfo const>>
      local_basis_set_cluster_info;

  /// Parameters for LocalOrbitCompositionCalculator
  std::map<std::string,
           std::shared_ptr<LocalOrbitCompositionCalculatorData const>>
      local_orbit_composition_calculator_data;

  /// Data used to construct clexulator::LocalClusterExpansion. Contains:
  /// - local_basis_set_name
  /// - clexulator::SparseCoefficients
  std::map<std::string, LocalClexData> local_clex_data;

  /// Data used to construct clexulator::LocalMultiClusterExpansion. Contains:
  /// - local_basis_set_name
  /// - std::vector<clexulator::SparseCoefficients>
  std::map<std::string, LocalMultiClexData> local_multiclex_data;

  // --- Monte Carlo events

  /// Single swap types for canonical Monte Carlo events
  ///
  /// Defaults to all possible canonical swaps
  std::vector<monte::OccSwap> canonical_swaps;

  /// Single swap types for semi-grand canonical Monte Carlo events
  ///
  /// Defaults to all possible semi-grand canonical swaps
  std::vector<monte::OccSwap> semigrand_canonical_swaps;

  /// Multiple swap types for semi-grand canonical Monte Carlo events
  ///
  /// Defaults to empty
  std::vector<monte::MultiOccSwap> semigrand_canonical_multiswaps;

  // --- KMC events

  /// KMC events index definitions
  std::shared_ptr<occ_events::OccSystem> event_system;

  /// KMC event symgroup representation
  std::vector<occ_events::OccEventRep> occevent_symgroup_rep;

  /// KMC events
  std::map<std::string, OccEventTypeData> event_type_data;

  /// Linear list of events associated with the origin unit cell, including
  /// each type of event, each equivalent event, and the forward and reverse
  /// events (if distinct)
  std::vector<PrimEventData> prim_event_list;

  // --- Supercells

  /// Supercells
  std::shared_ptr<config::SupercellSet> supercells;

  /// Supercell specific formation energy calculation data and methods (using
  /// transformation_matrix_to_super as key).
  std::map<Eigen::Matrix3l, SupercellSystemData, Matrix3lCompare>
      supercell_data;

  // -- Additional parameters --

  jsonParser additional_params;
};

/// \brief Data structure for holding supercell-specific Monte Carlo calculation
///     data and methods that should only exist once, and should be accessible
///     by sampling functions - occupation DoF
struct SupercellSystemData {
  /// \brief Constructor
  SupercellSystemData(System const &system,
                      Eigen::Matrix3l const &transformation_matrix_to_super);

  /// Number of unit cells in the supercell
  Index n_unitcells;

  // --- Index conversions and occupation tracking

  /// Performs index conversions in supercell
  monte::Conversions convert;

  /// List of unique pairs of (asymmetric unit index, species index)
  monte::OccCandidateList occ_candidate_list;

  // --- Cluster expansion

  /// SuperNeighborList, used for evaluating correlations in a particular
  /// supercell
  std::shared_ptr<clexulator::SuperNeighborList> supercell_neighbor_list;

  // --- Order parameter

  /// Order parameter calculators
  std::map<std::string, std::shared_ptr<clexulator::OrderParameter>>
      order_parameters;

  // --- Cluster expansion

  /// CASM::monte correlation calculators - calculate all correlations
  std::map<std::string, std::shared_ptr<clexulator::Correlations>> corr;

  /// CASM::monte local correlation calculators - calculate all correlations
  std::map<std::string, std::shared_ptr<clexulator::LocalCorrelations>>
      local_corr;

  /// CASM::monte compatible cluster expansion calculators. Contains:
  /// -  clexulator::Correlations - calculate non-zero eci correlations
  /// -  clexulator::SparseCoefficients
  std::map<std::string, std::shared_ptr<clexulator::ClusterExpansion>> clex;

  /// CASM::monte compatible cluster expansion calculators. Contains:
  /// -  clexulator::Correlations
  /// -  clexulator::SparseCoefficients
  std::map<std::string, std::shared_ptr<clexulator::MultiClusterExpansion>>
      multiclex;

  /// CASM::monte compatible local cluster expansion calculators. Contains:
  /// -  clexulator::LocalCorrelations
  /// -  clexulator::SparseCoefficients
  std::map<std::string, std::shared_ptr<clexulator::LocalClusterExpansion>>
      local_clex;

  /// CASM::monte compatible local cluster expansion calculators. Contains:
  /// -  clexulator::LocalCorrelations
  /// -  clexulator::SparseCoefficients
  std::map<std::string, std::shared_ptr<clexulator::MultiLocalClusterExpansion>>
      local_multiclex;
};

// ---
// The following are used to construct a common interface between "System"
// data, in this case System, and templated CASM::clexmonte methods such as
// sampling function factory methods
// ---

/// \brief Helper to get SupercellSystemData,
///     constructing as necessary
SupercellSystemData &get_supercell_data(
    System &system, Eigen::Matrix3l const &transformation_matrix_to_super);

/// \brief Helper to get SupercellSystemData,
///     constructing as necessary
SupercellSystemData &get_supercell_data(System &system,
                                        state_type const &state);

/// \brief Helper to get std::shared_ptr<config::Prim const>
std::shared_ptr<config::Prim const> const &get_prim_info(System const &system);

/// \brief Helper to get std::shared_ptr<xtal::BasicStructure const>
std::shared_ptr<xtal::BasicStructure const> const &get_prim_basicstructure(
    System const &system);

/// \brief Helper to get prim basis
std::vector<xtal::Site> const &get_basis(System const &system);

/// \brief Helper to get basis size
Index get_basis_size(System const &system);

/// \brief Helper to get composition::CompositionConverter
composition::CompositionConverter const &get_composition_converter(
    System const &system);

/// \brief Helper to get composition::CompositionCalculator
composition::CompositionCalculator const &get_composition_calculator(
    System const &system);

/// \brief Get the mol_composition from the conditions, assuming valid and
/// consistent conditions
Eigen::VectorXd get_mol_composition(System const &system,
                                    monte::ValueMap const &conditions);

/// \brief Get the param_composition from the conditions, assuming valid and
/// consistent conditions
Eigen::VectorXd get_param_composition(System const &system,
                                      monte::ValueMap const &conditions);

/// \brief Get or make a supercell
std::shared_ptr<config::Supercell const> get_supercell(
    System &system, Eigen::Matrix3l const &transformation_matrix_to_super);

/// \brief Helper to make the default configuration in prim basis
Configuration make_default_configuration(
    System const &system,
    Eigen::Matrix3l const &transformation_matrix_to_super);

/// \brief Convert configuration from standard basis to prim basis
Configuration from_standard_values(
    System const &system, Configuration const &configuration_in_standard_basis);

/// \brief Convert configuration from prim basis to standard basis
Configuration to_standard_values(
    System const &system, Configuration const &configuration_in_prim_basis);

/// \brief Helper to make the default configuration in prim basis
state_type make_default_state(
    System const &system,
    Eigen::Matrix3l const &transformation_matrix_to_super);

/// \brief Helper to make the Conditions object
std::shared_ptr<Conditions> make_conditions(System const &system,
                                            state_type const &state);

/// \brief Convert configuration from standard basis to prim basis
Configuration from_standard_values(
    System const &system, Configuration const &configuration_in_standard_basis);

/// \brief Convert configuration from prim basis to standard basis
Configuration to_standard_values(
    System const &system, Configuration const &configuration_in_prim_basis);

/// \brief Check for DoFSpace
bool is_dof_space(System const &system, std::string const &key);

/// \brief Helper to get the prim neighbor list for a system
std::shared_ptr<clexulator::PrimNeighborList> get_prim_neighbor_list(
    System &system);

/// \brief Check for basis set (Clexulator)
bool is_basis_set(System const &system, std::string const &key);

/// \brief Check for local basis set (LocalClexulator)
bool is_local_basis_set(System const &system, std::string const &key);

/// \brief Helper to get the Clexulator
std::shared_ptr<clexulator::Clexulator> get_basis_set(System const &system,
                                                      std::string const &key);

/// \brief Helper to get BasisSetClusterInfo
std::shared_ptr<BasisSetClusterInfo const> get_basis_set_cluster_info(
    System const &system, std::string const &key);

/// \brief Helper to get the local Clexulator
std::shared_ptr<std::vector<clexulator::Clexulator>> get_local_basis_set(
    System const &system, std::string const &key);

/// \brief Helper to get LocalBasisSetClusterInfo
std::shared_ptr<LocalBasisSetClusterInfo const>
get_local_basis_set_cluster_info(System const &system, std::string const &key);

/// \brief Check for ClexData
bool is_clex_data(System const &system, std::string const &key);

/// \brief Check for MultiClexData
bool is_multiclex_data(System const &system, std::string const &key);

/// \brief Check for LocalClexData
bool is_local_clex_data(System const &system, std::string const &key);

/// \brief Check for LocalMultiClexData
bool is_local_multiclex_data(System const &system, std::string const &key);

/// \brief Helper to get ClexData
ClexData const &get_clex_data(System const &system, std::string const &key);

/// \brief Helper to get MultiClexData
MultiClexData const &get_multiclex_data(System const &system,
                                        std::string const &key);

/// \brief Helper to get LocalClexData
LocalClexData const &get_local_clex_data(System const &system,
                                         std::string const &key);

/// \brief Helper to get LocalMultiClexData
LocalMultiClexData const &get_local_multiclex_data(System const &system,
                                                   std::string const &key);

/// \brief Construct impact tables
std::set<xtal::UnitCellCoord> get_required_update_neighborhood(
    System const &system, LocalClexData const &local_clex_data,
    Index equivalent_index);

/// \brief Construct impact tables
std::set<xtal::UnitCellCoord> get_required_update_neighborhood(
    System const &system, LocalMultiClexData const &local_multiclex_data,
    Index equivalent_index);

/// \brief Construct impact tables
std::set<xtal::UnitCellCoord> get_required_update_neighborhood(
    System const &system, LocalMultiClexData const &local_multiclex_data,
    Index equivalent_index, std::string const &key);

/// \brief Single swap types for canonical Monte Carlo events
std::vector<monte::OccSwap> const &get_canonical_swaps(System const &system);

/// \brief Single swap types for semi-grand canonical Monte Carlo events
std::vector<monte::OccSwap> const &get_semigrand_canonical_swaps(
    System const &system);

/// \brief Multiple swap types for semi-grand canonical Monte Carlo events
std::vector<monte::MultiOccSwap> const &get_semigrand_canonical_multiswaps(
    System const &system);

/// \brief KMC events index definitions
std::shared_ptr<occ_events::OccSystem> get_event_system(System const &system);

/// \brief KMC event symgroup representation
std::vector<occ_events::OccEventRep> const &get_occevent_symgroup_rep(
    System const &system);

/// \brief KMC events
std::map<std::string, OccEventTypeData> const &get_event_type_data(
    System const &system);

/// \brief KMC events
OccEventTypeData const &get_event_type_data(System const &system,
                                            std::string const &key);

/// \brief Random alloy correlation matching
CorrCalculatorFunction get_random_alloy_corr_f(System const &system);

// --- Supercell-specific

/// \brief Helper to get the correct clexulator::Correlations for a
///     particular state's supercell, constructing as necessary
std::shared_ptr<clexulator::Correlations> get_corr(System &system,
                                                   state_type const &state,
                                                   std::string const &key);

/// \brief Helper to get the correct clexulator::LocalCorrelations for a
///     particular state's supercell, constructing as necessary
std::shared_ptr<clexulator::LocalCorrelations> get_local_corr(
    System &system, state_type const &state, std::string const &key);

/// \brief Helper to get the correct clexulator::ClusterExpansion for a
///     particular state's supercell, constructing as necessary
std::shared_ptr<clexulator::ClusterExpansion> get_clex(System &system,
                                                       state_type const &state,
                                                       std::string const &key);

/// \brief Helper to get the correct clexulator::MultiClusterExpansion for a
///     particular state's supercell, constructing as necessary
std::shared_ptr<clexulator::MultiClusterExpansion> get_multiclex(
    System &system, state_type const &state, std::string const &key);

/// \brief Helper to get the correct clexulator::LocalClusterExpansion for a
///     particular state's supercell, constructing as necessary
std::shared_ptr<clexulator::LocalClusterExpansion> get_local_clex(
    System &system, state_type const &state, std::string const &key);

/// \brief Helper to get the correct clexulator::MultiLocalClusterExpansion for
/// a
///     particular state's supercell, constructing as necessary
std::shared_ptr<clexulator::MultiLocalClusterExpansion> get_local_multiclex(
    System &system, state_type const &state, std::string const &key);

/// \brief Helper to get the supercell neighbor list for a
///     particular state's supercell, constructing as necessary
std::shared_ptr<clexulator::SuperNeighborList> get_supercell_neighbor_list(
    System &system, state_type const &state);

/// \brief Helper to get the correct order parameter calculators for a
///     particular state's supercell, constructing as necessary
std::shared_ptr<clexulator::OrderParameter> get_order_parameter(
    System &system, state_type const &state, std::string const &key);

/// \brief Helper to get supercell index conversions
monte::Conversions const &get_index_conversions(System &system,
                                                state_type const &state);

/// \brief Helper to get unique pairs of (asymmetric unit index, species index)
monte::OccCandidateList const &get_occ_candidate_list(System &system,
                                                      state_type const &state);

/// \brief Make temporary monte::OccLocation if necessary
void make_temporary_if_necessary(state_type const &state,
                                 monte::OccLocation *&occ_location,
                                 std::unique_ptr<monte::OccLocation> &tmp,
                                 System &system, bool update_atoms,
                                 bool save_atom_info);

}  // namespace clexmonte
}  // namespace CASM

#endif
