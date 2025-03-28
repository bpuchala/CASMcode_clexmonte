#include "casm/clexmonte/system/System.hh"

#include "casm/clexmonte/state/Conditions.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexulator/ConfigDoFValuesTools_impl.hh"
#include "casm/monte/events/OccLocation.hh"

namespace CASM {
namespace clexmonte {

namespace {

std::map<std::string, std::shared_ptr<clexulator::OrderParameter>>
make_order_parameters(
    std::map<std::string, std::shared_ptr<clexulator::DoFSpace const>> const
        &dof_spaces,
    Eigen::Matrix3l const &transformation_matrix_to_super,
    xtal::UnitCellCoordIndexConverter const &supercell_index_converter) {
  std::map<std::string, std::shared_ptr<clexulator::OrderParameter>>
      order_parameters;
  for (auto const &pair : dof_spaces) {
    auto res = order_parameters.emplace(
        pair.first, std::make_shared<clexulator::OrderParameter>(*pair.second));
    clexulator::OrderParameter &order_parameter = *res.first->second;
    order_parameter.update(transformation_matrix_to_super,
                           supercell_index_converter);
  }
  return order_parameters;
}

}  // namespace

/// \brief Constructor
///
/// \param _shared_prim The prim
/// \param _composition_converter The composition axes
/// \param _n_dimensions Number of dimensions to use when calculating
///     properties such as kinetic coefficients. Does not actually restrict
///     calculations to a certain number of dimensions.
///
System::System(std::shared_ptr<xtal::BasicStructure const> const &_shared_prim,
               composition::CompositionConverter const &_composition_converter,
               Index _n_dimensions)
    : prim(std::make_shared<config::Prim const>(_shared_prim)),
      n_dimensions(_n_dimensions),
      composition_converter(_composition_converter),
      composition_calculator(composition_converter.components(),
                             xtal::allowed_molecule_names(*_shared_prim)),
      convert(*prim->basicstructure, Eigen::Matrix3l::Identity()),
      occevent_symgroup_rep(occ_events::make_occevent_symgroup_rep(
          prim->sym_info.unitcellcoord_symgroup_rep,
          prim->sym_info.occ_symgroup_rep,
          prim->sym_info.atom_position_symgroup_rep)),
      supercells(std::make_shared<config::SupercellSet>(prim)) {
  monte::OccCandidateList occ_candidate_list(convert);
  canonical_swaps = monte::make_canonical_swaps(convert, occ_candidate_list);
  semigrand_canonical_swaps =
      monte::make_semigrand_canonical_swaps(convert, occ_candidate_list);
}

/// \brief Constructor
SupercellSystemData::SupercellSystemData(
    System const &system, Eigen::Matrix3l const &transformation_matrix_to_super)
    : convert(*system.prim->basicstructure, transformation_matrix_to_super),
      occ_candidate_list(convert),
      order_parameters(make_order_parameters(
          system.dof_spaces, convert.transformation_matrix_to_super(),
          convert.index_converter())) {
  // make supercell_neighbor_list
  if (system.prim_neighbor_list != nullptr) {
    supercell_neighbor_list = std::make_shared<clexulator::SuperNeighborList>(
        transformation_matrix_to_super, *system.prim_neighbor_list);
  }

  // make corr
  for (auto const &pair : system.basis_sets) {
    if (supercell_neighbor_list == nullptr) {
      throw std::runtime_error(
          "Error constructing SupercellSystemData: Cannot construct corr with "
          "empty neighbor list");
    }
    auto const &key = pair.first;
    auto const &_clexulator = pair.second;
    auto _corr = std::make_shared<clexulator::Correlations>(
        supercell_neighbor_list, _clexulator);
    corr.emplace(key, _corr);
  }

  // make local_corr
  for (auto const &pair : system.local_basis_sets) {
    if (supercell_neighbor_list == nullptr) {
      throw std::runtime_error(
          "Error constructing SupercellSystemData: Cannot construct local_corr "
          "with empty neighbor list");
    }
    auto const &key = pair.first;
    auto const &_clexulator = pair.second;
    auto _corr = std::make_shared<clexulator::LocalCorrelations>(
        supercell_neighbor_list, _clexulator);
    local_corr.emplace(key, _corr);
  }

  // make clex
  for (auto const &pair : system.clex_data) {
    if (supercell_neighbor_list == nullptr) {
      throw std::runtime_error(
          "Error constructing SupercellSystemData: Cannot construct clex with "
          "empty neighbor list");
    }
    auto const &key = pair.first;
    auto const &data = pair.second;
    auto _clexulator = get_basis_set(system, data.basis_set_name);
    auto _clex = std::make_shared<clexulator::ClusterExpansion>(
        supercell_neighbor_list, _clexulator, data.coefficients);
    clex.emplace(key, _clex);
  }

  // make multiclex
  for (auto const &pair : system.multiclex_data) {
    if (supercell_neighbor_list == nullptr) {
      throw std::runtime_error(
          "Error constructing SupercellSystemData: Cannot construct multiclex "
          "with empty neighbor list");
    }
    auto const &key = pair.first;
    auto const &data = pair.second;
    auto _clexulator = get_basis_set(system, data.basis_set_name);
    auto _multiclex = std::make_shared<clexulator::MultiClusterExpansion>(
        supercell_neighbor_list, _clexulator, data.coefficients);
    multiclex.emplace(key, _multiclex);
  }

  // make local_clex
  for (auto const &pair : system.local_clex_data) {
    if (supercell_neighbor_list == nullptr) {
      throw std::runtime_error(
          "Error constructing SupercellSystemData: Cannot construct local_clex "
          "with empty neighbor list");
    }
    auto const &key = pair.first;
    auto const &data = pair.second;
    auto _local_clexulator =
        get_local_basis_set(system, data.local_basis_set_name);
    auto _local_clex = std::make_shared<clexulator::LocalClusterExpansion>(
        supercell_neighbor_list, _local_clexulator, data.coefficients);
    local_clex.emplace(key, _local_clex);
  }

  // make local_multiclex
  for (auto const &pair : system.local_multiclex_data) {
    if (supercell_neighbor_list == nullptr) {
      throw std::runtime_error(
          "Error constructing SupercellSystemData: Cannot construct "
          "local_multiclex with empty neighbor list");
    }
    auto const &key = pair.first;
    auto const &data = pair.second;
    auto _local_clexulator =
        get_local_basis_set(system, data.local_basis_set_name);
    auto _local_multiclex =
        std::make_shared<clexulator::MultiLocalClusterExpansion>(
            supercell_neighbor_list, _local_clexulator, data.coefficients);
    local_multiclex.emplace(key, _local_multiclex);
  }

  // make order_parameters
  for (auto const &pair : system.dof_spaces) {
    auto const &key = pair.first;
    auto const &definition = *pair.second;
    auto _order_parameter =
        std::make_shared<clexulator::OrderParameter>(definition);
    _order_parameter->update(convert.transformation_matrix_to_super(),
                             convert.index_converter());
    order_parameters.emplace(key, _order_parameter);
  }
}

// --- The following are used to construct a common interface between "System"
// data, in this case System, and templated CASM::clexmonte methods such as
// sampling function factory methods ---

/// \brief Helper to get SupercellSystemData,
///     constructing as necessary
SupercellSystemData &get_supercell_data(
    System &system, Eigen::Matrix3l const &transformation_matrix_to_super) {
  auto it = system.supercell_data.find(transformation_matrix_to_super);
  if (it == system.supercell_data.end()) {
    system.supercells->insert(transformation_matrix_to_super);
    it = system.supercell_data
             .emplace(
                 std::piecewise_construct,
                 std::forward_as_tuple(transformation_matrix_to_super),
                 std::forward_as_tuple(system, transformation_matrix_to_super))
             .first;
  }
  return it->second;
}

/// \brief Helper to get SupercellSystemData,
///     constructing as necessary
SupercellSystemData &get_supercell_data(System &system,
                                        state_type const &state) {
  auto const &T = get_transformation_matrix_to_super(state);
  return get_supercell_data(system, T);
}

/// \brief Helper to get std::shared_ptr<config::Prim const>
std::shared_ptr<config::Prim const> const &get_prim_info(System const &system) {
  return system.prim;
}

/// \brief Helper to get std::shared_ptr<xtal::BasicStructure const>
std::shared_ptr<xtal::BasicStructure const> const &get_prim_basicstructure(
    System const &system) {
  return system.prim->basicstructure;
}

/// \brief Helper to get prim basis
std::vector<xtal::Site> const &get_basis(System const &system) {
  return system.prim->basicstructure->basis();
}

/// \brief Helper to get basis size
Index get_basis_size(System const &system) {
  return system.prim->basicstructure->basis().size();
}

/// \brief Helper to get composition::CompositionConverter
composition::CompositionConverter const &get_composition_converter(
    System const &system) {
  return system.composition_converter;
}

/// \brief Helper to get composition::CompositionCalculator
composition::CompositionCalculator const &get_composition_calculator(
    System const &system) {
  return system.composition_calculator;
}

/// \brief Get the mol_composition from the conditions, assuming valid and
/// consistent conditions
Eigen::VectorXd get_mol_composition(System const &system,
                                    monte::ValueMap const &conditions) {
  if (conditions.vector_values.count("mol_composition")) {
    return conditions.vector_values.at("mol_composition");
  } else if (conditions.vector_values.count("param_composition")) {
    return get_composition_converter(system).mol_composition(
        conditions.vector_values.at("param_composition"));
  } else {
    throw std::runtime_error(
        "Error in get_mol_composition: "
        "conditions must have either \"mol_composition\" or "
        "\"param_composition\"");
  }
}

/// \brief Get the param_composition from the conditions, assuming valid and
/// consistent conditions
Eigen::VectorXd get_param_composition(System const &system,
                                      monte::ValueMap const &conditions) {
  if (conditions.vector_values.count("param_composition")) {
    return conditions.vector_values.at("param_composition");
  } else if (conditions.vector_values.count("mol_composition")) {
    return get_composition_converter(system).param_composition(
        conditions.vector_values.at("mol_composition"));
  } else {
    throw std::runtime_error(
        "Error in get_param_composition: "
        "conditions must have either \"mol_composition\" or "
        "\"param_composition\"");
  }
}

/// \brief Get or make a supercell
std::shared_ptr<config::Supercell const> get_supercell(
    System &system, Eigen::Matrix3l const &transformation_matrix_to_super) {
  return system.supercells->insert(transformation_matrix_to_super)
      .first->supercell;
}

/// \brief Helper to make the default configuration in a supercell
Configuration make_default_configuration(
    System const &system,
    Eigen::Matrix3l const &transformation_matrix_to_super) {
  return Configuration(system.supercells->insert(transformation_matrix_to_super)
                           .first->supercell);
}

/// \brief Convert configuration from standard basis to prim basis
Configuration from_standard_values(
    System const &system,
    Configuration const &configuration_in_standard_basis) {
  auto supercell = configuration_in_standard_basis.supercell;
  Eigen::Matrix3l const &T =
      supercell->superlattice.transformation_matrix_to_super();
  return Configuration(
      supercell,
      clexulator::from_standard_values(
          configuration_in_standard_basis.dof_values,
          system.prim->basicstructure->basis().size(), T.determinant(),
          system.prim->global_dof_info, system.prim->local_dof_info));
}

/// \brief Convert configuration from prim basis to standard basis
Configuration to_standard_values(
    System const &system, Configuration const &configuration_in_prim_basis) {
  auto supercell = configuration_in_prim_basis.supercell;
  Eigen::Matrix3l const &T =
      supercell->superlattice.transformation_matrix_to_super();
  return Configuration(
      supercell,
      clexulator::to_standard_values(
          configuration_in_prim_basis.dof_values,
          system.prim->basicstructure->basis().size(), T.determinant(),
          system.prim->global_dof_info, system.prim->local_dof_info));
}

/// \brief Helper to make the default configuration in prim basis
state_type make_default_state(
    System const &system,
    Eigen::Matrix3l const &transformation_matrix_to_super) {
  return state_type(
      make_default_configuration(system, transformation_matrix_to_super));
}

/// \brief Helper to make the Conditions object
std::shared_ptr<Conditions> make_conditions(System const &system,
                                            state_type const &state) {
  return std::make_shared<Conditions>(make_conditions_from_value_map(
      state.conditions, *get_prim_basicstructure(system),
      get_composition_converter(system), get_random_alloy_corr_f(system),
      CASM::TOL /*TODO*/));
}

/// \brief Convert configuration from standard basis to prim basis
state_type from_standard_values(System const &system,
                                state_type const &state_in_standard_basis) {
  state_type state_in_prim_basis{state_in_standard_basis};
  state_in_prim_basis.configuration =
      from_standard_values(system, state_in_standard_basis.configuration);
  return state_in_prim_basis;
}

/// \brief Convert configuration from prim basis to standard basis
state_type to_standard_values(System const &system,
                              state_type const &state_in_prim_basis) {
  state_type state_in_standard_basis{state_in_prim_basis};
  state_in_standard_basis.configuration =
      to_standard_values(system, state_in_prim_basis.configuration);
  return state_in_standard_basis;
}

template <typename MapType>
typename MapType::mapped_type &_verify(MapType &m, std::string const &key,
                                       std::string const &name) {
  auto it = m.find(key);
  if (it == m.end()) {
    std::stringstream msg;
    msg << "System error: '" << name << "' does not contain required '" << key
        << "'." << std::endl;
    throw std::runtime_error(msg.str());
  }
  return it->second;
}

template <typename MapType>
typename MapType::mapped_type const &_verify(MapType const &m,
                                             std::string const &key,
                                             std::string const &name) {
  auto it = m.find(key);
  if (it == m.end()) {
    std::stringstream msg;
    msg << "System error: '" << name << "' does not contain required '" << key
        << "'." << std::endl;
    throw std::runtime_error(msg.str());
  }
  return it->second;
}

/// \brief Check for DoFSpace
bool is_dof_space(System const &system, std::string const &key) {
  return system.dof_spaces.find(key) != system.dof_spaces.end();
}

/// \brief Helper to get the prim neighbor list for a system
std::shared_ptr<clexulator::PrimNeighborList> get_prim_neighbor_list(
    System &system) {
  return system.prim_neighbor_list;
}

/// \brief Check for basis set (Clexulator)
bool is_basis_set(System const &system, std::string const &key) {
  return system.basis_sets.find(key) != system.basis_sets.end();
}

/// \brief Check for local basis set (LocalClexulator)
bool is_local_basis_set(System const &system, std::string const &key) {
  return system.local_basis_sets.find(key) != system.local_basis_sets.end();
}

/// \brief Helper to get the Clexulator
std::shared_ptr<clexulator::Clexulator> get_basis_set(System const &system,
                                                      std::string const &key) {
  return _verify(system.basis_sets, key, "basis_sets");
}

/// \brief Helper to get BasisSetClusterInfo
std::shared_ptr<BasisSetClusterInfo const> get_basis_set_cluster_info(
    System const &system, std::string const &key) {
  return _verify(system.basis_set_cluster_info, key, "basis_set_cluster_info");
}

/// \brief Helper to get the local Clexulator
std::shared_ptr<std::vector<clexulator::Clexulator>> get_local_basis_set(
    System const &system, std::string const &key) {
  return _verify(system.local_basis_sets, key, "local_basis_sets");
}

/// \brief Helper to get LocalBasisSetClusterInfo
std::shared_ptr<LocalBasisSetClusterInfo const>
get_local_basis_set_cluster_info(System const &system, std::string const &key) {
  return _verify(system.local_basis_set_cluster_info, key,
                 "local_basis_set_cluster_info");
}

/// \brief Check for ClexData
bool is_clex_data(System const &system, std::string const &key) {
  return system.clex_data.find(key) != system.clex_data.end();
}

/// \brief Check for MultiClexData
bool is_multiclex_data(System const &system, std::string const &key) {
  return system.multiclex_data.find(key) != system.multiclex_data.end();
}

/// \brief Check for LocalClexData
bool is_local_clex_data(System const &system, std::string const &key) {
  return system.local_clex_data.find(key) != system.local_clex_data.end();
}

/// \brief Check for LocalMultiClexData
bool is_local_multiclex_data(System const &system, std::string const &key) {
  return system.local_multiclex_data.find(key) !=
         system.local_multiclex_data.end();
}

/// \brief Helper to get ClexData
///
/// \relates System
ClexData const &get_clex_data(System const &system, std::string const &key) {
  return _verify(system.clex_data, key, "clex");
}

/// \brief Helper to get MultiClexData
///
/// \relates System
MultiClexData const &get_multiclex_data(System const &system,
                                        std::string const &key) {
  return _verify(system.multiclex_data, key, "multiclex");
}

/// \brief Helper to get LocalClexData
///
/// \relates System
LocalClexData const &get_local_clex_data(System const &system,
                                         std::string const &key) {
  return _verify(system.local_clex_data, key, "local_clex");
}

/// \brief Helper to get LocalMultiClexData
///
/// \relates System
LocalMultiClexData const &get_local_multiclex_data(System const &system,
                                                   std::string const &key) {
  return _verify(system.local_multiclex_data, key, "local_multiclex");
}

/// \brief Construct impact tables
std::set<xtal::UnitCellCoord> get_required_update_neighborhood(
    System const &system, LocalClexData const &local_clex_data,
    Index equivalent_index) {
  auto const &clexulator =
      *_verify(system.local_basis_sets, local_clex_data.local_basis_set_name,
               "local_basis_sets");

  auto const &coeff = local_clex_data.coefficients;
  auto begin = coeff.index.data();
  auto end = begin + coeff.index.size();
  return clexulator[equivalent_index].site_neighborhood(begin, end);
}

/// \brief Construct impact tables
std::set<xtal::UnitCellCoord> get_required_update_neighborhood(
    System const &system, LocalMultiClexData const &local_multiclex_data,
    Index equivalent_index) {
  auto const &clexulator =
      *_verify(system.local_basis_sets,
               local_multiclex_data.local_basis_set_name, "local_basis_sets");

  std::set<xtal::UnitCellCoord> nhood;
  for (auto const &coeff : local_multiclex_data.coefficients) {
    auto begin = coeff.index.data();
    auto end = begin + coeff.index.size();
    auto tmp = clexulator[equivalent_index].site_neighborhood(begin, end);
    nhood.insert(tmp.begin(), tmp.end());
  }
  return nhood;
}

/// \brief Construct impact tables
std::set<xtal::UnitCellCoord> get_required_update_neighborhood(
    System const &system, LocalMultiClexData const &local_multiclex_data,
    Index equivalent_index, std::string const &key) {
  auto const &clexulator =
      *_verify(system.local_basis_sets,
               local_multiclex_data.local_basis_set_name, "local_basis_sets");

  if (!local_multiclex_data.coefficients_glossary.count(key)) {
    std::stringstream msg;
    msg << "Error: local_multiclex_data does not contain required "
        << "coefficients '" << key << "'." << std::endl;
    throw std::runtime_error(msg.str());
  }

  Index i_coeff = local_multiclex_data.coefficients_glossary.at(key);

  std::set<xtal::UnitCellCoord> nhood;
  auto const &coeff = local_multiclex_data.coefficients[i_coeff];
  auto begin = coeff.index.data();
  auto end = begin + coeff.index.size();
  auto tmp = clexulator[equivalent_index].site_neighborhood(begin, end);
  nhood.insert(tmp.begin(), tmp.end());
  return nhood;
}

/// \brief Single swap types for canonical Monte Carlo events
std::vector<monte::OccSwap> const &get_canonical_swaps(System const &system) {
  return system.canonical_swaps;
}

/// \brief Single swap types for semi-grand canonical Monte Carlo events
std::vector<monte::OccSwap> const &get_semigrand_canonical_swaps(
    System const &system) {
  return system.semigrand_canonical_swaps;
}

/// \brief Multiple swap types for semi-grand canonical Monte Carlo events
std::vector<monte::MultiOccSwap> const &get_semigrand_canonical_multiswaps(
    System const &system) {
  return system.semigrand_canonical_multiswaps;
}

/// \brief KMC events index definitions
std::shared_ptr<occ_events::OccSystem> get_event_system(System const &system) {
  if (system.event_system == nullptr) {
    std::stringstream msg;
    msg << "System error: event_system not provided" << std::endl;
    throw std::runtime_error(msg.str());
  }
  return system.event_system;
}

/// \brief KMC event symgroup representation
std::vector<occ_events::OccEventRep> const &get_occevent_symgroup_rep(
    System const &system) {
  return system.occevent_symgroup_rep;
}

/// \brief KMC events
std::map<std::string, OccEventTypeData> const &get_event_type_data(
    System const &system) {
  return system.event_type_data;
}

/// \brief KMC events
OccEventTypeData const &get_event_type_data(System const &system,
                                            std::string const &key) {
  return _verify(system.event_type_data, key, "kmc_events");
}

/// \brief Random alloy correlation matching
CorrCalculatorFunction get_random_alloy_corr_f(System const &system) {
  // TODO - this is a placeholder, need to implement actual function
  return [=](std::vector<Eigen::VectorXd> const &sublattice_prob) {
    throw std::runtime_error(
        "Error: random_alloy_corr_matching_pot is not yet implemented");
    return Eigen::VectorXd::Zero(1);
  };
}

// --- Supercell-specific

/// \brief Helper to get the correct clexulator::Correlations for a
///     particular state's supercell, constructing as necessary
///
/// Note:
/// - The resulting clexulator::Correlations is set to evaluate
///   all correlations for the specified state
/// - To only evaluate non-zero eci correlations, instead use
///   get_clex(...).correlations().
std::shared_ptr<clexulator::Correlations> get_corr(System &system,
                                                   state_type const &state,
                                                   std::string const &key) {
  auto corr = _verify(get_supercell_data(system, state).corr, key, "corr");
  corr->set(&get_dof_values(state));
  return corr;
}

/// \brief Helper to get the correct clexulator::Correlations for a
///     particular state's supercell, constructing as necessary
///
/// Note:
/// - The resulting clexulator::LocalCorrelations is set to
///   evaluate all correlations for the specified state
/// - To only evaluate non-zero eci correlations, instead use
///   get_local_clex(...).correlations().
std::shared_ptr<clexulator::LocalCorrelations> get_local_corr(
    System &system, state_type const &state, std::string const &key) {
  auto local_corr =
      _verify(get_supercell_data(system, state).local_corr, key, "local_corr");
  local_corr->set(&get_dof_values(state));
  return local_corr;
}

/// \brief Helper to get the correct clexulator::ClusterExpansion for a
///     particular state, constructing as necessary
///
/// \relates System
std::shared_ptr<clexulator::ClusterExpansion> get_clex(System &system,
                                                       state_type const &state,
                                                       std::string const &key) {
  auto clex = _verify(get_supercell_data(system, state).clex, key, "clex");

  set(*clex, state);
  return clex;
}

/// \brief Helper to get the correct clexulator::ClusterExpansion for a
///     particular state, constructing as necessary
///
/// Notes:
/// - The resulting object contains a clexulator::Correlations set
///   to evaluate all correlations which have non-zero eci for at
///   least one of the cluster expansions
///
/// \relates System
std::shared_ptr<clexulator::MultiClusterExpansion> get_multiclex(
    System &system, state_type const &state, std::string const &key) {
  auto clex =
      _verify(get_supercell_data(system, state).multiclex, key, "multiclex");
  set(*clex, state);
  return clex;
}

/// \brief Helper to get the correct clexulator::LocalClusterExpansion for a
///     particular state's supercell, constructing as necessary
std::shared_ptr<clexulator::LocalClusterExpansion> get_local_clex(
    System &system, state_type const &state, std::string const &key) {
  auto clex =
      _verify(get_supercell_data(system, state).local_clex, key, "local_clex");
  set(*clex, state);
  return clex;
}

/// \brief Helper to get the correct clexulator::LocalClusterExpansion for a
///     particular state's supercell, constructing as necessary
///
/// Notes:
/// - The resulting object contains a clexulator::LocalCorrelations
///   set to evaluate all correlations which have non-zero eci for at
///   least one of the local-cluster expansions
///
std::shared_ptr<clexulator::MultiLocalClusterExpansion> get_local_multiclex(
    System &system, state_type const &state, std::string const &key) {
  auto clex = _verify(get_supercell_data(system, state).local_multiclex, key,
                      "local_multiclex");
  set(*clex, state);
  return clex;
}

/// \brief Helper to get the supercell neighbor list for a
///     particular state's supercell, constructing as necessary
std::shared_ptr<clexulator::SuperNeighborList> get_supercell_neighbor_list(
    System &system, state_type const &state) {
  return get_supercell_data(system, state).supercell_neighbor_list;
}

/// \brief Helper to get the correct order parameter calculators for a
///     particular configuration, constructing as necessary
///
/// \relates System
std::shared_ptr<clexulator::OrderParameter> get_order_parameter(
    System &system, state_type const &state, std::string const &key) {
  auto order_parameter =
      _verify(get_supercell_data(system, state).order_parameters, key,
              "order_parameters");
  order_parameter->set(&get_dof_values(state));
  return order_parameter;
}

/// \brief Helper to get supercell index conversions
monte::Conversions const &get_index_conversions(System &system,
                                                state_type const &state) {
  return get_supercell_data(system, state).convert;
}

/// \brief Helper to get unique pairs of (asymmetric unit index, species index)
monte::OccCandidateList const &get_occ_candidate_list(System &system,
                                                      state_type const &state) {
  return get_supercell_data(system, state).occ_candidate_list;
}

/// \brief Make temporary monte::OccLocation if necessary
///
/// \param occ_location Reference-to-pointer. Use the pointed to OccLocation if
///     not nullptr. If nullptr, construct a temporary monte::OccLocation and
///     set `occ_location` to point at it.
/// \param tmp Where to construct temporary monte::OccLocation if
///     `occ_location` is nullptr.
/// \param calculation Where to get data needed to
///     construct temporary monte::OccLocation
/// \param update_atoms If True, construct OccLocation to track atom
///     movement. If False, do not.
/// \param save_atom_info If True, construct OccLocation to save atom initial
///     / final info. If False, do not.
void make_temporary_if_necessary(state_type const &state,
                                 monte::OccLocation *&occ_location,
                                 std::unique_ptr<monte::OccLocation> &tmp,
                                 System &system, bool update_atoms,
                                 bool save_atom_info) {
  if (!occ_location) {
    monte::Conversions const &convert = get_index_conversions(system, state);
    monte::OccCandidateList const &occ_candidate_list =
        get_occ_candidate_list(system, state);

    tmp = std::make_unique<monte::OccLocation>(convert, occ_candidate_list,
                                               update_atoms, save_atom_info);
    tmp->initialize(get_occupation(state));
    occ_location = tmp.get();
  }
}

}  // namespace clexmonte
}  // namespace CASM
