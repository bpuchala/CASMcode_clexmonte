#include "casm/clexmonte/system/System.hh"

#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexulator/ConfigDoFValuesTools_impl.hh"
// #include "casm/configuration/clusterography/impact_neighborhood.hh"
// #include "casm/configuration/clusterography/orbits.hh"
// #include "casm/configuration/occ_events/OccEventRep.hh"

namespace CASM {
namespace clexmonte {

namespace {

std::map<std::string, std::shared_ptr<clexulator::OrderParameter>>
make_order_parameters(
    std::map<std::string, clexulator::DoFSpace> const
        &order_parameter_definitions,
    Eigen::Matrix3l const &transformation_matrix_to_super,
    xtal::UnitCellCoordIndexConverter const &supercell_index_converter) {
  std::map<std::string, std::shared_ptr<clexulator::OrderParameter>>
      order_parameters;
  for (auto const &pair : order_parameter_definitions) {
    auto res = order_parameters.emplace(
        pair.first, std::make_shared<clexulator::OrderParameter>(pair.second));
    clexulator::OrderParameter &order_parameter = *res.first->second;
    order_parameter.update(transformation_matrix_to_super,
                           supercell_index_converter);
  }
  return order_parameters;
}

}  // namespace

/// \brief Constructor
System::System(std::shared_ptr<xtal::BasicStructure const> const &_shared_prim,
               composition::CompositionConverter const &_composition_converter)
    : prim(std::make_shared<config::Prim const>(_shared_prim)),
      composition_converter(_composition_converter),
      composition_calculator(composition_converter.components(),
                             xtal::allowed_molecule_names(*_shared_prim)),
      occevent_symgroup_rep(occ_events::make_occevent_symgroup_rep(
          prim->sym_info.unitcellcoord_symgroup_rep,
          prim->sym_info.occ_symgroup_rep,
          prim->sym_info.atom_position_symgroup_rep)) {}

/// \brief Constructor
SupercellSystemData::SupercellSystemData(
    System const &system, Eigen::Matrix3l const &transformation_matrix_to_super)
    : convert(*system.prim->basicstructure, transformation_matrix_to_super),
      occ_candidate_list(convert),
      order_parameters(
          make_order_parameters(system.order_parameter_definitions,
                                convert.transformation_matrix_to_super(),
                                convert.index_converter())) {
  // make supercell_neighbor_list
  if (system.prim_neighbor_list != nullptr) {
    supercell_neighbor_list = std::make_shared<clexulator::SuperNeighborList>(
        transformation_matrix_to_super, *system.prim_neighbor_list);
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
  for (auto const &pair : system.order_parameter_definitions) {
    auto const &key = pair.first;
    auto const &definition = pair.second;
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

namespace {

/// \brief Helper to get SupercellSystemData,
///     constructing as necessary
SupercellSystemData &get_supercell_data(
    System &system, Eigen::Matrix3l const &transformation_matrix_to_super) {
  auto it = system.supercell_data.find(transformation_matrix_to_super);
  if (it == system.supercell_data.end()) {
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
SupercellSystemData &get_supercell_data(
    System &system, monte::State<Configuration> const &state) {
  auto const &T = get_transformation_matrix_to_super(state);
  return get_supercell_data(system, T);
}

}  // namespace

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

/// \brief Helper to make the default configuration in a supercell
Configuration make_default_configuration(
    System const &system,
    Eigen::Matrix3l const &transformation_matrix_to_super) {
  return Configuration(
      transformation_matrix_to_super,
      clexulator::make_default_config_dof_values(
          system.prim->basicstructure->basis().size(),
          transformation_matrix_to_super.determinant(),
          system.prim->global_dof_info, system.prim->local_dof_info));
}

/// \brief Convert configuration from standard basis to prim basis
Configuration from_standard_values(
    System const &system,
    Configuration const &configuration_in_standard_basis) {
  Eigen::Matrix3l const &T =
      configuration_in_standard_basis.transformation_matrix_to_super;
  return Configuration(
      T, clexulator::from_standard_values(
             configuration_in_standard_basis.dof_values,
             system.prim->basicstructure->basis().size(), T.determinant(),
             system.prim->global_dof_info, system.prim->local_dof_info));
}

/// \brief Convert configuration from prim basis to standard basis
Configuration to_standard_values(
    System const &system, Configuration const &configuration_in_prim_basis) {
  Eigen::Matrix3l const &T =
      configuration_in_prim_basis.transformation_matrix_to_super;
  return Configuration(
      T, clexulator::to_standard_values(
             configuration_in_prim_basis.dof_values,
             system.prim->basicstructure->basis().size(), T.determinant(),
             system.prim->global_dof_info, system.prim->local_dof_info));
}

/// \brief Helper to make the default configuration in prim basis
monte::State<Configuration> make_default_state(
    System const &system,
    Eigen::Matrix3l const &transformation_matrix_to_super) {
  return monte::State<Configuration>(
      make_default_configuration(system, transformation_matrix_to_super));
}

/// \brief Convert configuration from standard basis to prim basis
monte::State<Configuration> from_standard_values(
    System const &system,
    monte::State<Configuration> const &state_in_standard_basis) {
  monte::State<Configuration> state_in_prim_basis{state_in_standard_basis};
  state_in_prim_basis.configuration =
      from_standard_values(system, state_in_standard_basis.configuration);
  return state_in_prim_basis;
}

/// \brief Convert configuration from prim basis to standard basis
monte::State<Configuration> to_standard_values(
    System const &system,
    monte::State<Configuration> const &state_in_prim_basis) {
  monte::State<Configuration> state_in_standard_basis{state_in_prim_basis};
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

/// \brief Helper to get the Clexulator
std::shared_ptr<clexulator::Clexulator> get_basis_set(System const &system,
                                                      std::string const &key) {
  return _verify(system.basis_sets, key, "basis_sets");
}

/// \brief Helper to get the local Clexulator
std::shared_ptr<std::vector<clexulator::Clexulator>> get_local_basis_set(
    System const &system, std::string const &key) {
  return _verify(system.local_basis_sets, key, "local_basis_sets");
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

/// \brief KMC events index definitions
std::shared_ptr<occ_events::OccSystem> get_event_system(System const &system) {
  return system.event_system;
}

/// \brief KMC events
std::map<std::string, OccEventTypeData> const &get_event_type_data(
    System const &system) {
  return system.event_type_data;
}

/// \brief KMC events
OccEventTypeData const &get_event_type_data(System const &system,
                                            std::string const &key) {
  return _verify(system.event_type_data, key, "events");
}

/// \brief Helper to get the correct clexulator::ClusterExpansion for a
///     particular state, constructing as necessary
///
/// \relates System
std::shared_ptr<clexulator::ClusterExpansion> get_clex(
    System &system, monte::State<Configuration> const &state,
    std::string const &key) {
  auto clex = _verify(get_supercell_data(system, state).clex, key, "clex");

  set(*clex, state);
  return clex;
}

/// \brief Helper to get the correct clexulator::ClusterExpansion for a
///     particular state, constructing as necessary
///
/// \relates System
std::shared_ptr<clexulator::MultiClusterExpansion> get_multiclex(
    System &system, monte::State<Configuration> const &state,
    std::string const &key) {
  auto clex =
      _verify(get_supercell_data(system, state).multiclex, key, "multiclex");
  set(*clex, state);
  return clex;
}

/// \brief Helper to get the correct clexulator::LocalClusterExpansion for a
///     particular state's supercell, constructing as necessary
std::shared_ptr<clexulator::LocalClusterExpansion> get_local_clex(
    System &system, std::string const &key,
    monte::State<Configuration> const &state) {
  auto clex =
      _verify(get_supercell_data(system, state).local_clex, key, "local_clex");
  set(*clex, state);
  return clex;
}

/// \brief Helper to get the correct clexulator::LocalClusterExpansion for a
///     particular state's supercell, constructing as necessary
std::shared_ptr<clexulator::MultiLocalClusterExpansion> get_local_multiclex(
    System &system, std::string const &key,
    monte::State<Configuration> const &state) {
  auto clex = _verify(get_supercell_data(system, state).local_multiclex, key,
                      "local_multiclex");
  set(*clex, state);
  return clex;
}

/// \brief Helper to get the correct order parameter calculators for a
///     particular configuration, constructing as necessary
///
/// \relates System
std::shared_ptr<clexulator::OrderParameter> get_order_parameter(
    System &system, monte::State<Configuration> const &state,
    std::string const &key) {
  auto order_parameter =
      _verify(get_supercell_data(system, state).order_parameters, key,
              "order_parameters");
  order_parameter->set(&get_dof_values(state));
  return order_parameter;
}

/// \brief Helper to get supercell index conversions
monte::Conversions const &get_index_conversions(
    System &system, monte::State<Configuration> const &state) {
  return get_supercell_data(system, state).convert;
}

/// \brief Helper to get unique pairs of (asymmetric unit index, species index)
monte::OccCandidateList const &get_occ_candidate_list(
    System &system, monte::State<Configuration> const &state) {
  return get_supercell_data(system, state).occ_candidate_list;
}

}  // namespace clexmonte
}  // namespace CASM