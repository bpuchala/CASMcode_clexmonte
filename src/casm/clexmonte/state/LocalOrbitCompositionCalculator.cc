#include "casm/clexmonte/state/LocalOrbitCompositionCalculator.hh"

#include "casm/clexulator/ConfigDoFValues.hh"

// #include "casm/clexmonte/definitions.hh"
// #include "casm/clexmonte/state/Configuration.hh"
// #include "casm/clexmonte/system/System.hh"
// #include "casm/configuration/Configuration.hh"
// #include "casm/monte/run_management/State.hh"

namespace CASM {
namespace clexmonte {

// LocalOrbitCompositionCalculator::LocalOrbitCompositionCalculator(
//     std::shared_ptr<system_type> _system, std::string _event_type_name,
//     std::set<int> _orbits_to_calculate)
//     : m_system(_system),
//       m_event_type_name(_event_type_name),
//       m_orbits_to_calculate(_orbits_to_calculate) {
//   // Make m_occ_index_to_component_index_converter
//   auto const &composition_calculator = get_composition_calculator(*m_system);
//   m_occ_index_to_component_index_converter =
//       composition::make_occ_index_to_component_index_converter(
//           composition_calculator.components(),
//           composition_calculator.allowed_occs());
//
//   // Setup m_num_each_component_by_orbit and validate orbits_to_calculate
//   auto cluster_info =
//       get_local_basis_set_cluster_info(*m_system, m_event_type_name);
//   int n_orbits = 0;
//   if (cluster_info->orbits.size() > 0) {
//     n_orbits = cluster_info->orbits[0].size();
//   }
//
//   for (int orbit_index : m_orbits_to_calculate) {
//     if (orbit_index < 0 || orbit_index >= n_orbits) {
//       std::stringstream msg;
//       msg << "Error in LocalOrbitCompositionCalculator: "
//           << "orbit_to_calculate=" << orbit_index << " out of range [0,"
//           << n_orbits << ").";
//       throw std::runtime_error(msg.str());
//     }
//   }
//
//   m_num_each_component_by_orbit.resize(
//       composition_calculator.components().size(),
//       m_orbits_to_calculate.size());
//   m_num_each_component_by_orbit.setZero();
// }
//
///// \brief Reset pointer to state currently being calculated
// void LocalOrbitCompositionCalculator::set(state_type const *state) {
//   // supercell-specific
//   m_state = state;
//   if (m_state == nullptr) {
//     throw std::runtime_error(
//         "Error setting LocalOrbitCompositionCalculator state: state is "
//         "empty");
//   }
//
//   // set shell composition calculation data
//   auto cluster_info =
//       get_local_basis_set_cluster_info(*m_system, m_event_type_name);
//
//   // Make m_local_orbits_neighbor_indices:
//   m_local_orbits_neighbor_indices.clear();
//   m_supercell_nlist = get_supercell_neighbor_list(*m_system, *m_state);
//   auto const &convert = get_index_conversions(*m_system, *m_state);
//   auto const &supercell_index_converter = convert.index_converter();
//   for (Index equivalent_index = 0;
//        equivalent_index < cluster_info->orbits.size(); ++equivalent_index) {
//     std::vector<std::set<std::pair<int, int>>> _neighbor_indices_by_orbit;
//     for (auto const &orbit : cluster_info->orbits[equivalent_index]) {
//       std::set<std::pair<int, int>> _neighbor_indices;
//       for (auto const &cluster : orbit) {
//         for (auto const &site : cluster) {
//           Index site_index = supercell_index_converter(site);
//           _neighbor_indices.emplace(
//               m_supercell_nlist->neighbor_index(site_index),
//               site.sublattice());
//         }
//       }
//       _neighbor_indices_by_orbit.emplace_back(std::move(_neighbor_indices));
//     }
//     m_local_orbits_neighbor_indices.emplace_back(
//         std::move(_neighbor_indices_by_orbit));
//   }
// }

/// \brief Constructor - for a single supercell
///
/// \param _orbits The cluster orbits, in order matching a Clexulator, by
///     equivalent index:
///     - The cluster `orbits[equivalent_index][orbit_index][j]` is `j`-th
///       cluster equivalent to the prototype cluster
///       `orbits[equivalent_index][orbit_index][0]` around the
///       `equivalent_index`-th equivalent phenomenal cluster, in the
///       `orbit_index`-th orbit.
/// \param _orbits_to_calculate Orbits to calculate
/// \param _combine_orbits If true, calculate the number of each component for
///      the union of the orbits in `_orbits_to_calculate`. If false, calculate
///      the number of each component for each orbit in `_orbits_to_calculate`
///      individually. If true, the resulting value will be a matrix with a
///      single column, if false, the value will be a matrix with a column for
///      each orbit.
/// \param supercell_nlist Supercell neighbor list
/// \param supercell_index_converter Converter from linear site index to
///     unitcell index and sublattice index
/// \param _dof_values Pointer to the configuration to be calculated (optional).
///     If not provided, the configuration must be set before calling
///     `calculate_num_each_component`.
LocalOrbitCompositionCalculator::LocalOrbitCompositionCalculator(
    std::vector<std::vector<std::set<clust::IntegralCluster>>> const &_orbits,
    std::set<int> _orbits_to_calculate, bool _combine_orbits,
    std::shared_ptr<clexulator::SuperNeighborList> _supercell_nlist,
    xtal::UnitCellCoordIndexConverter const &_supercell_index_converter,
    composition::CompositionCalculator const &_composition_calculator,
    clexulator::ConfigDoFValues const *_dof_values)
    : m_orbits_to_calculate(_orbits_to_calculate),
      m_combine_orbits(_combine_orbits),
      m_supercell_nlist(_supercell_nlist) {
  // Make m_occ_index_to_component_index_converter
  m_occ_index_to_component_index_converter =
      composition::make_occ_index_to_component_index_converter(
          _composition_calculator.components(),
          _composition_calculator.allowed_occs());

  // Validate orbits_to_calculate:
  int n_orbits = 0;
  if (_orbits.size() > 0) {
    n_orbits = _orbits[0].size();
  }

  for (int orbit_index : m_orbits_to_calculate) {
    if (orbit_index < 0 || orbit_index >= n_orbits) {
      std::stringstream msg;
      msg << "Error in LocalOrbitCompositionCalculator: "
          << "orbit_to_calculate=" << orbit_index << " out of range [0,"
          << n_orbits << ").";
      throw std::runtime_error(msg.str());
    }
  }

  // Setup m_num_each_component_by_orbit
  m_num_each_component_by_orbit.resize(
      _composition_calculator.components().size(),
      m_orbits_to_calculate.size());
  m_num_each_component_by_orbit.setZero();

  // Make m_local_orbits_neighbor_indices:
  m_local_orbits_neighbor_indices.clear();
  for (Index equivalent_index = 0; equivalent_index < _orbits.size();
       ++equivalent_index) {
    std::vector<std::set<std::pair<int, int>>> _neighbor_indices_by_orbit;
    for (auto const &orbit : _orbits[equivalent_index]) {
      std::set<std::pair<int, int>> _neighbor_indices;
      for (auto const &cluster : orbit) {
        for (auto const &site : cluster) {
          Index site_index = _supercell_index_converter(site);
          _neighbor_indices.emplace(
              m_supercell_nlist->neighbor_index(site_index), site.sublattice());
        }
      }
      _neighbor_indices_by_orbit.emplace_back(std::move(_neighbor_indices));
    }
    m_local_orbits_neighbor_indices.emplace_back(
        std::move(_neighbor_indices_by_orbit));
  }

  // (Optional) Set configuration to be calculated
  if (_dof_values) {
    set(_dof_values);
  }
}

/// \brief Reset pointer to configuration currently being calculated
void LocalOrbitCompositionCalculator::set(
    clexulator::ConfigDoFValues const *dof_values) {
  m_dof_values = dof_values;
  if (m_dof_values == nullptr) {
    throw std::runtime_error(
        "Error setting LocalOrbitCompositionCalculator dof_values: "
        "dof_values is empty");
  }
}

/// \brief Value at particular unit cell and phenomenal cluster
Eigen::MatrixXi const &LocalOrbitCompositionCalculator::value(
    Index unitcell_index, Index equivalent_index) {
  Eigen::VectorXi const &occupation = m_dof_values->occupation;

  std::vector<Index> const &neighbor_index_to_linear_site_index =
      m_supercell_nlist->sites(unitcell_index);

  // indices[orbit_index] = std::set<std::pair<int, int>>
  std::vector<std::set<std::pair<int, int>>> const &indices =
      m_local_orbits_neighbor_indices[equivalent_index];

  m_num_each_component_by_orbit.setZero();
  int col = 0;
  for (int orbit_index : m_orbits_to_calculate) {
    for (auto const &pair : indices[orbit_index]) {
      int neighbor_index = pair.first;
      int sublattice_index = pair.second;
      int site_index = neighbor_index_to_linear_site_index[neighbor_index];
      int occ_index = occupation(site_index);
      int component_index =
          m_occ_index_to_component_index_converter[sublattice_index][occ_index];
      m_num_each_component_by_orbit.col(col)(component_index) += 1;
    }
    ++col;
  }

  return m_num_each_component_by_orbit;
}

}  // namespace clexmonte
}  // namespace CASM
