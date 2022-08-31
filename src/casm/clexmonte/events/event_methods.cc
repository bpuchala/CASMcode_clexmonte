#include "casm/clexmonte/events/event_methods.hh"

#include "casm/monte/Conversions.hh"
#include "casm/monte/events/OccLocation.hh"

namespace CASM {
namespace clexmonte {

/// \brief Sets a monte::OccEvent consistent with the PrimEventData and
/// OccLocation
///
/// Notes:
/// - This doesn't need the current occupation state, just unchanging indices
///   into OccLocation, so monte::OccEvent can be set once per supercell and
///   does not need to be updated after an event occurs.
monte::OccEvent &set_event(monte::OccEvent &event,
                           PrimEventData const &prim_event_data,
                           xtal::UnitCell const &translation,
                           monte::OccLocation const &occ_location) {
  for (auto const &traj : prim_event_data.event) {
    for (auto const &pos : traj.position) {
      if (pos.is_in_resevoir) {
        throw std::runtime_error(
            "Error: KMC events exchanging with the resevoir is not allowed.");
      }
      if (!pos.is_atom) {
        throw std::runtime_error(
            "Error: KMC event trajectories must describe individual atoms.");
      }
    }
  }

  Index n_sites = prim_event_data.sites.size();
  Index n_atoms = prim_event_data.event.size();
  monte::Conversions const &convert = occ_location.convert();
  auto const &unitcellcoord_index_converter = convert.index_converter();

  // set e.new_occ --- specify new site occupation
  event.new_occ = prim_event_data.occ_final;

  // set e.linear_site_index --- specify sites being transformed
  event.linear_site_index.resize(n_sites);
  for (Index i = 0; i < n_sites; ++i) {
    event.linear_site_index[i] =
        unitcellcoord_index_converter(prim_event_data.sites[i] + translation);
  }

  // set e.occ_transform --- specify change in occupation variable
  event.occ_transform.resize(n_sites);
  for (Index i = 0; i < n_sites; ++i) {
    Index l = event.linear_site_index[i];
    Index asym = convert.l_to_asym(l);

    monte::OccTransform &transform = event.occ_transform[i];
    transform.mol_id = occ_location.l_to_mol_id(l);
    transform.l = l;
    transform.asym = asym;
    transform.from_species =
        convert.species_index(asym, prim_event_data.occ_init[i]);
    transform.to_species =
        convert.species_index(asym, prim_event_data.occ_final[i]);
  }

  // set e.atom_traj --- specify atom motion
  event.atom_traj.resize(n_atoms);
  for (Index i = 0; i < n_atoms; ++i) {
    occ_events::OccTrajectory const &occ_traj = prim_event_data.event[i];
    if (occ_traj.position.size() != 2) {
      throw std::runtime_error("Error: KMC event trajectories must be size 2.");
    }
    xtal::UnitCellCoord const &from_site =
        occ_traj.position[0].integral_site_coordinate + translation;
    xtal::UnitCell const &from_unitcell = from_site.unitcell();
    xtal::UnitCellCoord const &to_site =
        occ_traj.position[1].integral_site_coordinate + translation;
    xtal::UnitCell const &to_unitcell = to_site.unitcell();

    monte::AtomTraj &atom_traj = event.atom_traj[i];
    atom_traj.from.l = convert.bijk_to_l(from_site);
    atom_traj.from.mol_id = occ_location.l_to_mol_id(atom_traj.from.l);
    atom_traj.from.mol_comp = occ_traj.position[0].atom_position_index;

    atom_traj.to.l = convert.bijk_to_l(to_site);
    atom_traj.to.mol_id = occ_location.l_to_mol_id(atom_traj.to.l);
    atom_traj.to.mol_comp = occ_traj.position[1].atom_position_index;

    atom_traj.delta_ijk = to_unitcell - from_unitcell;
  }

  return event;
}

}  // namespace clexmonte
}  // namespace CASM