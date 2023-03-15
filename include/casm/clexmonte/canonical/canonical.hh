#ifndef CASM_clexmonte_canonical
#define CASM_clexmonte_canonical

#include <random>

#include "casm/clexmonte/definitions.hh"
#include "casm/clexulator/ClusterExpansion.hh"
#include "casm/monte/MethodLog.hh"
#include "casm/monte/RandomNumberGenerator.hh"

namespace CASM {
namespace clexmonte {
namespace canonical {

/// \brief Implements potential for canonical Monte Carlo
class CanonicalPotential {
 public:
  CanonicalPotential(std::shared_ptr<system_type> _system);

  /// \brief Reset pointer to state currently being calculated
  void set(state_type const *state, std::shared_ptr<Conditions> conditions);

  /// \brief Pointer to current state
  state_type const *state() const;

  /// \brief Pointer to current conditions
  std::shared_ptr<Conditions> const &conditions() const;

  /// \brief Calculate (extensive) cluster expansion value
  double extensive_value();

  /// \brief Calculate change in (extensive) cluster expansion value due to a
  ///     series of occupation changes
  double occ_delta_extensive_value(std::vector<Index> const &linear_site_index,
                                   std::vector<int> const &new_occ);

 private:
  /// System pointer
  std::shared_ptr<system_type> m_system;

  /// State to use
  state_type const *m_state;

  /// Conditions, depends on current state
  std::shared_ptr<Conditions> m_conditions;

  /// Formation energy cluster expansion calculator;
  std::shared_ptr<clexulator::ClusterExpansion> m_formation_energy_clex;
};

typedef CanonicalPotential potential_type;

/// \brief Helper for making a conditions ValueMap for canonical Monte
///     Carlo calculations
monte::ValueMap make_conditions(
    double temperature,
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> comp);

/// \brief Helper for making a conditions ValueMap for canonical Monte
///     Carlo calculations
monte::ValueMap make_conditions_increment(
    double temperature,
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> comp);

/// \brief Implements canonical Monte Carlo calculations
template <typename EngineType>
struct Canonical {
  explicit Canonical(std::shared_ptr<system_type> _system,
                     std::shared_ptr<EngineType> _random_number_engine =
                         std::shared_ptr<EngineType>());

  /// System data
  std::shared_ptr<system_type> system;

  /// Random number generator
  monte::RandomNumberGenerator<EngineType> random_number_generator;

  /// Update species in monte::OccLocation tracker?
  bool update_species = false;

  /// Current state
  state_type const *state;

  /// Current supercell
  Eigen::Matrix3l transformation_matrix_to_super;

  /// Occupant tracker
  monte::OccLocation const *occ_location;

  /// The current state's conditions in efficient-to-use form
  std::shared_ptr<clexmonte::Conditions> conditions;

  /// \brief Perform a single run, evolving current state
  void run(state_type &state, monte::OccLocation &occ_location,
           run_manager_type &run_manager);

  /// \brief Construct functions that may be used to sample various quantities
  /// of
  ///     the Monte Carlo calculation as it runs
  static std::map<std::string, state_sampling_function_type>
  standard_sampling_functions(
      std::shared_ptr<Canonical<EngineType>> const &calculation);

  /// \brief Construct functions that may be used to analyze Monte Carlo
  ///     calculation results
  static std::map<std::string, results_analysis_function_type>
  standard_analysis_functions(
      std::shared_ptr<Canonical<EngineType>> const &calculation);

  /// \brief Construct functions that may be used to modify states
  static std::map<std::string, state_modifying_function_type>
  standard_modifying_functions(
      std::shared_ptr<Canonical<EngineType>> const &calculation);
};

/// \brief Explicitly instantiated Canonical calculator
typedef Canonical<std::mt19937_64> Canonical_mt19937_64;

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM

#endif