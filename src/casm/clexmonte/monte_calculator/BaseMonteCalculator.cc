#include "casm/clexmonte/monte_calculator/BaseMonteCalculator.hh"

#include "casm/clexmonte/misc/check_params.hh"

namespace CASM {
namespace clexmonte {

BaseMonteCalculator::BaseMonteCalculator(
    std::string _calculator_name, std::set<std::string> _required_basis_set,
    std::set<std::string> _required_local_basis_set,
    std::set<std::string> _required_clex,
    std::set<std::string> _required_multiclex,
    std::set<std::string> _required_local_clex,
    std::set<std::string> _required_local_multiclex,
    std::set<std::string> _required_dof_spaces,
    std::set<std::string> _required_params,
    std::set<std::string> _optional_params, bool _time_sampling_allowed,
    bool _update_atoms, bool _track_unique_atom_ids, bool _save_atom_info,
    bool _is_multistate_method)
    : calculator_name(_calculator_name),
      required_basis_set(_required_basis_set),
      required_local_basis_set(_required_local_basis_set),
      required_clex(_required_clex),
      required_multiclex(_required_multiclex),
      required_local_clex(_required_local_clex),
      required_local_multiclex(_required_local_multiclex),
      required_dof_spaces(_required_dof_spaces),
      required_params(_required_params),
      optional_params(_optional_params),
      time_sampling_allowed(_time_sampling_allowed),
      update_atoms(_update_atoms),
      track_unique_atom_ids(_track_unique_atom_ids),
      save_atom_info(_save_atom_info),
      is_multistate_method(_is_multistate_method) {
  // Use the RandomNumberGenerator default constructor to make a
  // RandomNumberEngine seeded by std::random_device.
  engine = monte::RandomNumberGenerator<engine_type>().engine;
}

BaseMonteCalculator::~BaseMonteCalculator() {}

/// \brief Standardized check for whether system has required data
void BaseMonteCalculator::_check_system() const {
  for (std::string const &key : this->required_basis_set) {
    if (!is_basis_set(*this->system, key)) {
      std::stringstream msg;
      msg << "Error preparing MonteCalculator " << this->calculator_name
          << ": no '" << key << "' basis_set.";
      throw std::runtime_error(msg.str());
    }
  }
  for (std::string const &key : this->required_local_basis_set) {
    if (!is_local_basis_set(*this->system, key)) {
      std::stringstream msg;
      msg << "Error preparing MonteCalculator " << this->calculator_name
          << ": no '" << key << "' local_basis_set.";
      throw std::runtime_error(msg.str());
    }
  }
  for (std::string const &key : this->required_clex) {
    if (!is_clex_data(*this->system, key)) {
      std::stringstream msg;
      msg << "Error preparing MonteCalculator " << this->calculator_name
          << ": no '" << key << "' clex.";
      throw std::runtime_error(msg.str());
    }
  }
  for (std::string const &key : this->required_multiclex) {
    if (!is_multiclex_data(*this->system, key)) {
      std::stringstream msg;
      msg << "Error preparing MonteCalculator " << this->calculator_name
          << ": no '" << key << "' multiclex.";
      throw std::runtime_error(msg.str());
    }
  }
  for (std::string const &key : this->required_local_clex) {
    if (!is_multiclex_data(*this->system, key)) {
      std::stringstream msg;
      msg << "Error preparing MonteCalculator " << this->calculator_name
          << ": no '" << key << "' local_clex.";
      throw std::runtime_error(msg.str());
    }
  }
  for (std::string const &key : this->required_local_multiclex) {
    if (!is_local_multiclex_data(*this->system, key)) {
      std::stringstream msg;
      msg << "Error preparing MonteCalculator " << this->calculator_name
          << ": no '" << key << "' local_multiclex.";
      throw std::runtime_error(msg.str());
    }
  }
  for (std::string const &key : this->required_dof_spaces) {
    if (!is_dof_space(*this->system, key)) {
      std::stringstream msg;
      msg << "Error preparing MonteCalculator " << this->calculator_name
          << ": no '" << key << "' dof_space.";
      throw std::runtime_error(msg.str());
    }
  }
}

/// \brief Validate calculator_params (for top-level key existence only)
///
/// Notes:
/// - Keys that start with "_" are ignored
/// - A warning is printed for keys that are not in the
///   required_params or optional_params sets
///
void BaseMonteCalculator::_check_params() const {
  check_params(params, this->required_params, this->optional_params);
}

/// \brief Clone the BaseMonteCalculator
std::unique_ptr<BaseMonteCalculator> BaseMonteCalculator::clone() const {
  return std::unique_ptr<BaseMonteCalculator>(this->_clone());
}

}  // namespace clexmonte
}  // namespace CASM
