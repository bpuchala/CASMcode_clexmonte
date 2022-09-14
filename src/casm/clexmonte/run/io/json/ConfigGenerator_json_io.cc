#include "casm/clexmonte/run/io/json/ConfigGenerator_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/misc/polymorphic_method_json_io.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/clexulator/io/json/ConfigDoFValues_json_io.hh"
#include "casm/monte/state/FixedConfigGenerator.hh"

namespace CASM {
namespace clexmonte {

/// \brief Construct ConfigGenerator from JSON
///
/// A configuration generation method generates a configuration given a set of
/// conditions and results from previous runs. It may be a way to customize a
/// state generation method.
///
/// Expected:
///   method: string (required)
///     The name of the chosen config generation method. Currently, the only
///     option is:
///     - "fixed": monte::FixedConfigGenerator
///
///   kwargs: dict (optional, default={})
///     Method-specific options. See documentation for particular methods:
///     - "fixed": `parse(InputParser<monte::FixedConfigGenerator> &, ...)`
void parse(
    InputParser<config_generator_type> &parser,
    MethodParserMap<config_generator_type> const &config_generator_methods) {
  parse_polymorphic_method(parser, config_generator_methods);
}

/// \brief Construct FixedConfigGenerator from JSON
///
/// Requires:
/// - `Configuration from_standard_values(
///        system_type const &system,
///        Configuration const &configuration)`
/// - `Configuration make_default_configuration(
///        system_type const &system,
///        Eigen::Matrix3l const &transformation_matrix_to_super)`
void parse(InputParser<monte::FixedConfigGenerator<Configuration>> &parser,
           std::shared_ptr<system_type> const &system) {
  Eigen::Matrix3l T;
  parser.require(T, "transformation_matrix_to_super");
  if (!parser.valid()) {
    return;
  }

  // TODO: validation of dof types and dimensions?
  // Note: expect "dof" to be in standard basis
  std::unique_ptr<clexulator::ConfigDoFValues> standard_dof_values =
      parser.optional<clexulator::ConfigDoFValues>("dof");

  if (parser.valid()) {
    if (standard_dof_values != nullptr) {
      parser.value =
          notstd::make_unique<monte::FixedConfigGenerator<Configuration>>(
              from_standard_values(*system,
                                   Configuration(T, *standard_dof_values)));
    } else {
      parser.value =
          notstd::make_unique<monte::FixedConfigGenerator<Configuration>>(
              make_default_configuration(*system, T));
    }
  }
}

}  // namespace clexmonte
}  // namespace CASM
