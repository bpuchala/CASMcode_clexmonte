#include "casm/clexmonte/canonical/io/json/InputData_json_io.hh"

#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/canonical/io/InputData.hh"
#include "casm/clexmonte/canonical/io/json/StateGenerator_json_io.hh"
#include "casm/clexmonte/canonical/sampling_functions.hh"
#include "casm/clexmonte/system/OccSystem.hh"
#include "casm/clexmonte/system/io/json/OccSystem_json_io.hh"
#include "casm/composition/io/json/CompositionConverter_json_io.hh"
#include "casm/crystallography/io/BasicStructureIO.hh"
#include "casm/monte/results/io/json/jsonResultsIO.hh"

namespace CASM {
namespace clexmonte {
namespace canonical {

/// \brief Parse canonical Monte Carlo input file
///
/// Input file summary:
/// \code
/// {
///   "method": "canonical"
///   "kwargs": {
///     "system": {
///       "prim": <xtal::BasicStructure or file path>
///           Specifies the primitive crystal structure and allowed DoF. Must
///           be the prim used to generate the cluster expansion.
///       "composition_axes": <composition::CompositionConverter>
///           Specifies composition axes
///       "formation_energy": <Clex>
///           Input specifies a CASM cluster expansion basis set source file,
///           coefficients, and compilation settings.
///     },
///     "state_generation": <monte::StateGenerator>
///         Specifies a "path" of input states at which to run Monte Carlo
///         calculations. Each state is an initial configuration and set of
///         thermodynamic conditions (temperature, chemical potential,
///         composition, etc.).
///     "random_number_generator": <monte::RandomNumberGenerator>
///         (Future) Options controlling the random number generator.
///     "sampling": <monte::SamplingParams>
///         Options controlling which quantities are sampled and how often
///         sampling is performed.
///     "completion_check": <monte::CompletionCheck>
///         Controls when a single Monte Carlo run is complete. Options include
///         convergence of sampled quantiies, min/max number of samples, min/
///         max number of passes, etc.
///     "results_io": <monte::ResultsIO>
///         Options controlling results output.
///   }
/// }
void parse(InputParser<InputData> &parser) {
  // Parse canonical MC calculation data. Includes input:
  // - "prim"
  // - "composition_axes"
  // - "formation_energy"
  auto system_data_subparser = parser.subparse<system_type>("system");
  if (!system_data_subparser->valid()) {
    return;
  }
  std::shared_ptr<system_type> system_data =
      std::move(system_data_subparser->value);

  // Make state sampling functions, with current supercell-specific info
  monte::StateSamplingFunctionMap<config_type> sampling_functions =
      make_sampling_functions(system_data, canonical_tag());

  // Construct state generator
  auto state_generator_subparser = parser.subparse<state_generator_type>(
      "state_generation", system_data, sampling_functions, canonical_tag());

  // // Read sampling params
  // std::unique_ptr<monte::SamplingParams> sampling_params =
  //     parser.require<monte::SamplingParams>("sampling");

  // // Read completion check params
  // std::unique_ptr<monte::CompletionCheckParams> completion_check_params =
  //     parser.require<monte::CompletionCheckParams>("completion_check");

  // // Construct results I/O instance
  // auto results_io_subparser =
  //     parser.subparse<ResultsIO>("results_io");

  // Construct random number generator
  MTRand random_number_generator;

  if (!parser.valid()) {
    parser.value = std::make_unique<InputData>(
        system_data, std::move(state_generator_subparser->value),
        sampling_functions,
        monte::SamplingParams(),         // *sampling_params,
        monte::CompletionCheckParams(),  // *completion_check_params,
        std::make_unique<monte::jsonResultsIO<
            config_type>>(),  // std::move(results_io_subparser->value),
        random_number_generator);
  }
}

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM