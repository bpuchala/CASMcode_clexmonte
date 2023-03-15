#ifndef CASM_clexmonte_parse_and_run_series
#define CASM_clexmonte_parse_and_run_series

#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/run/functions.hh"
#include "casm/clexmonte/run/io/RunParams.hh"
#include "casm/clexmonte/run/io/json/RunParams_json_io.hh"
#include "casm/clexmonte/system/io/json/System_json_io.hh"

namespace CASM {
namespace clexmonte {

template <typename CalculationType>
void parse_and_run_series(fs::path system_json_file,
                          fs::path run_params_json_file);

// --- Implementation ---

template <typename CalculationType>
void parse_and_run_series(fs::path system_json_file,
                          fs::path run_params_json_file) {
  /// Parse and construct system
  if (!fs::exists(system_json_file)) {
    std::stringstream msg;
    msg << "Error in CASM::clexmonte::parse_and_run_series: system_json_file "
           "does not exist: "
        << system_json_file;
    throw std::runtime_error(msg.str());
  }
  jsonParser system_json(system_json_file);
  InputParser<clexmonte::System> system_parser(system_json);
  std::runtime_error system_error_if_invalid{
      "Error reading Monte Carlo system JSON input"};
  report_and_throw_if_invalid(system_parser, CASM::log(),
                              system_error_if_invalid);

  std::shared_ptr<clexmonte::System> system(system_parser.value.release());

  // read run_params file
  jsonParser run_params_json(run_params_json_file);

  // parse and construct calculation options
  jsonParser calculation_options_json = jsonParser::object();
  if (run_params_json.contains("calculation_options")) {
    calculation_options_json = run_params_json["calculation_options"];
  }
  InputParser<CalculationType> calculation_parser(calculation_options_json,
                                                  system);
  std::shared_ptr<CalculationType> calculation(
      calculation_parser.value.release());

  /// Make state sampling & analysis functions
  auto sampling_functions =
      CalculationType::standard_sampling_functions(calculation);
  auto analysis_functions =
      CalculationType::standard_analysis_functions(calculation);
  auto modifying_functions =
      CalculationType::standard_modifying_functions(calculation);

  /// Make config generator / state generator / results_io JSON parsers
  auto config_generator_methods =
      clexmonte::standard_config_generator_methods(calculation->system);
  auto state_generator_methods = clexmonte::standard_state_generator_methods(
      calculation->system, modifying_functions, config_generator_methods);
  auto results_io_methods = clexmonte::standard_results_io_methods(
      sampling_functions, analysis_functions);

  /// Parse and construct run parameters
  if (!fs::exists(run_params_json_file)) {
    std::stringstream msg;
    msg << "Error in CASM::clexmonte::parse_and_run_series: "
           "run_params_json_file does not exist: "
        << run_params_json_file;
    throw std::runtime_error(msg.str());
  }
  InputParser<clexmonte::RunParams> run_params_parser(
      run_params_json, sampling_functions, analysis_functions,
      state_generator_methods, results_io_methods);
  std::runtime_error run_params_error_if_invalid{
      "Error reading Monte Carlo run parameters JSON input"};
  report_and_throw_if_invalid(run_params_parser, CASM::log(),
                              run_params_error_if_invalid);

  clexmonte::RunParams &run_params = *run_params_parser.value;

  clexmonte::run_series(*calculation, *run_params.state_generator,
                        run_params.run_manager_params,
                        run_params.sampling_fixture_params);
}

}  // namespace clexmonte
}  // namespace CASM

#endif
