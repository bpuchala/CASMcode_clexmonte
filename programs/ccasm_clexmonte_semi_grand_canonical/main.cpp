#include "casm/casm_io/Log.hh"
#include "casm/clexmonte/run/io/json/parse_and_run_series.hh"
#include "casm/clexmonte/semi_grand_canonical.hh"

using namespace CASM;

// ///////////////////////////////////////
// ccasm-clexmonte-semi-grand-canonical main:

void print_help() {
  int verbosity = Log::standard;
  bool show_clock = false;
  int indent_space = 4;
  Log log(std::cout, verbosity, show_clock, indent_space);
  log.set_width(80);

  log.paragraph(
      "usage: ccasm-clexmonte-semi-grand-canonical [-h] [-V] system.json "
      "run_params.json");
  log << std::endl;

  log.paragraph(
      "ccasm-clexmonte-semi-grand-canonical is a program for running "
      "semi-grand canonical Monte Carlo calculations using cluster expansions "
      "generated by CASM.");
  log << std::endl;

  log << "Options:" << std::endl << std::endl;

  // ## positional arguments
  log << "positional arguments:" << std::endl;
  log.increase_indent();

  // ### systemfile
  log.indent() << "system.json" << std::endl;
  log.increase_indent();
  log.paragraph("JSON formatted file specifying the Monte Carlo system");
  log.decrease_indent();
  log << std::endl;

  // ### run_params_json_file
  log.indent() << "run_params.json" << std::endl;
  log.increase_indent();
  log.paragraph("JSON formatted file specifying Monte Carlo run parameters");
  log.decrease_indent();
  log << std::endl;

  log.decrease_indent();
  log << std::endl;

  // ## optional arguments
  log << "optional arguments:" << std::endl;
  log.increase_indent();

  // ### -h
  log.indent() << "-h, --help" << std::endl;
  log.increase_indent();
  log.paragraph("Print help message and exit");
  log.decrease_indent();
  log << std::endl;

  // ### -V
  log.indent() << "-V, --version" << std::endl;
  log.increase_indent();
  log.paragraph("Print version number and exit");
  log.decrease_indent();
  log << std::endl;
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    print_help();
    return 1;
  }

  std::string param = argv[1];

  if (param == "-h" || param == "--help") {
    print_help();
    return 0;
  } else if (param == "-V" || param == "--version") {
    log() << "2.0.0-alpha" << std::endl;
    return 0;
  } else if (argc != 3) {
    print_help();
    return 1;
  }

  fs::path system_json_file = argv[1];
  fs::path run_params_json_file = argv[2];

  if (!fs::exists(system_json_file)) {
    log() << "Error: file does not exist: " << system_json_file << std::endl;
    return 1;
  }
  if (!fs::exists(run_params_json_file)) {
    log() << "Error: file does not exist: " << run_params_json_file
          << std::endl;
    return 1;
  }

  using namespace CASM::clexmonte;
  using namespace CASM::clexmonte::semi_grand_canonical;
  try {
    parse_and_run_series<SemiGrandCanonical_mt19937_64>(system_json_file,
                                                        run_params_json_file);
  } catch (std::exception &e) {
    log() << e.what() << std::endl;
  }

  return 0;
}
