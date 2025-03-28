#include "casm/clexmonte/system/io/json/System_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/events/event_methods.hh"
#include "casm/clexmonte/misc/parse_array.hh"
#include "casm/clexmonte/misc/subparse_from_file.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/clexmonte/system/io/json/system_data_json_io.hh"
#include "casm/clexulator/NeighborList.hh"
#include "casm/clexulator/io/json/Clexulator_json_io.hh"
#include "casm/clexulator/io/json/DoFSpace_json_io.hh"
#include "casm/clexulator/io/json/SparseCoefficients_json_io.hh"
#include "casm/composition/io/json/CompositionConverter_json_io.hh"
#include "casm/configuration/clusterography/io/json/IntegralCluster_json_io.hh"
#include "casm/configuration/clusterography/orbits.hh"
#include "casm/configuration/occ_events/io/json/OccEvent_json_io.hh"
#include "casm/configuration/occ_events/io/json/OccSystem_json_io.hh"
#include "casm/crystallography/io/BasicStructureIO.hh"
#include "casm/monte/events/io/OccCandidate_json_io.hh"

namespace CASM {
namespace clexmonte {

namespace {

/// \brief Parse
template <typename ParserType, typename RequiredType>
bool parse_from_files_object(ParserType &parser, fs::path option,
                             std::vector<fs::path> search_path,
                             std::vector<RequiredType> &vec,
                             std::map<std::string, Index> &glossary) {
  auto obj_it = parser.self.find(option);
  if (obj_it == parser.self.end() || !obj_it->is_obj()) {
    parser.insert_error(option, "Missing required JSON object");
    return false;
  }
  Index i = 0;
  for (auto it = obj_it->begin(); it != obj_it->end(); ++it) {
    auto subparser = subparse_from_file<RequiredType>(
        parser, option / std::to_string(i), search_path);
    if (!subparser->valid()) {
      return false;
    }
    vec.push_back(std::move(*subparser->value));
    glossary.emplace(it.name(), i);
    ++i;
  }
  return true;
}

/// \brief Parse from file
template <typename ParserType, typename RequiredType, typename... Args>
bool parse_from_file(ParserType &parser, fs::path option,
                     std::vector<fs::path> search_path, RequiredType &value,
                     Args &&...args) {
  auto subparser = subparse_from_file<RequiredType>(
      parser, option, search_path, std::forward<Args>(args)...);
  if (!subparser->valid()) {
    return false;
  }
  value = std::move(*subparser->value);
  return true;
}

/// \brief Parse and validate "basis_set" or "local_basis_set" (name)
template <typename ParserType, typename BasisSetMapType>
bool parse_and_validate_basis_set_name(ParserType &parser, fs::path option,
                                       std::string &basis_set_name,
                                       BasisSetMapType const &basis_sets) {
  parser.require(basis_set_name, option);

  // validate basis_set_name exists
  auto basis_set_it = basis_sets.find(basis_set_name);
  if (basis_set_it == basis_sets.end()) {
    parser.insert_error(option, "No basis set with matching name");
    return false;
  }
  return true;
}

/// \brief Parse "kmc_events"/<event_name> from JSON
///
/// TODO: document format (see
/// tests/unit/clexmonte/data/kmc/system_template.json)
/// - event : occ_events::OccEvent file path
/// - local_basis_set: name (matching one in System::local_basis_sets /
/// equivalents_info)
/// - coefficients/kra: SparseCoefficients file path
/// - coefficients/freq: SparseCoefficients file path
template <typename ParserType>
bool parse_event(
    ParserType &parser, fs::path option, std::vector<fs::path> search_path,
    std::map<std::string, OccEventTypeData> &event_type_data,
    std::map<std::string, LocalMultiClexData> &local_multiclex,
    xtal::BasicStructure const &prim,
    std::map<std::string,
             std::shared_ptr<std::vector<clexulator::Clexulator>>> const
        &local_basis_sets,
    std::map<std::string, std::shared_ptr<LocalBasisSetClusterInfo const>>
        &local_basis_set_cluster_info,
    std::map<std::string, EquivalentsInfo> const &equivalents_info,
    std::vector<occ_events::OccEventRep> const &occevent_symgroup_rep,
    occ_events::OccSystem const &event_system) {
  auto &log = CASM::log();
  log.begin_section<Log::verbose>();

  std::string event_name = option.filename();

  log.indent() << "- Getting \"coefficients\"..." << std::endl;
  auto coeffs_it = parser.self.find_at(option / "coefficients");
  if (coeffs_it == parser.self.end()) {
    parser.insert_error(option / "coefficients",
                        "Missing coefficients info for this local basis set");
    log.indent() << "- Error: no \"coefficients\" input" << std::endl;
    log.indent() << "- Not adding " << event_name << std::endl;
    log.end_section();
    return false;
  }
  if (!coeffs_it->is_obj()) {
    parser.insert_error(
        option / "coefficients",
        "Coefficients info must be a JSON object of <key>:<file path>");
    log.indent() << "- Error: \"coefficients\" input must be a JSON object of "
                    "<key>:<file path> pairs"
                 << std::endl;
    log.indent() << "- Not adding " << event_name << std::endl;
    log.end_section();
    return false;
  }
  Index coeffs_size = coeffs_it->size();
  log.indent() << "- Found " << coeffs_size << " coefficients for this event"
               << std::endl;

  // --- Create a local multi-cluster expansion for the event type ---
  LocalMultiClexData curr_local_multiclex;
  curr_local_multiclex.coefficients.resize(coeffs_size);

  // parse option / "coefficients" / <name> : <coefficients file path>
  Index i = 0;
  for (auto it = coeffs_it->begin(); it != coeffs_it->end(); ++it) {
    log.indent() << "- Getting \"" << it.name() << "\" coefficients..."
                 << std::endl;
    std::string key = it.name();
    if (!parse_from_file(parser, option / "coefficients" / key, search_path,
                         curr_local_multiclex.coefficients[i])) {
      log.indent() << "- Getting \"" << it.name() << "\" coefficients: failed"
                   << std::endl;
      log.indent() << "- Not adding " << event_name << std::endl;
      log.end_section();
      return false;
    } else {
      log.indent() << "- Getting \"" << it.name()
                   << "\" coefficients: succeeded" << std::endl;
    }
    curr_local_multiclex.coefficients_glossary.emplace(key, i);
    ++i;
  }

  // parse and validate "local_basis_set"
  log.indent() << "- Getting \"local_basis_set\"..." << std::endl;
  if (!parse_and_validate_basis_set_name(
          parser, option / "local_basis_set",
          curr_local_multiclex.local_basis_set_name, local_basis_sets)) {
    log.indent() << "- Invalid local basis set name." << std::endl;
    log.indent() << "- Not adding " << event_name << std::endl;
    log.end_section();
    return false;
  } else {
    log.indent() << "- Local basis set: "
                 << curr_local_multiclex.local_basis_set_name << std::endl;
  }
  if (equivalents_info.find(curr_local_multiclex.local_basis_set_name) ==
      equivalents_info.end()) {
    parser.insert_error(option / "local_basis_set",
                        "Missing equivalents_info for this local basis set");
    log.indent() << "- Missing equivalents_info for this local basis set."
                 << std::endl;
    log.indent() << "- Not adding " << event_name << std::endl;
    log.end_section();
    return false;
  }
  // get local basis set cluster info if available
  {
    auto it = local_basis_set_cluster_info.find(
        curr_local_multiclex.local_basis_set_name);
    if (it != local_basis_set_cluster_info.end()) {
      curr_local_multiclex.cluster_info = it->second;
      log.indent() << "- Added cluster_info" << std::endl;
    } else {
      log.indent() << "- No cluster_info found" << std::endl;
    }
  }

  // Save the LocalMultiClexData
  log.indent() << "- Adding local_multiclex data for " << event_name
               << std::endl;
  local_multiclex.emplace(event_name, curr_local_multiclex);

  /// --- Parse "event", construct and validate equivalents ---

  // parse "event" (from file)
  log.indent() << "- Getting \"event\"..." << std::endl;
  auto event_subparser = subparse_from_file<occ_events::OccEvent>(
      parser, option / "event", search_path, event_system);
  if (!event_subparser->valid()) {
    log.indent() << "- Getting \"event\": failed" << std::endl;
    log.indent() << "- Not adding " << event_name << std::endl;
    log.end_section();
    return false;
  } else {
    log.indent() << "- Getting \"event\": succeeded" << std::endl;
  }

  // construct event_data using default constructor
  OccEventTypeData curr_event_type_data;
  curr_event_type_data.prototype_event = *event_subparser->value;
  curr_event_type_data.local_multiclex_name = event_name;

  // generate equivalent events
  log.indent() << "- Generating equivalent events..." << std::endl;
  EquivalentsInfo const &info =
      equivalents_info.at(curr_local_multiclex.local_basis_set_name);
  curr_event_type_data.events = make_equivalents(
      curr_event_type_data.prototype_event, info, occevent_symgroup_rep);
  log.indent() << "- Generated " << curr_event_type_data.events.size()
               << " equivalent events" << std::endl;

  // double-check consistency of event phenomenal clusters and
  // the local basis set phenomenal clusters
  if (!is_same_phenomenal_clusters(curr_event_type_data.events, info)) {
    parser.insert_error(
        option / "local_basis_set",
        "Error generating equivalent events. There is a mismatch between the "
        "event and the phenomenal clusters of the local basis set.");

    log.indent() << "- There is a mismatch between the event and the "
                    "phenomenal clusters of the local basis set."
                 << std::endl;
    log.indent() << "- Not adding " << event_name << std::endl;
    log.end_section();
    return false;
  }

  event_type_data.emplace(event_name, curr_event_type_data);
  log.indent() << "- Added event " << event_name << std::endl;
  log.end_section();
  return true;
}

}  // namespace

/// \brief Parse System from JSON
///
/// Expected format:
/// \code
///   "prim": <xtal::BasicStructure or file path>
///       Specifies the primitive crystal structure and allowed DoF. Must
///       be the prim used to generate the cluster expansion.
///
///   "n_dimensions": int = 3
///       Number of dimensions to use when calculating properties such as
///       kinetic coefficients. Does not actually restrict calculations
///       to a certain number of dimensions.
///
///   "composition_axes": <composition::CompositionConverter>
///       Specifies composition axes
///
///   "basis_sets": object (optional)
///       Input specifies one or more CASM cluster expansion basis sets. A JSON
///       object containing one or more of:
///
///           "<name>": {
///             "source": "<path>",
///             "basis": "<path>"
///           }
///
///        where "<name>" is the basis set name, and "source" is the path to a
///        CASM clexulator source file, i.e.
///
///            "/path/to/basis_sets/bset.energy/ZrO_Clexulator_energy.cc"
///
///        and "basis" is the path to a CASM "basis.json" file for the basis
///        set, i.e.
///
///            "/path/to/basis_sets/bset.energy/basis.json"
///
///        The "basis" input is required for the "formation_energy" and
///        "entropy" basis sets when performing kinetic Monte Carlo
///        simulations.
///
///
///   "local_basis_sets": object (optional)
///       Input specifies one or more CASM local-cluster expansion basis sets.
///.      A JSON object containing one or more of:
///
///           "<name>": {
///             "source": "<path>",
///             "equivalents_info": "<path">,
///             "basis": "<path>",
///             "local_orbit_composition": {
///               "<key>": {
///                 "event": "<event name>",
///                 "orbits_to_calculate": [orbit indices... ],
///                 "combine_orbits": bool,
///                 "max_size": int
///               },
///               ...
///             }
///           }
///
///        where "<name>" is the local-cluster basis set name, and "source" is
///        the path to a CASM clexulator source file, i.e.
///
///            "/path/to/basis_sets/bset.local/ZrO_Clexulator_local.cc",
///
///        and "equivalents_info" is the path to a "equivalents_info.json" file,
///        i.e.
///
///            "/path/to/basis_sets/bset.local/equivalents_info.json"
///
///        and "basis" is the path to a "basis.json" file,
///        i.e.
///
///            "/path/to/basis_sets/bset.local/basis.json"
///
///
///        The "local_orbit_composition" input is optional and specifies how to
///        calculate the local orbit composition for each event. Both
///        "equivalents_info" and "basis" are required for
///        "local_orbit_composition". The key is used to specify which local
///        orbit composition calculator to collect data from during a Monte
///        Carlo simulation. The value is a JSON object containing:
///
///            "event": string
///                The name of the event type to calculate the local
///                composition for.
///            "orbits_to_calculate": array of int
///                The indices of the local orbit to calculate the composition
///                of.
///            "combine_orbits": bool
///                If true, calculate the number of each component for the sites
///                in the union of the orbits_to_calculate; else, calculate the
///                number of each component for the sites in each of the
///                orbits_to_calculate individually. If true, the output is a
///                matrix with a single column, where each row corresponds to
///                a component in the order defined by the system's composition
///                calculator. If false, the output is a matrix with a column
///                for each orbit in orbits_to_calculate.
///            "max_size": int
///                The maximum number of distinct local compositions to track.
///                Once this number is reached, the count of new compositions
///                is stored in the "out-of-range" bin.
///
///         For each KMC event, several default local orbit composition
///         calculator functions are constructed: one for each point cluster
///         orbit, one for all orbits calculated individually, and one for all
///         orbits calculated with sites combined.
///
///   "clex": object (optional)
///       Input specifies one or more cluster expansions. A JSON object
///       containing one or more of:
///
///           "<name>": {
///             "basis_set": "<basis set name>",
///             "coefficients": "<path>"
///           }
///
///        where "<name>" is the cluster expansion name, "<basis set name>" is
///        the name of the basis set for the cluster expansion (matching a key
///        in "basis_sets"), and "coefficients" is the path to a CASM
///        cluster expansion coefficients file.
///
///   "multiclex": object (optional)
///       Input specifies one or more cluster expansions which share basis sets.
///       Similar to "clex", but "coefficients" is a JSON object containing
///       one or more coefficients files:
///
///           "<name>": {
///             "basis_set": "<basis set name>",
///             "coefficients": {
///               "<key>": "<path>"
///               "<key>": "<path>",
///               ...
///           }
///
///   "local_clex": object (optional)
///       Input specifies one or more local-cluster expansions. A JSON object
///       containing one or more of:
///
///           "<name>": {
///             "local_basis_set": "<local basis set name>",
///             "coefficients": "<path>"
///           }
///
///        where "<name>" is the cluster expansion name, "<local basis set
///        name>" is the name of the local-cluster basis set for the
///        local-cluster expansion (matching a key in "local_basis_sets"), and
///        "coefficients" is the path to a CASM cluster expansion coefficients
///        file.
///
///   "local_multiclex": object (optional)
///       Input specifies one or more local-cluster expansions which share
///       basis sets. Similar to "local_clex", but "coefficients" is a JSON
///       object containing one or more coefficients files:
///
///           "<name>": {
///             "local_basis_set": "<local basis set name>",
///             "coefficients": {
///               "<key>": "<path>"
///               "<key>": "<path>",
///               ...
///           }
///
///   "kmc_events": object (optional)
///       Input specifies KMC events. A JSON object specifiying one or
///       more events:
///
///           "<event name>": {
///             "event": "<event description file path>",
///             "local_basis_set": "<local basis set name>",
///             "coefficients": {
///               "kra": "<kra path>",
///               "freq": "<freq path>"
///             }
///           }
///
///       where "<event name>" is a name given to the event,
///       "<event description file path>" is the path to a file defining a KMC
///       event (CASM::occ_events::OccEvent in JSON format), "<local basis set
///       name>" is the name of the local-cluster basis set for the
///       local-cluster expansion (matching a key in "local_basis_sets"), "<kra
///       path>" is the path to a cluster expansion coefficients file for the
///       event KRA, and
///       "<freq path>" is the path to a cluster expansion coefficients file
///       for the event attempt frequency. A local multi-clex with name matching
///       "<event name>" will be generated.
///
///   "event_system": string (optional)
///       Input specifies the path to a file specifying how some KMC events are
///       defined (a CASM::occ_events::OccSystem in JSON format).
///
///   "dof_spaces": object (optional)
///        Key:value pairs, where the keys are DoFSpace names and the value is
///        a DoFSpace or path to a file containing a DoFSpace.
///
///   "dof_subspaces": object (optional)
///        Key:value pairs, where the keys are DoFSpace names, and the values
///        arrays of arrays of int. The inner-most arrays are indices of
///        DoFSpace basis vectors forming subspaces in which order parameter
///        magnitudes are to be calculated.
/// \endcode
///
void parse(InputParser<System> &parser, std::vector<fs::path> search_path,
           bool verbose) {
  auto &log = CASM::log();
  int initial_verbosity = log.verbosity();
  if (verbose) {
    log.set_verbosity(Log::verbose);
  }

  log.begin<Log::verbose>("Reading System from JSON");

  // Parse "prim"
  log.indent() << "Parsing required \"prim\"..." << std::endl;
  std::shared_ptr<xtal::BasicStructure const> shared_prim =
      parser.require<xtal::BasicStructure>("prim", TOL);

  // Parse "composition_axes"
  log.indent() << "Parsing required \"composition_axes\"..." << std::endl;
  std::unique_ptr<composition::CompositionConverter> composition_axes =
      parser.require<composition::CompositionConverter>("composition_axes");

  if (!parser.valid()) {
    log.indent() << "Errors encountered while parsing System from JSON."
                 << std::endl
                 << std::endl;
    log.end_section();
    log.set_verbosity(initial_verbosity);
    return;
  }

  // Construct System
  parser.value = std::make_unique<System>(shared_prim, *composition_axes);
  System &system = *parser.value;

  // Parse "n_dimensions"
  log.indent() << "Parsing \"n_dimensions\": Optional[int] = 3..." << std::endl;
  parser.optional(system.n_dimensions, "n_dimensions");
  log.indent() << "- n_dimensions: " << system.n_dimensions << std::endl;

  // Parse "basis_sets"
  log << "Parsing \"basis_sets\": Optional[dict]..." << std::endl;
  if (parser.self.contains("basis_sets")) {
    auto &prim = *system.prim;
    auto &basis_sets = system.basis_sets;
    auto &basis_set_cluster_info = system.basis_set_cluster_info;
    auto &prim_neighbor_list = system.prim_neighbor_list;

    auto begin = parser.self["basis_sets"].begin();
    auto end = parser.self["basis_sets"].end();
    for (auto it = begin; it != end; ++it) {
      log.indent() << "- " << it.name() << "..." << std::endl;
      log.increase_indent();

      // parse "basis_sets"/<name>/"source"
      log.indent() << "- Getting \"source\"..." << std::endl;
      auto subparser = parser.subparse<clexulator::Clexulator>(
          fs::path("basis_sets") / it.name(), prim_neighbor_list, search_path);
      if (subparser->valid()) {
        auto clexulator = std::make_shared<clexulator::Clexulator>(
            std::move(*subparser->value));
        basis_sets.emplace(it.name(), clexulator);
        log.indent() << "- Getting \"source\": succeeded" << std::endl;
      } else {
        log.indent() << "- Getting \"source\": failed" << std::endl;
        log.decrease_indent();
        continue;
      }

      // "basis_sets"/<name>/"basis" (BasisSetClusterInfo, optional)
      if (parser.self.find_at(fs::path("basis_sets") / it.name() / "basis") !=
          parser.self.end()) {
        log.indent() << "- Getting \"basis\"..." << std::endl;
        BasisSetClusterInfo cluster_info;
        if (parse_from_file(parser,
                            fs::path("basis_sets") / it.name() / "basis",
                            search_path, cluster_info, prim)) {
          basis_set_cluster_info.emplace(
              it.name(), std::make_shared<BasisSetClusterInfo const>(
                             std::move(cluster_info)));
        }
      } else {
        log.indent() << "- No \"basis\" input found..." << std::endl;
      }
      log.decrease_indent();
    }
  } else {
    log.indent() << "- No \"basis_sets\" input found..." << std::endl;
  }

  // Parse "local_basis_sets"
  log << "Parsing \"local_basis_sets\": Optional[dict]..." << std::endl;
  if (parser.self.contains("local_basis_sets")) {
    auto &local_basis_sets = system.local_basis_sets;
    auto &local_basis_set_cluster_info = system.local_basis_set_cluster_info;
    auto &equivalents_info = system.equivalents_info;
    auto &prim_neighbor_list = system.prim_neighbor_list;
    auto &prim = *system.prim;

    // construct local Clexulator
    auto begin = parser.self["local_basis_sets"].begin();
    auto end = parser.self["local_basis_sets"].end();
    for (auto it = begin; it != end; ++it) {
      std::string local_basis_set_name = it.name();

      log.indent() << "- " << local_basis_set_name << "..." << std::endl;
      log.increase_indent();

      // parse "local_basis_sets"/<name>/"source"
      log.indent() << "- Getting \"source\"..." << std::endl;
      auto subparser = parser.subparse<std::vector<clexulator::Clexulator>>(
          fs::path("local_basis_sets") / local_basis_set_name,
          prim_neighbor_list, search_path);
      if (subparser->valid()) {
        auto local_clexulator =
            std::make_shared<std::vector<clexulator::Clexulator>>(
                std::move(*subparser->value));
        local_basis_sets.emplace(local_basis_set_name, local_clexulator);

        log.indent() << "- Getting \"source\": succeeded" << std::endl;
      } else {
        log.indent() << "- Getting \"source\": failed" << std::endl;
        log.decrease_indent();
        continue;
      }

      fs::path opt;
      opt = fs::path("local_basis_sets") / local_basis_set_name /
            "equivalents_info";

      // parse "local_basis_sets"/<name>/"equivalents_info"
      log.indent() << "- Getting \"equivalents_info\"..." << std::endl;
      auto info_subparser =
          subparse_from_file<EquivalentsInfo>(parser, opt, search_path, prim);
      if (info_subparser->valid()) {
        equivalents_info.emplace(local_basis_set_name,
                                 std::move(*info_subparser->value));
        log.indent() << "- Getting \"equivalents_info\": succeeded"
                     << std::endl;
      } else {
        log.indent() << "- Getting \"equivalents_info\": failed" << std::endl;
        log.decrease_indent();
        continue;
      }

      // "local_basis_sets"/<name>/"basis" (LocalBasisSetClusterInfo, optional)
      opt = fs::path("local_basis_sets") / local_basis_set_name / "basis";
      if (parser.self.find_at(opt) != parser.self.end()) {
        log.indent() << "- Getting \"basis\"..." << std::endl;
        LocalBasisSetClusterInfo cluster_info;
        if (parse_from_file(parser, opt, search_path, cluster_info, prim,
                            equivalents_info.at(local_basis_set_name))) {
          local_basis_set_cluster_info.emplace(
              local_basis_set_name,
              std::make_shared<LocalBasisSetClusterInfo const>(
                  std::move(cluster_info)));
        }
      } else {
        log.indent() << "- No \"basis\" input found..." << std::endl;
      }

      // "local_basis_sets"/<name>/"local_orbit_composition" (optional)
      opt = fs::path("local_basis_sets") / local_basis_set_name /
            "local_orbit_composition";
      if (parser.self.find_at(opt) != parser.self.end()) {
        log.indent() << "- Getting \"local_orbit_composition\"..." << std::endl;
        auto key_begin = parser.self.at(opt).begin();
        auto key_end = parser.self.at(opt).end();
        for (auto key_it = key_begin; key_it != key_end; ++key_it) {
          log.increase_indent();
          log.indent() << "- " << key_it.name() << "..." << std::endl;

          // parse local-orbit composition calculator data
          auto subparser = parser.subparse<LocalOrbitCompositionCalculatorData>(
              opt / key_it.name(), local_basis_set_name);
          if (!subparser->valid()) {
            log.indent() << "- " << key_it.name() << ": failed" << std::endl;
            log.decrease_indent();
            continue;
          } else {
            log.indent() << "- " << key_it.name() << ": succeeded" << std::endl;
          }
          system.local_orbit_composition_calculator_data.emplace(
              key_it.name(),
              std::shared_ptr<LocalOrbitCompositionCalculatorData>(
                  std::move(subparser->value)));
          log.decrease_indent();
        }
      } else {
        log.indent() << "- No \"local_orbit_composition\" input found..."
                     << std::endl;
      }

      log.decrease_indent();
    }
  } else {
    log.indent() << "- No \"local_basis_sets\" input found..." << std::endl;
  }

  // Parse "clex"
  log << "Parsing \"clex\": Optional[dict]..." << std::endl;
  if (parser.self.contains("clex")) {
    auto &prim = *system.prim;
    auto &clex_data = system.clex_data;
    auto const &basis_sets = system.basis_sets;

    auto begin = parser.self["clex"].begin();
    auto end = parser.self["clex"].end();
    for (auto it = begin; it != end; ++it) {
      log.indent() << "- " << it.name() << "..." << std::endl;
      log.increase_indent();

      ClexData curr;
      fs::path clex_path = fs::path("clex") / it.name();

      // "clex"/<name>/"basis_set"
      log.indent() << "- Getting \"basis_set\"..." << std::endl;
      if (!parse_and_validate_basis_set_name(parser, clex_path / "basis_set",
                                             curr.basis_set_name, basis_sets)) {
        log.indent() << "- invalid basis set name" << std::endl;
        log.decrease_indent();

        continue;
      } else {
        log.indent() << "- valid basis set name: " << curr.basis_set_name
                     << std::endl;
      }

      // get basis set cluster info if available
      log.indent() << "- Adding cluster_info (if " << curr.basis_set_name
                   << " has \"basis\" input)..." << std::endl;
      auto find_it = system.basis_set_cluster_info.find(curr.basis_set_name);
      if (find_it != system.basis_set_cluster_info.end()) {
        curr.cluster_info = find_it->second;
        log.indent() << "- Added cluster_info" << std::endl;
      } else {
        log.indent() << "- No cluster_info found" << std::endl;
      }

      // "clex"/<name>/"coefficients" (SparseCoefficients)
      log.indent() << "- Getting \"coefficients\"..." << std::endl;
      if (!parse_from_file(parser, clex_path / "coefficients", search_path,
                           curr.coefficients)) {
        log.indent() << "- Getting \"coefficients\": failed" << std::endl;
        log.decrease_indent();
        continue;
      } else {
        log.indent() << "- Getting \"coefficients\": succeeded" << std::endl;
      }

      clex_data.emplace(it.name(), curr);

      log.decrease_indent();
    }
  } else {
    log.indent() << "- No \"clex\" input found..." << std::endl;
  }

  // Parse "multiclex"
  log << "Parsing \"multiclex\": Optional[dict]..." << std::endl;
  if (parser.self.contains("multiclex")) {
    auto &prim = *system.prim;
    auto &multiclex_data = system.multiclex_data;
    auto const &basis_sets = system.basis_sets;

    auto begin = parser.self["multiclex"].begin();
    auto end = parser.self["multiclex"].end();
    for (auto it = begin; it != end; ++it) {
      log.indent() << "- " << it.name() << "..." << std::endl;
      log.increase_indent();

      MultiClexData curr;
      fs::path clex_path = fs::path("multiclex") / it.name();

      // "multiclex"/<name>/"basis_set"
      log.indent() << "- Getting \"basis_set\"..." << std::endl;
      if (!parse_and_validate_basis_set_name(parser, clex_path / "basis_set",
                                             curr.basis_set_name, basis_sets)) {
        log.indent() << "- Invalid basis set name" << std::endl;
        log.decrease_indent();
        continue;
      } else {
        log.indent() << "- Valid basis set name: " << curr.basis_set_name
                     << std::endl;
      }

      // get basis set cluster info if available
      log.indent() << "- Adding cluster_info (if " << curr.basis_set_name
                   << " has \"basis\" input)..." << std::endl;
      auto find_it = system.basis_set_cluster_info.find(curr.basis_set_name);
      if (find_it != system.basis_set_cluster_info.end()) {
        curr.cluster_info = find_it->second;
        log.indent() << "- Added cluster_info" << std::endl;
      } else {
        log.indent() << "- No cluster_info" << std::endl;
      }

      // "multiclex"/<name>/"coefficients"
      log.indent() << "- Getting \"coefficients\"..." << std::endl;
      if (!parse_from_files_object(parser, clex_path / "coefficients",
                                   search_path, curr.coefficients,
                                   curr.coefficients_glossary)) {
        log.indent() << "- Getting \"coefficients\": failed" << std::endl;
        log.decrease_indent();
        continue;
      } else {
        log.indent() << "- Getting \"coefficients\": succeeded" << std::endl;
      }

      multiclex_data.emplace(it.name(), curr);
      log.decrease_indent();
    }
  } else {
    log.indent() << "- No \"multiclex\" input found..." << std::endl;
  }

  // Parse "local_clex"
  log << "Parsing \"local_clex\": Optional[dict]..." << std::endl;
  if (parser.self.contains("local_clex")) {
    auto &local_clex_data = system.local_clex_data;
    auto const &local_basis_sets = system.local_basis_sets;

    auto begin = parser.self["local_clex"].begin();
    auto end = parser.self["local_clex"].end();
    for (auto it = begin; it != end; ++it) {
      log.indent() << "- " << it.name() << "..." << std::endl;
      log.increase_indent();

      LocalClexData curr;
      fs::path clex_path = fs::path("local_clex") / it.name();

      // "local_clex"/<name>/"local_basis_set"
      log.indent() << "- Getting \"local_basis_set\"..." << std::endl;
      if (!parse_and_validate_basis_set_name(
              parser, clex_path / "local_basis_set", curr.local_basis_set_name,
              local_basis_sets)) {
        log.indent() << "- Invalid local basis set name" << std::endl;
        log.decrease_indent();
        continue;
      } else {
        log.indent() << "- Valid local basis set name: "
                     << curr.local_basis_set_name << std::endl;
      }

      // get local basis set cluster info if available
      log.indent() << "- Adding cluster_info (if " << curr.local_basis_set_name
                   << " has \"basis\" input)..." << std::endl;
      auto find_it =
          system.local_basis_set_cluster_info.find(curr.local_basis_set_name);
      if (find_it != system.local_basis_set_cluster_info.end()) {
        curr.cluster_info = find_it->second;
        log.indent() << "- Added cluster_info" << std::endl;
      } else {
        log.indent() << "- No cluster_info" << std::endl;
      }

      // "local_clex"/<name>/"coefficients"
      log.indent() << "- Getting \"coefficients\"..." << std::endl;
      if (!parse_from_file(parser, clex_path / "coefficients", search_path,
                           curr.coefficients)) {
        log.indent() << "- Getting \"coefficients\": failed" << std::endl;
        log.decrease_indent();
        continue;
      } else {
        log.indent() << "- Getting \"coefficients\": succeeded" << std::endl;
      }

      local_clex_data.emplace(it.name(), curr);
      log.decrease_indent();
    }
  } else {
    log.indent() << "- No \"local_clex\" input found..." << std::endl;
  }

  // Parse "local_multiclex"
  log << "Parsing \"local_multiclex\": Optional[dict]..." << std::endl;
  if (parser.self.contains("local_multiclex")) {
    auto &local_multiclex_data = system.local_multiclex_data;
    auto const &local_basis_sets = system.local_basis_sets;

    auto begin = parser.self["local_multiclex"].begin();
    auto end = parser.self["local_multiclex"].end();
    for (auto it = begin; it != end; ++it) {
      log.indent() << "- " << it.name() << "..." << std::endl;
      log.increase_indent();

      LocalMultiClexData curr;
      fs::path clex_path = fs::path("local_multiclex") / it.name();

      // "local_multiclex"/<name>/"local_basis_set"
      log.indent() << "- Getting \"local_basis_set\"..." << std::endl;
      if (!parse_and_validate_basis_set_name(
              parser, clex_path / "local_basis_set", curr.local_basis_set_name,
              local_basis_sets)) {
        log.indent() << "- Invalid local basis set name" << std::endl;
        log.decrease_indent();
        continue;
      } else {
        log.indent() << "- Valid local basis set name: "
                     << curr.local_basis_set_name << std::endl;
      }

      // get local basis set cluster info if available
      log.indent() << "- Adding cluster_info (if " << curr.local_basis_set_name
                   << " has \"basis\" input)..." << std::endl;
      auto find_it =
          system.local_basis_set_cluster_info.find(curr.local_basis_set_name);
      if (find_it != system.local_basis_set_cluster_info.end()) {
        curr.cluster_info = find_it->second;
        log.indent() << "- Added cluster_info" << std::endl;
      } else {
        log.indent() << "- No cluster_info" << std::endl;
      }

      // "local_multiclex"/<name>/"coefficients"
      log.indent() << "- Getting \"coefficients\"..." << std::endl;
      if (!parse_from_files_object(parser, clex_path / "coefficients",
                                   search_path, curr.coefficients,
                                   curr.coefficients_glossary)) {
        log.indent() << "- Getting \"coefficients\": failed" << std::endl;
        log.decrease_indent();
        continue;
      } else {
        log.indent() << "- Getting \"coefficients\": succeeded" << std::endl;
      }

      local_multiclex_data.emplace(it.name(), curr);
      log.decrease_indent();
    }
  } else {
    log.indent() << "- No \"local_multiclex\" input found..." << std::endl;
  }

  // Parse "event_system"
  log << "Parsing \"event_system\": Optional[str]..." << std::endl;
  if (parser.self.contains("event_system")) {
    log.indent() << "- Parsing \"event_system\" from file..." << std::endl;
    auto const &basicstructure = system.prim->basicstructure;
    auto event_system_subparser = subparse_from_file<occ_events::OccSystem>(
        parser, "event_system", search_path, basicstructure);
    if (event_system_subparser->valid()) {
      system.event_system = std::make_shared<occ_events::OccSystem>(
          std::move(*event_system_subparser->value));
      log.indent() << "- Parsing \"event_system\": succeeded" << std::endl;
    } else {
      log.indent() << "- Parsing \"event_system\": failed" << std::endl;
    }
  } else {
    log.indent() << "- No \"event_system\" input found..." << std::endl;
  }

  // Parse "canonical_swaps", "semigrand_canonical_swaps",
  // "semigrand_canonical_multiswaps"

  // Canonical swaps
  log << "Parsing \"canonical_swaps\": Optional[list[dict]]..." << std::endl;
  log.indent() << "- There are " << system.canonical_swaps.size()
               << " \"canonical_swaps\" by default" << std::endl;
  if (parser.self.contains("canonical_swaps")) {
    auto subparser = parser.subparse_with<std::vector<monte::OccSwap>>(
        parse_array<monte::OccSwap, monte::Conversions const &>,
        "canonical_swaps", system.convert);
    if (subparser->valid()) {
      system.canonical_swaps = std::move(*subparser->value);
      log.indent() << "- Parsing \"canonical_swaps\": succeeded" << std::endl;
      log.indent() << "- Restricting system to "
                   << system.canonical_swaps.size() << "\"canonical_swaps\""
                   << std::endl;
    } else {
      log.indent() << "- Parsing \"canonical_swaps\": failed" << std::endl;
    }
  } else {
    log.indent() << "- No \"canonical_swaps\" input found..." << std::endl;
  }
  log.indent() << "- System has " << system.semigrand_canonical_swaps.size()
               << " \"canonical_swaps\"" << std::endl;

  // Semi-grand canonical swaps
  log << "Parsing \"semigrand_canonical_swaps\": Optional[list[dict]]..."
      << std::endl;
  log.indent() << "- There are " << system.semigrand_canonical_swaps.size()
               << " \"semigrand_canonical_swaps\" by default" << std::endl;
  if (parser.self.contains("semigrand_canonical_swaps")) {
    auto subparser = parser.subparse_with<std::vector<monte::OccSwap>>(
        parse_array<monte::OccSwap, monte::Conversions const &>,
        "semigrand_canonical_swaps", system.convert);
    if (subparser->valid()) {
      system.semigrand_canonical_swaps = std::move(*subparser->value);
      log.indent() << "- Parsing \"semigrand_canonical_swaps\": succeeded"
                   << std::endl;
      log.indent() << "- Restricting system to "
                   << system.semigrand_canonical_swaps.size()
                   << "\"semigrand_canonical_swaps\"" << std::endl;
    } else {
      log.indent() << "- Parsing \"semigrand_canonical_swaps\": failed"
                   << std::endl;
    }
  } else {
    log.indent() << "- No \"semigrand_canonical_swaps\" input found..."
                 << std::endl;
    log.indent() << "- System has " << system.semigrand_canonical_swaps.size()
                 << " \"semigrand_canonical_swaps\"" << std::endl;
  }

  // Semi-grand canonical multiswaps
  log << "Parsing \"semigrand_canonical_multiswaps\": Union[list[dict], int, "
         "None]..."
      << std::endl;
  log.indent() << "- There are " << system.semigrand_canonical_multiswaps.size()
               << " \"semigrand_canonical_multiswaps\" by default" << std::endl;
  if (parser.self.contains("semigrand_canonical_multiswaps")) {
    jsonParser const &tjson = parser.self["semigrand_canonical_multiswaps"];
    if (tjson.is_array()) {
      auto subparser = parser.subparse_with<std::vector<monte::MultiOccSwap>>(
          parse_array<monte::MultiOccSwap, monte::Conversions const &>,
          "semigrand_canonical_multiswaps", system.convert);
      if (subparser->valid()) {
        system.semigrand_canonical_multiswaps = std::move(*subparser->value);
        log.indent()
            << "- Parsing \"semigrand_canonical_multiswaps\": succeeded"
            << std::endl;
      } else {
        log.indent() << "- Parsing \"semigrand_canonical_swaps\": failed"
                     << std::endl;
      }
    } else if (tjson.is_int()) {
      int max_total_count = tjson.get<int>();

      log.indent() << "- Generating semi-grand canonical multiswaps with up to "
                   << max_total_count << " single semi-grand canonical swaps..."
                   << std::endl;
      system.semigrand_canonical_multiswaps = monte::make_multiswaps(
          system.semigrand_canonical_swaps, max_total_count);
      log.indent() << "- System has "
                   << system.semigrand_canonical_multiswaps.size()
                   << " semi-grand canonical multiswaps" << std::endl;
    } else if (tjson.is_null()) {
      log.indent() << "- Clearing semi-grand canonical multiswaps..."
                   << std::endl;
      system.semigrand_canonical_multiswaps.clear();
    } else {
      log.indent() << "- Invalid \"semigrand_canonical_multiswaps\" input type"
                   << std::endl;
      parser.insert_error(
          "semigrand_canonical_multiswaps",
          "If present, must be an array of MultiOccSwap, an integer indicating "
          "multiswaps should be generated using this value for the maximum "
          "number of single swaps per multiswap, or null.");
    }
  } else {
    log.indent() << "- No \"semigrand_canonical_multiswaps\" input found..."
                 << std::endl;
  }
  log.indent() << "- System has "
               << system.semigrand_canonical_multiswaps.size()
               << " \"semigrand_canonical_multiswaps\"" << std::endl;

  // Parse "kmc_events"
  log << "Parsing \"kmc_events\": Optional[dict]..." << std::endl;
  if (parser.self.contains("kmc_events")) {
    if (system.event_system == nullptr) {
      log << "- Error: \"event_system\" is required to parse \"kmc_events\""
          << std::endl;
      parser.insert_error("kmc_events",
                          "event_system is required to parse events");
    } else {
      log.indent()
          << "- Getting clex \"formation_energy\" (required for KMC)..."
          << std::endl;
      auto const &clex_data = get_clex_data(system, "formation_energy");
      if (!clex_data.cluster_info) {
        log.indent() << "- Error: no \"basis\" input for '"
                     << clex_data.basis_set_name << "' (required for KMC)"
                     << std::endl;
        std::stringstream ss;
        ss << "Warning: no \"basis\" input for basis_set '"
           << clex_data.basis_set_name << "' (required for KMC)" << std::endl;
        parser.insert_warning("kmc_events", ss.str());
      } else {
        log.indent() << "- Found cluster_info from '"
                     << clex_data.basis_set_name
                     << "' \"basis\" input (required for KMC)" << std::endl;
      }

      // check that if "entropy" is provided it also has the `basis` input
      {
        auto it = system.clex_data.find("entropy");
        if (it != system.clex_data.end() && !it->second.cluster_info) {
          log.indent() << "- Error: no \"basis\" input for '"
                       << it->second.basis_set_name << "' (required for KMC)"
                       << std::endl;
          std::stringstream ss;
          ss << "Warning: no \"basis\" input for basis_set '"
             << it->second.basis_set_name << "' (required for KMC)"
             << std::endl;
          parser.insert_warning("kmc_events", ss.str());
        }
      }

      // parse "kmc_events"/<name>
      auto const &basicstructure = system.prim->basicstructure;
      auto begin = parser.self["kmc_events"].begin();
      auto end = parser.self["kmc_events"].end();
      for (auto it = begin; it != end; ++it) {
        log.indent() << "- " << it.name() << "..." << std::endl;
        log.increase_indent();
        parse_event(parser, fs::path("kmc_events") / it.name(), search_path,
                    system.event_type_data, system.local_multiclex_data,
                    *basicstructure, system.local_basis_sets,
                    system.local_basis_set_cluster_info,
                    system.equivalents_info, system.occevent_symgroup_rep,
                    *system.event_system);
        log.decrease_indent();
      }
    }
    log.indent() << "- Found " << system.event_type_data.size()
                 << " event types" << std::endl;

    // For each event, construct default local orbit composition calculator data
    log.indent()
        << "- Construct default local_orbit_composition collecting functions..."
        << std::endl;
    for (auto const &event : system.event_type_data) {
      std::string const &event_type_name = event.first;
      OccEventTypeData const &event_data = event.second;

      log.indent() << "- " << event_type_name << "..." << std::endl;
      log.increase_indent();

      log.indent() << "- Getting local_multiclex '"
                   << event_data.local_multiclex_name << "' data..."
                   << std::endl;
      LocalMultiClexData const &local_multiclex_data =
          system.local_multiclex_data.at(event_data.local_multiclex_name);
      std::string const &local_basis_set_name =
          local_multiclex_data.local_basis_set_name;
      if (local_multiclex_data.cluster_info == nullptr) {
        throw std::runtime_error(
            "Error: missing cluster_info for local basis set '" +
            local_basis_set_name + "'");
      }
      LocalBasisSetClusterInfo const &cluster_info =
          *local_multiclex_data.cluster_info;
      auto const &orbits = cluster_info.orbits;

      if (!orbits.size()) {
        log.indent() << "- Warning: no local orbits for event '"
                     << event_type_name +
                            "', so local orbit composition calculators cannot "
                            "be constructed"
                     << std::endl;

        std::stringstream ss;
        ss << "Warning: no local orbits for event '" << event_type_name
           << "', so local orbit composition calculators cannot be constructed";
        parser.insert_warning("kmc_events", ss.str());
        log.decrease_indent();
        continue;
      }

      std::set<int> point_cluster_orbits =
          get_point_cluster_orbit_indices(orbits);

      if (point_cluster_orbits.size() == 0) {
        log.indent() << "- No point-cluster orbits: not constructing any "
                     << event_type_name + "-<orbit_index> collecting functions"
                     << std::endl;
      } else {
        // Local-orbit composition calculator data for each point-cluster orbit:
        // - Named <event_type_name>-<orbit_index>
        for (int i_orbit : point_cluster_orbits) {
          log.indent() << "- Construct "
                       << event_type_name + "-" + std::to_string(i_orbit) +
                              "..."
                       << std::endl;

          system.local_orbit_composition_calculator_data.emplace(
              event_type_name + "-" + std::to_string(i_orbit),
              std::make_shared<LocalOrbitCompositionCalculatorData>(
                  event_type_name, local_basis_set_name,
                  std::set<int>({i_orbit}) /* orbits_to_calculate */,
                  false /* combine_orbits */, 1000 /* max_size */));
        }
      }

      std::set<int> all_orbits = get_all_orbit_indices(orbits);
      if (!all_orbits.empty()) {
        log.indent() << "- Construct " << event_type_name + "-all..."
                     << std::endl;

        // Local-orbit composition calculator data for all orbits
        // - Named <event_type_name>-all
        system.local_orbit_composition_calculator_data.emplace(
            event_type_name + "-all",
            std::make_shared<LocalOrbitCompositionCalculatorData>(
                event_type_name, local_basis_set_name,
                all_orbits /* orbits_to_calculate */,
                false /* combine_orbits */, 10000 /* max_size */));

        log.indent() << "- Construct " << event_type_name + "-all-combined..."
                     << std::endl;

        // Local-orbit composition calculator data for the union of sites in all
        // orbits
        // - Named <event_type_name>-all-combined
        system.local_orbit_composition_calculator_data.emplace(
            event_type_name + "-all-combined",
            std::make_shared<LocalOrbitCompositionCalculatorData>(
                event_type_name, local_basis_set_name,
                all_orbits /* orbits_to_calculate */, true /* combine_orbits */,
                10000 /* max_size */));
      } else {
        log.indent() << "- All orbits are empty: not constructing "
                     << event_type_name + "-all"
                     << " or " << event_type_name + "-all-combined"
                     << std::endl;
      }

      log.decrease_indent();
    }

    // Construct prim event list
    log.indent() << "- Constructing prim_event_list..." << std::endl;
    system.prim_event_list = make_prim_event_list(system.event_type_data);
    log.indent() << "- prim_event_list size: " << system.prim_event_list.size()
                 << std::endl;
    if (system.prim_event_list.size() > 0) {
      log.indent() << "- prim_event_list: " << std::endl;
      log.increase_indent();
      for (auto const &prim_event_data : system.prim_event_list) {
        log.indent() << "- " << prim_event_data.prim_event_index << ": "
                     << name(prim_event_data) << std::endl;
      }
      log.decrease_indent();
    }
  } else {
    log.indent() << "- No \"kmc_events\" input found..." << std::endl;
  }

  // Parse "dof_spaces"
  log << "Parsing \"dof_spaces\": Optional[dict]..." << std::endl;
  if (parser.self.contains("dof_spaces")) {
    std::string key = "dof_spaces";
    auto begin = parser.self[key].begin();
    auto end = parser.self[key].end();
    for (auto it = begin; it != end; ++it) {
      log.indent() << "- " << it.name() << "..." << std::endl;
      log.increase_indent();

      std::string label = it.name();
      std::shared_ptr<InputParser<clexulator::DoFSpace>> subparser;
      if (it->is_obj()) {
        log.indent() << "- Parsing DoFSpace from JSON object..." << std::endl;
        subparser = parser.template subparse<clexulator::DoFSpace>(
            fs::path(key) / label, system.prim->basicstructure);
      } else if (it->is_string()) {
        log.indent() << "- Parsing DoFSpace from file..." << std::endl;
        subparser = subparse_from_file<clexulator::DoFSpace>(
            parser, fs::path(key) / label, search_path,
            system.prim->basicstructure);
      } else {
        log.indent() << "- Parsing DoFSpace: invalid type, must be a file name "
                        "or JSON object "
                     << std::endl;
        parser.insert_error(fs::path(key) / label,
                            "Error: must be a file name or JSON object");
        log.decrease_indent();
        continue;
      }
      if (subparser->valid()) {
        log.indent() << "- Parsing DoFSpace: succeeded" << std::endl;
        system.dof_spaces.emplace(
            label,
            std::make_shared<clexulator::DoFSpace const>(*subparser->value));
      } else {
        log.indent() << "- Parsing DoFSpace: failed" << std::endl;
      }
      log.decrease_indent();
    }
  } else {
    log.indent() << "- No \"dof_spaces\" input found..." << std::endl;
  }

  // Parse "dof_subspaces"
  log << "Parsing \"dof_subspaces\": Optional[dict[str,list[list[int]]]]..."
      << std::endl;
  parser.optional(system.dof_subspaces, "dof_subspaces");
  log.indent() << "- Found " << system.dof_subspaces.size() << " dof_subspaces"
               << std::endl;

  // Parse "additional_params"
  log << "Parsing \"additional_params\": Optional[dict[str,list[list[int]]]]..."
      << std::endl;
  system.additional_params.put_obj();
  parser.optional(system.additional_params, "additional_params");
  if (!system.additional_params.is_obj()) {
    log.indent() << "- Error:  if present, \"additional_params\" must be a "
                    "JSON object or null"
                 << std::endl;
    parser.insert_error("additional_params",
                        "Error: if present, must be a JSON object or null");
  }
  log.indent() << "- Found " << system.additional_params.size()
               << " additional parameters" << std::endl;
  if (system.additional_params.size() > 0) {
    log.indent() << "- Additional parameters:" << std::endl;
    log.increase_indent();

    auto begin = system.additional_params.begin();
    auto end = system.additional_params.end();
    for (auto it = begin; it != end; ++it) {
      log.indent() << "- " << it.name() << std::endl;
    }
    log.decrease_indent();
  }

  log.end_section();
  log.set_verbosity(initial_verbosity);
}

}  // namespace clexmonte
}  // namespace CASM
