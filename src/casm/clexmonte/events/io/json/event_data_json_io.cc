#include "casm/clexmonte/events/io/json/event_data_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/events/event_data.hh"
#include "casm/crystallography/io/UnitCellCoordIO.hh"

namespace CASM {
namespace clexmonte {

jsonParser &to_json(EventState const &event_state, jsonParser &json) {
  json["is_allowed"] = event_state.is_allowed;
  if (event_state.is_allowed) {
    if (event_state.formation_energy_delta_corr) {
      to_json(*event_state.formation_energy_delta_corr,
              json["formation_energy_delta_corr"], jsonParser::as_array());
    }
    if (event_state.local_corr) {
      to_json(*event_state.local_corr, json["local_corr"],
              jsonParser::as_array());
    }
    json["is_normal"] = event_state.is_normal;
    json["dE_final"] = event_state.dE_final;
    json["Ef_kra"] = event_state.Ef_kra;
    json["dE_activated"] = event_state.dE_activated;
    json["freq"] = event_state.freq;
    json["rate"] = event_state.rate;
    json["dEf_final"] = event_state.dEf_final;
    json["dEf_activated"] = event_state.dEf_activated;
    json["reverse_rate"] = event_state.reverse_rate;
    json["dS_final"] = event_state.dS_final;
    json["S_kra"] = event_state.S_kra;
    json["dS_activated"] = event_state.dS_activated;
  }
  return json;
}

jsonParser &to_json(EventState const &event_state, jsonParser &json,
                    PrimEventData const &prim_event_data) {
  to_json(prim_event_data, json);
  to_json(event_state, json);
  return json;
}

jsonParser &to_json(EventState const &event_state, jsonParser &json,
                    EventData const &event_data,
                    PrimEventData const &prim_event_data) {
  json["unitcell_index"] = event_data.unitcell_index;
  json["linear_site_index"] = event_data.event.linear_site_index;
  to_json(event_state, json, prim_event_data);
  return json;
}

jsonParser &to_json(EventData const &event_data, jsonParser &json) {
  json["unitcell_index"] = event_data.unitcell_index;
  json["linear_site_index"] = event_data.event.linear_site_index;
  return json;
}

jsonParser &to_json(EventData const &event_data, jsonParser &json,
                    PrimEventData const &prim_event_data) {
  to_json(event_data, json);
  to_json(prim_event_data, json);
  return json;
}

jsonParser &to_json(PrimEventData const &prim_event_data, jsonParser &json) {
  json["prim_event_index"] = prim_event_data.prim_event_index;
  json["event_type_name"] = prim_event_data.event_type_name;
  json["equivalent_index"] = prim_event_data.equivalent_index;
  json["is_forward"] = prim_event_data.is_forward;
  json["occ_init"] = prim_event_data.occ_init;
  json["occ_final"] = prim_event_data.occ_final;
  return json;
}

jsonParser &to_json(
    clexmonte::PrimEventData const &data, jsonParser &json,
    std::optional<std::reference_wrapper<occ_events::OccSystem const>>
        event_system,
    occ_events::OccEventOutputOptions const &options) {
  json.put_obj();
  json["event_type_name"] = data.event_type_name;
  json["equivalent_index"] = data.equivalent_index;
  json["is_forward"] = data.is_forward;
  json["prim_event_index"] = data.prim_event_index;
  to_json(data.event, json["event"], event_system, options);
  json["sites"] = data.sites;
  json["occ_init"] = data.occ_init;
  json["occ_final"] = data.occ_final;
  return json;
}

jsonParser &to_json(EventImpactInfo const &impact, jsonParser &json) {
  json.put_obj();
  json["phenomenal_sites"] = impact.phenomenal_sites;
  json["required_update_neighborhood"] = impact.required_update_neighborhood;
  return json;
}

jsonParser &to_json(clexmonte::EventID const &event_id, jsonParser &json) {
  json["unitcell_index"] = event_id.unitcell_index;
  json["prim_event_index"] = event_id.prim_event_index;
  return json;
}

void parse(InputParser<clexmonte::EventID> &parser) {
  auto ptr = std::make_unique<clexmonte::EventID>();
  clexmonte::EventID &id = *ptr;
  parser.require(id.unitcell_index, "unitcell_index");
  parser.require(id.prim_event_index, "prim_event_index");
  if (parser.valid()) {
    parser.value = std::move(ptr);
  }
}

void from_json(clexmonte::EventID &event_id, jsonParser const &json) {
  InputParser<clexmonte::EventID> parser{json};
  std::stringstream ss;
  ss << "Error: Invalid clexmonte::EventID object";
  report_and_throw_if_invalid(parser, err_log(), std::runtime_error{ss.str()});
  event_id = std::move(*parser.value);
}

jsonParser &to_json(clexmonte::EventFilterGroup const &filter,
                    jsonParser &json) {
  json.put_obj();
  json["unitcell_index"] = filter.unitcell_index;
  json["include_by_default"] = filter.include_by_default;
  json["prim_event_index"] = filter.prim_event_index;
  return json;
}

void parse(InputParser<clexmonte::EventFilterGroup> &parser) {
  auto ptr = std::make_unique<clexmonte::EventFilterGroup>();
  clexmonte::EventFilterGroup &filter = *ptr;
  parser.require(filter.unitcell_index, "unitcell_index");
  parser.require(filter.include_by_default, "include_by_default");
  parser.require(filter.prim_event_index, "prim_event_index");
  if (parser.valid()) {
    parser.value = std::move(ptr);
  }
}

void from_json(clexmonte::EventFilterGroup &filter, jsonParser const &json) {
  InputParser<clexmonte::EventFilterGroup> parser{json};
  std::stringstream ss;
  ss << "Error: Invalid clexmonte::EventFilterGroup object";
  report_and_throw_if_invalid(parser, err_log(), std::runtime_error{ss.str()});
  filter = std::move(*parser.value);
}

jsonParser &to_json(
    clexmonte::SelectedEvent const &selected_event, jsonParser &json,
    std::optional<std::reference_wrapper<occ_events::OccSystem const>>
        event_system,
    occ_events::OccEventOutputOptions const &options) {
  json.put_obj();
  json["event_id"] = selected_event.event_id;
  json["event_index"] = selected_event.event_index;
  json["total_rate"] = selected_event.total_rate;
  json["time"] = selected_event.time;
  json["time_increment"] = selected_event.time_increment;

  if (selected_event.prim_event_data) {
    to_json(*selected_event.prim_event_data, json["prim_event_data"],
            event_system, options);
  }

  if (selected_event.event_data) {
    to_json(*selected_event.event_data, json["event_data"]);
  }

  if (selected_event.event_state) {
    to_json(*selected_event.event_state, json["event_state"]);
  }

  json["group"] = selected_event.group;
  json["group_state"] = selected_event.group_state;

  return json;
}

}  // namespace clexmonte
}  // namespace CASM
