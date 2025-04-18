#include "casm/clexmonte/monte_calculator/MonteEventData.hh"

#include "casm/monte/misc/memory_used.hh"
#include "casm/monte/sampling/io/json/SelectedEventFunctions_json_io.hh"

namespace CASM {
namespace clexmonte {

EventTypeStats::EventTypeStats(
    std::vector<std::string> const &_partion_names_by_type,
    std::vector<std::string> const &_partion_names_by_equivalent_index,
    double _initial_begin, double _bin_width, bool _is_log, Index _max_size)
    : n_total(0),
      min(0.0),
      max(0.0),
      sum(0.0),
      mean(0.0),
      hist_by_type(_partion_names_by_type, _initial_begin, _bin_width, _is_log,
                   _max_size),
      hist_by_equivalent_index(_partion_names_by_equivalent_index,
                               _initial_begin, _bin_width, _is_log, _max_size) {
}

void EventTypeStats::insert(int partition_by_type,
                            int partition_by_equivalent_index, double value) {
  n_total += 1;
  min = std::min(min, value);
  max = std::max(max, value);
  sum += value;
  mean = sum / n_total;
  hist_by_type.insert(partition_by_type, value);
  hist_by_equivalent_index.insert(partition_by_equivalent_index, value);
}

EventDataSummary::EventDataSummary(
    std::shared_ptr<StateData> const &_state_data,
    MonteEventData const &_event_data, double energy_bin_width,
    double freq_bin_width, double rate_bin_width)
    : state_data(_state_data),
      event_data(_event_data),
      prim_event_list(event_data.prim_event_list()) {
  for (Index i = 0; i < prim_event_list.size(); ++i) {
    auto type = type_key(i);
    auto equiv = equiv_key(i);

    all_types.insert(type);
    equiv_keys_by_type[type].insert(equiv);
    all_equiv_keys.insert(equiv);

    n_possible.by_type[type] = 0.0;
    n_possible.by_equivalent_index[equiv] = 0.0;

    n_allowed.by_type[type] = 0.0;
    n_allowed.by_equivalent_index[equiv] = 0.0;

    n_abnormal.by_type[type] = 0.0;
    n_abnormal.by_equivalent_index[equiv] = 0.0;

    rate.by_type[type] = 0.0;
    rate.by_equivalent_index[equiv] = 0.0;

    //    n_impact.by_type[type] = 0.0;
    //    n_impact.by_equivalent_index[equiv] = 0.0;

    auto &_impact_by_type = impact_table.by_type[type];
    auto &_impact_by_equiv = impact_table.by_equivalent_index[equiv];
    for (Index j = 0; j < prim_event_list.size(); ++j) {
      _impact_by_type[type_key(j)] = 0.0;
      _impact_by_equiv[equiv_key(j)] = 0.0;
    }
  }

  n_events_allowed = 0;
  n_events_possible = 0;
  n_abnormal_total = 0;
  event_list_size = event_data.event_list().size();
  total_rate = event_data.event_list().total_rate();
  mean_time_increment = 1.0 / total_rate;

  resident_bytes_used = memory_used(true);
  resident_MiB_used = memory_used_MiB(true);

  SelectedEventInfo info(event_data.prim_event_list());
  info.make_indices_by_type();
  to_event_type = *info.prim_event_index_to_index;
  event_type_names = info.partition_names;

  info.make_indices_by_equivalent_index();
  to_equivalent_index = *info.prim_event_index_to_index;
  equivalent_index_names = info.partition_names;

  stats_labels.push_back("dE_final");
  stats.emplace_back(event_type_names, equivalent_index_names,
                     0.0 /* initial_begin */, energy_bin_width /* bin_width */,
                     false /* is_log */);

  stats_labels.push_back("dE_activated");
  stats.emplace_back(event_type_names, equivalent_index_names,
                     0.0 /* initial_begin */, energy_bin_width /* bin_width */,
                     false /* is_log */);

  stats_labels.push_back("Ef_kra");
  stats.emplace_back(event_type_names, equivalent_index_names,
                     0.0 /* initial_begin */, energy_bin_width /* bin_width */,
                     false /* is_log */);

  stats_labels.push_back("freq");
  stats.emplace_back(event_type_names, equivalent_index_names,
                     0.0 /* initial_begin */, freq_bin_width /* bin_width */,
                     true /* is_log */);

  stats_labels.push_back("rate");
  stats.emplace_back(event_type_names, equivalent_index_names,
                     0.0 /* initial_begin */, rate_bin_width /* bin_width */,
                     true /* is_log */);

  for (EventID const &id : event_data.event_list()) {
    EventState const &state = event_data.event_state(id);
    _add_count(id, state);
    _add_impact(id, state);
    _add_stats(id, state);
  }
}

void EventDataSummary::_add_count(EventID const &id, EventState const &state) {
  auto type = type_key(id);
  auto equiv = equiv_key(id);

  // double increment = 1.0 / state_data->n_unitcells;
  Index increment = 1;

  if (prim_event_list[id.prim_event_index].is_forward) {
    n_possible.by_type[type] += increment;
    n_possible.by_equivalent_index[equiv] += increment;
    n_events_possible += increment;
  }
  if (state.is_allowed) {
    n_events_allowed += increment;
    n_allowed.by_type[type] += increment;
    n_allowed.by_equivalent_index[equiv] += increment;
    if (!state.is_normal) {
      n_abnormal.by_type[type] += increment;
      n_abnormal.by_equivalent_index[equiv] += increment;
      n_abnormal_total += increment;
    }
  }
  rate.by_type[type] += state.rate;
  rate.by_equivalent_index[equiv] += state.rate;
}

void EventDataSummary::_add_impact(EventID const &id, EventState const &state) {
  TypeKey type = type_key(id);
  EquivKey equiv = equiv_key(id);

  // double increment = 1.0 / state_data->n_unitcells;
  Index increment = 1;

  if (state.is_allowed) {
    //    n_impact.by_type[type] += increment;
    //    n_impact.by_equivalent_index[equiv] += increment;

    // -- impact_table.by_type --
    for (EventID const &impact_id : event_data.event_impact(id)) {
      TypeKey impact_type = type_key(impact_id);
      impact_table.by_type[type][impact_type] += increment;
    }

    // -- impact_table.by_equivalent_index --
    for (EventID const &impact_id : event_data.event_impact(id)) {
      EquivKey impact_equiv = equiv_key(impact_id);
      impact_table.by_equivalent_index[equiv][impact_equiv] += increment;
    }
  }

  // -- neighborhood_size_total --
  if (neighborhood_size_total.count(type) == 0) {
    PrimEventData const &prim_event_data =
        event_data.prim_event_list()[id.prim_event_index];
    System const &system = *event_data.system();
    OccEventTypeData const &event_type_data =
        get_event_type_data(system, prim_event_data.event_type_name);
    clust::IntegralCluster phenom = make_cluster(prim_event_data.event);

    std::set<xtal::UnitCellCoord> total_hood;
    ClexData const &clex_data = get_clex_data(system, "formation_energy");
    expand(phenom, total_hood, *clex_data.cluster_info, clex_data.coefficients);
    neighborhood_size_formation_energy[type] = total_hood.size();

    // local basis set dependence
    LocalMultiClexData const &local_multiclex_data =
        get_local_multiclex_data(system, event_type_data.local_multiclex_name);
    std::set<xtal::UnitCellCoord> kra_hood = get_required_update_neighborhood(
        system, local_multiclex_data, prim_event_data.equivalent_index, "kra");
    neighborhood_size_kra[type] = kra_hood.size();
    total_hood.insert(kra_hood.begin(), kra_hood.end());

    std::set<xtal::UnitCellCoord> freq_hood = get_required_update_neighborhood(
        system, local_multiclex_data, prim_event_data.equivalent_index, "freq");
    neighborhood_size_freq[type] = freq_hood.size();
    total_hood.insert(freq_hood.begin(), freq_hood.end());

    neighborhood_size_total[type] = total_hood.size();
  }
}

void EventDataSummary::_add_stats(EventID const &id, EventState const &state) {
  if (!state.is_allowed) {
    return;
  }
  int t = to_event_type[id.prim_event_index];
  int e = to_equivalent_index[id.prim_event_index];

  if (state.freq <= 0.0) {
    std::stringstream msg;
    std::string event_type_name = type_key(id);

    msg << "Error in EventDataSummary: ";
    msg << "state.is_allowed=true && state.freq <= 0.0, ";
    msg << "for event=" << event_type_name << "." << e << " ";
    msg << "(state.freq=" << state.freq << ") ";
    throw std::runtime_error(msg.str());
  }

  if (state.rate <= 0.0) {
    std::stringstream msg;
    std::string event_type_name = type_key(id);

    msg << "Error in EventDataSummary: ";
    msg << "state.is_allowed=true && state.rate <= 0.0, ";
    msg << "for event=" << event_type_name << "." << e << " ";
    msg << "(state.rate=" << state.rate << ")";
    throw std::runtime_error(msg.str());
  }

  // order determined by constructor
  int i = 0;
  stats[i++].insert(t, e, state.dE_final);
  stats[i++].insert(t, e, state.dE_activated);
  stats[i++].insert(t, e, state.Ef_kra);
  stats[i++].insert(t, e, state.freq);
  stats[i++].insert(t, e, state.rate);
}

}  // namespace clexmonte

jsonParser &to_json(clexmonte::EventTypeStats const &stats, jsonParser &json) {
  json.put_obj();
  json["n_total"] = stats.n_total;
  json["min"] = stats.min;
  json["max"] = stats.max;
  json["sum"] = stats.sum;
  json["mean"] = stats.mean;
  json["by_type"] = stats.hist_by_type;
  json["by_equivalent_index"] = stats.hist_by_equivalent_index;
  return json;
}

jsonParser &to_json(clexmonte::EventDataSummary::IntCountByType const &count,
                    jsonParser &json) {
  json.put_obj();
  {
    jsonParser &y = json["by_type"];
    for (auto const &pair : count.by_type) {
      y[pair.first] = pair.second;
    }
  }
  {
    jsonParser &y = json["by_equivalent_index"];
    for (auto const &pair : count.by_equivalent_index) {
      auto const &key = pair.first;
      y[key.first + "." + std::to_string(key.second)] = pair.second;
    }
  }
  return json;
}

jsonParser &to_json(clexmonte::EventDataSummary::FloatCountByType const &count,
                    jsonParser &json) {
  json.put_obj();
  {
    jsonParser &y = json["by_type"];
    for (auto const &pair : count.by_type) {
      y[pair.first] = pair.second;
    }
  }
  {
    jsonParser &y = json["by_equivalent_index"];
    for (auto const &pair : count.by_equivalent_index) {
      auto const &key = pair.first;
      y[key.first + "." + std::to_string(key.second)] = pair.second;
    }
  }
  return json;
}

jsonParser &to_json(clexmonte::EventDataSummary const &event_data_summary,
                    jsonParser &json) {
  auto const &x = event_data_summary;
  json.put_obj();

  // Number of unit cells
  json["n_unitcells"] = x.state_data->n_unitcells;

  // Number of each event type
  to_json(x.n_allowed, json["n_events"]);

  // Total number of allowed events
  json["n_events"]["total"] = x.n_events_allowed;

  // Event list size
  json["event_list_size"] = x.event_list_size;

  // Number of events without barrier, by type
  to_json(x.n_abnormal, json["n_abnormal_events"]);

  json["n_abnormal_events"]["total"] = x.n_abnormal_total;

  // Event rate sum, by type
  to_json(x.rate, json["rate"]);

  // Total event rate sum
  json["rate"]["total"] = x.total_rate;

  // Mean time increment, by type
  clexmonte::EventDataSummary::FloatCountByType mean_time_increment;
  for (auto const &pair : x.rate.by_type) {
    mean_time_increment.by_type[pair.first] = 1.0 / pair.second;
  }
  for (auto const &pair : x.rate.by_equivalent_index) {
    mean_time_increment.by_equivalent_index[pair.first] = 1.0 / pair.second;
  }
  to_json(mean_time_increment, json["mean_time_increment"]);

  // Mean time increment total
  json["mean_time_increment"]["total"] = x.mean_time_increment;

  // Memory used (in MiB)
  json["memory_used_MiB"] = x.resident_MiB_used;

  // Memory used (str)
  json["memory_used"] = convert_size(x.resident_bytes_used);

  // Impact neighborhood
  {
    jsonParser &y = json["impact_neighborhood"];
    for (auto const &pair : x.neighborhood_size_total) {
      std::string event_type_name = pair.first;
      jsonParser &z = json["impact_neighborhood"][event_type_name];
      z["total"] = pair.second;
      z["formation_energy"] =
          x.neighborhood_size_formation_energy.at(event_type_name);
      z["kra"] = x.neighborhood_size_kra.at(event_type_name);
      z["freq"] = x.neighborhood_size_freq.at(event_type_name);
    }
  }

  //  // Impact number
  //  to_json(x.n_impact, json["impact_number"]);

  // Impact table
  {
    jsonParser &y = json["impact_number"];

    y["type"] = jsonParser::array();
    for (auto const &type : x.all_types) {
      y["type"].push_back(type);
    }
    y["by_type"] = jsonParser::array();
    for (auto const &pair : x.impact_table.by_type) {
      jsonParser j = jsonParser::array();
      for (auto const &pair2 : pair.second) {
        j.push_back(pair2.second / x.n_allowed.by_type.at(pair.first));
      }
      y["by_type"].push_back(j);
    }

    y["equiv"] = jsonParser::array();
    for (auto &equiv : x.all_equiv_keys) {
      y["equiv"].push_back(equiv.first + "." + std::to_string(equiv.second));
    }
    y["by_equiv"] = jsonParser::array();
    for (auto const &pair : x.impact_table.by_equivalent_index) {
      jsonParser j = jsonParser::array();
      for (auto const &pair2 : pair.second) {
        j.push_back(pair2.second /
                    x.n_allowed.by_equivalent_index.at(pair.first));
      }
      y["by_equiv"].push_back(j);
    }
  }

  // Stats
  json["stats"] = jsonParser::object();
  for (int i = 0; i < x.stats_labels.size(); ++i) {
    json["stats"][x.stats_labels[i]] = x.stats[i];
  }

  return json;
}

}  // namespace CASM