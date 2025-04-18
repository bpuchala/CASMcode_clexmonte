#include "casm/clexmonte/events/io/stream/EventState_stream_io.hh"

#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/clexmonte/events/event_data.hh"

namespace CASM {
namespace clexmonte {

void print(Log &log, EventState const &event_state) {
  log.indent() << "is_allowed: " << std::boolalpha << event_state.is_allowed
               << std::endl;
  if (event_state.is_allowed) {
    log.indent() << "dE_activated: " << event_state.dE_activated << std::endl;
    log.indent() << "dE_final: " << event_state.dE_final << std::endl;
    log.indent() << "is_normal: " << std::boolalpha << event_state.is_normal
                 << std::endl;
    log.indent() << "Ef_kra: " << event_state.Ef_kra << std::endl;
    log.indent() << "freq: " << event_state.freq << std::endl;
    log.indent() << "rate: " << event_state.rate << std::endl;
  }
}

void print(Log &log, EventState const &event_state,
           PrimEventData const &prim_event_data) {
  log.indent() << "prim_event_index: " << prim_event_data.prim_event_index
               << std::endl;
  log.indent() << "event_type_name: " << prim_event_data.event_type_name
               << std::endl;
  log.indent() << "equivalent_index: " << prim_event_data.equivalent_index
               << std::endl;
  log.indent() << "is_forward: " << std::boolalpha << prim_event_data.is_forward
               << std::endl;
  log.indent() << "occ_init: " << prim_event_data.occ_init << std::endl;
  log.indent() << "occ_final: " << prim_event_data.occ_final << std::endl;
  print(log.indent(), event_state);
}

void print(Log &log, EventState const &event_state, EventData const &event_data,
           PrimEventData const &prim_event_data) {
  log.indent() << "prim_event_index: " << prim_event_data.prim_event_index
               << std::endl;
  log.indent() << "unitcell_index: " << event_data.unitcell_index << std::endl;
  log.indent() << "event_type_name: " << prim_event_data.event_type_name
               << std::endl;
  log.indent() << "equivalent_index: " << prim_event_data.equivalent_index
               << std::endl;
  log.indent() << "is_forward: " << std::boolalpha << prim_event_data.is_forward
               << std::endl;
  log.indent() << "linear_site_index: " << event_data.event.linear_site_index
               << std::endl;
  log.indent() << "occ_init: " << prim_event_data.occ_init << std::endl;
  log.indent() << "occ_final: " << prim_event_data.occ_final << std::endl;
  print(log.indent(), event_state);
}

}  // namespace clexmonte
}  // namespace CASM
