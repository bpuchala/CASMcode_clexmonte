import json
import pathlib

from libcasm.local_configuration import (
    LocalConfiguration,
    LocalConfigurationList,
    OccEventSymInfo,
)

from ._clexmonte_system import (
    System,
)


def make_system_event_info(
    system: System,
):
    """Make event info (:class:`libcasm.local_configuration.OccEventSymInfo`) for all
    events in a system

    Parameters
    ----------
    system : libcasm.clexmonte.System

    Returns
    -------
    event_info: dict[str, libcasm.local_configuration.OccEventSymInfo]
        A dictionary of event info
        (:class:`libcasm.local_configuration.OccEventSymInfo`) for each event type in
        the system, with the event type name as the key.
    """
    event_info = dict()
    event_system = system.event_system
    for event_type_name in system.event_type_names:
        (
            phenomenal_clusters,
            equivalent_generating_op_indices,
            translations,
        ) = system.equivalents_info(event_type_name)

        event_info[event_type_name] = OccEventSymInfo.init(
            prim=system.prim,
            system=event_system,
            prototype_event=system.prototype_event(event_type_name),
            phenomenal_clusters=phenomenal_clusters,
            equivalent_generating_op_indices=equivalent_generating_op_indices,
        )

    return event_info


def read_abnormal_events(
    path: pathlib.Path,
    system: System,
    event_info: dict[str, OccEventSymInfo],
):
    """Read local configurations and event data from one or more abnormal events jsonl
    files

    Parameters
    ----------
    path : list[pathlib.Path]
        The path to one or more jsonl files containing the local configurations and
        event state data.
    system : libcasm.clexmonte.System
        The system
    event_info: dict[str, libcasm.local_configuration.OccEventSymInfo]
        A dictionary of event info
        (:class:`libcasm.local_configuration.OccEventSymInfo`) for each event type in
        the system, with the event type name as the key.

    Returns
    -------
    n: Optional[int]
        The number of local configurations in the file, or None if the files
        do not exist. This is the sum of the number of local configurations of
        each event type.
    local_configurations : dict[str, libcasm.local_configuration.LocalConfigurationList]
        Local configurations
        (:class:`libcasm.local_configuration.LocalConfigurationList`)
        for by event type, with the event type name as the key.
    event_states : dict[str, list[dict]]
        Event states, for the local configurations, by event type, with the event type
        name as the key, in the same order as the local configurations. The

        Example output, for one local configuration:

        .. code-block:: Python

            {
              "event_state": {
                "Ef_kra": 0.7375,
                "dE_activated": 1.6666666666666665,
                "dE_final": 1.6666666666666665,
                "formation_energy_delta_corr": [
                    0.0, 0.0, 0.0, -0.8333333333333333, 0.5892558333333333, 0.0, 0.0,
                    0.0, 0.0,
                ],
                "freq": 10000000000000.0,
                "is_allowed": true,
                "is_normal": false,
                "local_corr": [1.0, 0.5, 0.0, 0.5, 0.0, 0.25, 0.0, 0.5, 0.0],
                "rate": 1000704.0785393054
              },
              "linear_site_index": [212, 122],
              "prim_event_data": {
                "equivalent_index": 2,
                "event_type_name": "B_Va_1NN",
                "is_forward": false,
                "occ_final": [1, 2],
                "occ_init": [2, 1],
                "prim_event_index": 17
              },
              "unitcell_index": 212
            }

    """
    local_configurations = dict()
    event_states = dict()
    for event_type_name in event_info:
        local_configurations[event_type_name] = LocalConfigurationList(
            event_info=event_info[event_type_name],
        )
        event_states[event_type_name] = list()

    n = None
    for _p in path:
        p = pathlib.Path(_p)
        if not p.exists():
            continue
        if n is None:
            n = 0

        with open(p, "r") as f:
            for line in f:
                data = json.loads(line)
                _event = data.get("event")
                event_type_name = _event.get("prim_event_data").get("event_type_name")
                _local_configuration = LocalConfiguration.from_dict(
                    data=data.get("local_configuration"),
                    supercells=system.supercells,
                    event_info=event_info[event_type_name],
                )
                n += 1
                local_configurations[event_type_name].append(_local_configuration)
                event_states[event_type_name].append(_event)

    return (n, local_configurations, event_states)
