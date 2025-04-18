import textwrap
from typing import Optional

from libcasm.local_configuration import (
    LocalConfiguration,
    LocalConfigurationList,
)

from ._clexmonte_monte_calculator import (
    MonteCalculatorCore,
)
from ._clexmonte_system import (
    System,
)
from ._system_methods import (
    make_system_event_info,
    read_abnormal_events,
)


class MonteCalculator(MonteCalculatorCore):
    def __init__(
        self,
        method: str,
        system: System,
        params: Optional[dict] = None,
    ):
        """
        .. rubric:: Constructor

        Parameters
        ----------
        method : str
            Monte Carlo method name. The options are:

            - "semigrand_canonical": Metropolis algorithm in the semi-grand
              canonical ensemble. Input states require `"temperature"` and
              `"param_chem_pot"` conditions.
            - "canonical": Metropolis algorithm in the canonical ensemble.
              Input states require `"temperature"` and one of
              `"param_composition"` or `"mol_composition"` conditions.
            - "kinetic": Kinetic Monte Carlo method. Input states require
              `"temperature"` and one of `"param_composition"` or
              `"mol_composition"` conditions.

        system : libcasm.clexmonte.System
            Cluster expansion model system data. The required data depends on
            the calculation method. See links under `method` for what system
            data is required for each method.

        params: Optional[dict] = None
            Monte Carlo calculation method parameters. Expected values
            depends on the calculation method.

        """
        super().__init__(
            method=method,
            system=system,
            params=params,
        )

        self.event_info = make_system_event_info(self.system)
        """dict[str, libcasm.local_configuration.OccEventSymInfo]: A dictionary of 
        event info used to construct and transform 
        :class:`~libcasm.local_configuration.LocalConfiguration`."""

    def make_local_configuration(self):
        """Make a :class:`~libcasm.local_configuration.LocalConfiguration` from the
        current state and current selected event.

        Returns
        -------
        local_configuration : libcasm.local_configuration.LocalConfiguration
            The local configuration.
        """
        selected_event = self.selected_event
        event_type_name = selected_event.prim_event_data.event_type_name
        unitcell_index = selected_event.event_data.unitcell_index
        equivalent_index = selected_event.prim_event_data.equivalent_index
        return LocalConfiguration(
            configuration=self.state_data.state.configuration,
            pos=(unitcell_index, equivalent_index),
            event_info=self.event_info[event_type_name],
        )

    def read_abnormal_events(
        self,
        which="all",
    ) -> tuple[Optional[int], dict[str, LocalConfigurationList], dict[str, list[dict]]]:
        """Read abnormal events file, if present.

        Notes
        -----
        - The `"kinetic"` MonteCalculator method has optional support for writing
          the local configuration and event states of abnormal events during a
          simulation. This is controlled by the `"abnormal_event_handling"`
          MonteCalculator parameter.
        - An "encountered_abnormal_events.jsonl" file is written for events that are
          encountered and calculated (but may not be selected).
        - A "selected_abnormal_events.jsonl" file is written for events that are
          selected.
        - The files are expected in `self.event_data.output_dir`.
        - This function reads the file(s) and returns the local configurations and
          event states.

        Examples
        --------

        Read all encountered and selected abnormal events, and print a summary:

        >>> n, local_configurations, event_data = mc.read_abnormal_events()

        Read only the encountered abnormal events:

        >>> n, local_configurations, event_data = mc.read_abnormal_events(
            which="encountered",
        )

        Read only the selected abnormal events:

        >>> n, local_configurations, event_data = mc.read_abnormal_events(
            which="selected",
        )

        Parameters
        ----------
        which : str = "all"
            The type of abnormal events to read. The options are:

            - "all": Read all abnormal events.
            - "encountered": Read only the encountered abnormal events from the
              "encountered_abnormal_events.jsonl" file, if it exists.
            - "selected": Read only the selected abnormal events from the
              "selected_abnormal_events.jsonl" file, if it exists.

        Returns
        -------
        n: Optional[int]
            The number of local configurations in the file, or None if the files
            do not exist. This is the sum of the number of local configurations of
            each event type.
        local_configurations : dict[str,
        libcasm.local_configuration.LocalConfigurationList]
            Local configurations
            (:class:`libcasm.local_configuration.LocalConfigurationList`)
            for by event type, with the event type name as the key.
        event_data : dict[str, list[dict]]
            Event info, for the local configurations, by event type, with the event type
            name as the key, in the same order as the local configurations. The

            Example event data, for one local configuration:

            .. code-block:: Python

                {
                  "event_state": {
                    "Ef_kra": 0.7375,
                    "dE_activated": 1.6666666666666665,
                    "dE_final": 1.6666666666666665,
                    "formation_energy_delta_corr": [
                        0.0, 0.0, 0.0, -0.8333333333333333, 0.5892558333333333, 0.0,
                        0.0, 0.0, 0.0,
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
        if which not in ["selected", "encountered", "all"]:
            raise ValueError(
                f"Invalid value for 'which': {which}. "
                "Must be one of ['selected', 'encountered', 'all']"
            )

        encountered_path = (
            self.event_data.output_dir / "encountered_abnormal_events.jsonl"
        )
        selected_path = self.event_data.output_dir / "selected_abnormal_events.jsonl"
        path = []
        if which in ["selected", "all"]:
            path.append(selected_path)
        if which in ["encountered", "all"]:
            path.append(encountered_path)
        return read_abnormal_events(
            path=path,
            system=self.system,
            event_info=self.event_info,
        )

    def print_selected_event_functions(self):
        """Print a summary of the selected event functions in a MonteCalculator"""

        all_functions = self.selected_event_functions

        def fill(text):
            return textwrap.fill(
                text,
                width=80,
                initial_indent="",
                subsequent_indent="    ",
            )

        def print_generic_functions(functions):
            for key, function in functions.items():
                print(key + ":")
                print(fill("  Description = " + function.description))
                print(
                    fill(
                        "  Requires event state = " + str(function.requires_event_state)
                    )
                )
                print(fill("  Default order = " + str(function.order)))
                print()

        def print_functions(functions):
            for key, function in functions.items():
                print(key + ":")
                print(fill("  Description = " + function.description))
                if hasattr(function, "shape"):
                    if len(function.shape) == 0:
                        print(fill("  Shape = [] (Scalar)"))
                    else:
                        print(fill("  Shape = " + str(function.shape)))
                        print(
                            fill("  Component names = " + str(function.component_names))
                        )
                if hasattr(function, "partition_names"):
                    print(fill("  Partition names = " + str(function.partition_names)))
                if hasattr(function, "value_labels"):
                    value_labels = function.value_labels()
                    if value_labels is not None:
                        labels = [x[1] for x in value_labels]
                        print(fill("  Value labels = " + str(labels)))
                print(
                    fill(
                        "  Requires event state = " + str(function.requires_event_state)
                    )
                )
                if hasattr(function, "is_log"):
                    if function.is_log:
                        print(fill("  Is log = " + str(function.is_log)))
                if hasattr(function, "initial_begin"):
                    print(
                        fill("  Default initial bin = " + str(function.initial_begin))
                    )
                    print(fill("  Default bin width = " + str(function.bin_width)))
                print(fill("  Default max size = " + str(function.max_size)))
                if hasattr(function, "tol"):
                    print(fill("  Default tol = " + str(function.tol)))
                print()

        functions = all_functions.generic_functions
        if len(functions):
            print("Selected event functions:\n")
            print_generic_functions(functions)

        int_functions = all_functions.discrete_vector_int_functions
        float_functions = all_functions.discrete_vector_float_functions
        continuous_1d_functions = all_functions.continuous_1d_functions

        if len(int_functions) + len(float_functions) + len(continuous_1d_functions):
            print("Selected event data functions:\n")
            print_functions(int_functions)
            print_functions(float_functions)
            print_functions(continuous_1d_functions)
