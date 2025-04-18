import contextlib
import io

import numpy as np
import pytest

import libcasm.clexmonte as clexmonte
from libcasm.clexmonte.misc.contexts import (
    working_dir,
)
from libcasm.local_configuration import (
    LocalConfiguration,
    LocalConfigurationList,
)

CalculatorTestRunner = pytest.helpers.CalculatorTestRunner


def setup_1(runner: CalculatorTestRunner):
    """Adds the following:

    - A custom event state calculation function for the A_Va_1NN and B_Va_1NN events,
      which in turn just prints a message, calls the default event state calculation
      function, and increments a counter.
    - A custom selected event function, which just prints a message and increments a
      counter.

    runner.kinetics_params: SamplingFixtureParams
        Sampling fixture parameters:
        - sample_by_step()
        - set_max_count(100)

    runner.n_event_calculations: int
        The count of the number of calls of the custom event state calculation.

    runner.n_selected_events: int
        The count of the number of calls of the custom selected event function.

    Parameters
    ----------
    runner: CalculatorTestRunner
        The runner object with the MonteCalculator.

    Returns
    -------
    runner: CalculatorTestRunner
        The runner object with the added attributes.
    """

    kinetics_params = runner.calculator.make_default_sampling_fixture_params(
        label="kinetics",
        output_dir=str(runner.output_dir / "kinetics"),
        write_observations=False,
    )
    kinetics_params.sample_by_pass()
    kinetics_params.set_max_count(100)
    runner.kinetics_params = kinetics_params

    calculator = runner.calculator
    # calculator.collect("selected_event.by_equivalent_index")
    # calculator.collect("dE_activated.by_type", bin_width=0.01)

    # A custom event state calculation function which just increments a counter, prints
    # a message, and calls the default event state calculation
    runner.n_event_calculations = 0
    runner.n_encountered_abnormal = 0

    def event_state_calculation_f(
        state: clexmonte.EventState,  # <-- this is a mutable reference
        event_state_calculator: clexmonte.EventStateCalculator,
    ):
        runner.n_event_calculations += 1

        # print(f"CUSTOM EVENT STATE CALCULATION {runner.n_event_calculations}")
        #
        # prim_event_data = event_state_calculator.curr_prim_event_data
        # print(
        #     "event_type_name:",
        #     prim_event_data.event_type_name,
        #     "unitcell_index:",
        #     event_state_calculator.curr_unitcell_index,
        #     "equivalent_index:",
        #     prim_event_data.equivalent_index,
        # )
        event_state_calculator.set_default_event_state(state)
        if not state.is_normal:
            runner.n_encountered_abnormal += 1
            # print(f"- NOT NORMAL {runner.n_encountered_abnormal}")
        # print(
        #     "dE_activated:",
        #     state.dE_activated,
        #     "dE_final:",
        #     state.dE_final,
        #     "is_normal:",
        #     state.is_normal,
        # )
        # print(
        #     "rate:",
        #     state.rate,
        # )
        # print()
        # sys.stdout.flush()

    calculator.event_data.set_custom_event_state_calculation(
        event_type_name="A_Va_1NN",
        function=event_state_calculation_f,
    )

    calculator.event_data.set_custom_event_state_calculation(
        event_type_name="B_Va_1NN",
        function=event_state_calculation_f,
    )

    # A custom selected event function which just increments a counter and prints
    # a message
    runner.n_selected_events = 0

    def print_selected_event_f():
        runner.n_selected_events += 1
        # print(f"Event was selected {runner.n_selected_events}")

    calculator.add_generic_function(
        name="print_selected_event_f",
        description="Print information about the selected event state",
        requires_event_state=False,
        function=print_selected_event_f,
        order=0,
    )
    calculator.evaluate("print_selected_event_f")

    return runner


def add_initial_state(
    runner: CalculatorTestRunner,
    temperature: float,
    mol_composition: np.ndarray,
    min_volume: int = 1000,
):
    """Adds the following:

    runner.initial_state: MonteCarloState
        The initial state is set to a canonical state with
        temperature=1200.0 and mol_composition=[0.899, 0.1, 0.001],
        and volume=1000.

    runner.state: MonteCarloState
        A copy of `runner.initial_state`.

    runner.motif: Configuration
        The initial state's motif (default configuration)

    Parameters
    ----------
    runner: CalculatorTestRunner
        The runner object with the MonteCalculator.

    Returns
    -------
    runner: CalculatorTestRunner
        The runner object with the added attributes.
    """
    calculator = runner.calculator
    initial_state, motif = clexmonte.make_canonical_initial_state(
        calculator=calculator,
        conditions={
            "temperature": temperature,
            # one of param/mol composition is needed
            # "param_composition": [0.0, 0.0],
            "mol_composition": mol_composition,
        },
        min_volume=min_volume,
    )
    runner.initial_state = initial_state
    runner.state = initial_state.copy()
    runner.motif = motif
    return runner


def summarize(local_configurations, event_data):
    for event_type_name in local_configurations:
        _event_data = event_data[event_type_name]
        print(f"Event type name: {event_type_name}")
        if len(_event_data) == 0:
            print("- No abnormal events")
        else:
            for i, entry in enumerate(_event_data):
                s = entry.get("event_state")
                print(
                    f"- {i}: ",
                    f"Ef_kra={s['Ef_kra']:.6f}",
                    f"dE_final={s['dE_final']:.6f}",
                    f"dE_activated= {s['dE_activated']:.6f}",
                )
        print()


def assert_no_encountered_warnings(output: str):
    """Asserts that no warnings were printed to stdout."""
    msg = "## WARNING: ENCOUNTERED ABNORMAL EVENT ##############"
    assert msg not in output


def assert_has_encountered_warnings(output: str):
    """Asserts that no warnings were printed to stdout."""
    msg = "## WARNING: ENCOUNTERED ABNORMAL EVENT ##############"
    assert msg in output


def assert_no_selected_warnings(output: str):
    """Asserts that no warnings were printed to stdout."""
    msg = "## WARNING: SELECTED ABNORMAL EVENT #################"
    assert msg not in output


def assert_has_selected_warnings(output: str):
    """Asserts that no warnings were printed to stdout."""
    msg = "## WARNING: SELECTED ABNORMAL EVENT #################"
    assert msg in output


def is_encountered_exception(msg: str):
    expected_msg = "Error: encountered abnormal event, which is not allowed."
    return expected_msg in msg


def is_selected_exception(msg: str):
    expected_msg = "Error: selected abnormal event, which is not allowed."
    return expected_msg in msg


def run_test(
    runner: CalculatorTestRunner,
    temperature: float,
    mol_composition: list,
    seed: int,
    warn_if_encountered_event_is_abnormal: bool,
    throw_if_encountered_event_is_abnormal: bool,
    disallow_if_encountered_event_is_abnormal: bool,
    n_write_if_encountered_event_is_abnormal: int,
    warn_if_selected_event_is_abnormal: bool,
    throw_if_selected_event_is_abnormal: bool,
    n_write_if_selected_event_is_abnormal: int,
    max_n_attempts: int = 10,
):
    """Runs the test once with the given parameters and checks if behavior is consistent
    with options (assumes there is no guarantee that abnormal events are encountered
    or selected)."""
    runner = setup_1(runner)
    runner.calculator.engine.seed(seed)
    runner = add_initial_state(
        runner=runner,
        temperature=temperature,
        mol_composition=mol_composition,
    )
    runner.calculator.set_state_and_potential(runner.state)
    runner.calculator.make_occ_location()

    print(f"Test output dir: {runner.output_dir}\n")

    n_encountered_exception = 0
    n_selected_exception = 0
    n_write_encountered = 0
    n_write_selected = 0

    n_attempts = 0
    while True:
        threw_encountered_exception = False
        threw_selected_exception = False

        # Run - expect an exception due to an abnormal event
        f = io.StringIO()
        with contextlib.redirect_stdout(f), contextlib.redirect_stderr(f):
            try:
                runner.calculator.run_fixture(
                    state=runner.state,
                    sampling_fixture_params=runner.kinetics_params,
                )

            except Exception as e:
                if is_encountered_exception(str(e)):
                    threw_encountered_exception = True
                    n_encountered_exception += 1
                elif is_selected_exception(str(e)):
                    threw_selected_exception = True
                    n_selected_exception += 1
        output = f.getvalue()
        # print("!!!!")
        # print(output)
        # print("!!!!")

        event_data = runner.calculator.event_data
        n_encountered_abnormal = event_data.n_encountered_abnormal
        n_encountered_abnormal_sum = sum(n_encountered_abnormal.values())
        n_selected_abnormal = event_data.n_selected_abnormal
        n_selected_abnormal_sum = sum(n_selected_abnormal.values())

        # -- Check warnings --
        if not warn_if_encountered_event_is_abnormal:
            assert_no_encountered_warnings(output)
        elif n_encountered_abnormal_sum > 0:
            assert_has_encountered_warnings(output)
        if not warn_if_selected_event_is_abnormal:
            assert_no_selected_warnings(output)
        elif n_selected_abnormal_sum > 0:
            assert_has_selected_warnings(output)

        # -- Check throws --
        if not throw_if_encountered_event_is_abnormal:
            assert not threw_encountered_exception
        elif n_encountered_abnormal_sum > 0:
            assert threw_encountered_exception
        if not throw_if_selected_event_is_abnormal:
            assert not threw_selected_exception
        elif n_selected_abnormal_sum > 0:
            assert threw_selected_exception

        # -- Check disallows --
        if disallow_if_encountered_event_is_abnormal:
            assert n_selected_abnormal_sum == 0
            assert not threw_selected_exception

        # --- Check for abnormal_events jsonl files ---
        # - Note: these files should be appended as multiple runs are done, but
        #   always only include LocalConfiguration with unique correlations
        file = runner.output_dir / "encountered_abnormal_events.jsonl"
        if n_write_if_encountered_event_is_abnormal > 0:
            if file.exists():
                lines = 0
                with open(file, "r") as f:
                    for line in f:
                        lines += 1
                n_write_encountered = lines

                assert lines <= n_write_if_encountered_event_is_abnormal
        else:
            assert not file.exists()

        file = runner.output_dir / "selected_abnormal_events.jsonl"
        if n_write_if_selected_event_is_abnormal > 0:
            if file.exists():
                with open(file, "r") as f:
                    lines = 0
                    for line in f:
                        lines += 1
                n_write_selected = lines

                assert lines <= n_write_if_selected_event_is_abnormal
        else:
            assert not file.exists()

        # --- Test the read_abnormal_events function ---

        def _checks(n, local_configurations, event_data, n_expected, which):
            if n_expected == 0:
                assert n is None
                return

            assert n == n_expected
            for event_type_name in local_configurations:
                assert event_type_name in event_data
                _local_config_list = local_configurations[event_type_name]
                assert isinstance(_local_config_list, LocalConfigurationList)
                _event_data = event_data[event_type_name]
                assert len(_local_config_list) == len(_event_data)

                # Check types, dict contents, and check for unique local corr
                # (but same local config may exist in the "encountered" and "selected"
                # lists" so skip which=="any")
                local_corr_list = list()
                for lc, s in zip(_local_config_list, _event_data):
                    assert isinstance(lc, LocalConfiguration)
                    assert isinstance(s, dict)
                    assert "event_state" in s
                    assert "prim_event_data" in s

                    if which == "any":
                        continue
                    local_corr_i = np.array(s["event_state"]["local_corr"])
                    for local_corr_j in local_corr_list:
                        assert not np.allclose(local_corr_i, local_corr_j, atol=1e-5)

                    local_corr_list.append(local_corr_i)

        # which=="all"
        calculator = runner.calculator
        n, local_configurations, event_data = calculator.read_abnormal_events()
        n_expected = n_write_encountered + n_write_selected
        _checks(n, local_configurations, event_data, n_expected, which="any")

        # which=="encountered"
        n, local_configurations, event_data = calculator.read_abnormal_events(
            which="encountered",
        )
        n_expected = n_write_encountered
        _checks(n, local_configurations, event_data, n_expected, which="encountered")

        # which=="selected"
        n, local_configurations, event_data = calculator.read_abnormal_events(
            which="selected",
        )
        n_expected = n_write_selected
        _checks(n, local_configurations, event_data, n_expected, which="selected")

        n_attempts += 1

        if n_attempts == max_n_attempts:
            break


def test_abnormal_event_handling_1a(FCCBinaryVacancy_kmc_System_2, tmp_path):
    """Test the default handling of abnormal events."""
    with working_dir(tmp_path):
        run_test(
            runner=CalculatorTestRunner(
                system=FCCBinaryVacancy_kmc_System_2,
                method="kinetic",
                params={
                    # "verbosity": "quiet",
                    # "abnormal_event_handling": None,
                },
                output_dir=tmp_path / "output",
            ),
            temperature=1200.0,
            mol_composition=[0.799, 0.2, 0.001],
            seed=0,
            warn_if_encountered_event_is_abnormal=True,
            throw_if_encountered_event_is_abnormal=True,
            disallow_if_encountered_event_is_abnormal=False,
            n_write_if_encountered_event_is_abnormal=100,
            warn_if_selected_event_is_abnormal=True,
            throw_if_selected_event_is_abnormal=True,
            n_write_if_selected_event_is_abnormal=100,
            max_n_attempts=10,
        )


def test_abnormal_event_handling_1b(FCCBinaryVacancy_kmc_System_2, tmp_path):
    """Test only throw on selected abnormal events

    Seed = 699, mol_composition=[0.799, 0.2, 0.001],
    for an initial state that has an abnormal event.
    """

    with working_dir(tmp_path):
        run_test(
            runner=CalculatorTestRunner(
                system=FCCBinaryVacancy_kmc_System_2,
                method="kinetic",
                params={
                    # "verbosity": "quiet",
                    "abnormal_event_handling": {
                        "encountered_events": {
                            "throw": False,
                        },
                    }
                },
                output_dir=tmp_path / "output",
            ),
            temperature=1200.0,
            mol_composition=[0.799, 0.2, 0.001],
            seed=0,
            warn_if_encountered_event_is_abnormal=True,
            throw_if_encountered_event_is_abnormal=False,
            disallow_if_encountered_event_is_abnormal=False,
            n_write_if_encountered_event_is_abnormal=100,
            warn_if_selected_event_is_abnormal=True,
            throw_if_selected_event_is_abnormal=True,
            n_write_if_selected_event_is_abnormal=100,
            max_n_attempts=10,
        )


def test_abnormal_event_handling_1c(FCCBinaryVacancy_kmc_System_2, tmp_path):
    """Test only throw on selected abnormal events

    with lower temperature - more likely to encounter but never select an abnormal event


    Seed = 699, mol_composition=[0.799, 0.2, 0.001],
    for an initial state that has an abnormal event.
    """
    with working_dir(tmp_path):
        # Test only throw on selected abnormal events:
        run_test(
            runner=CalculatorTestRunner(
                system=FCCBinaryVacancy_kmc_System_2,
                method="kinetic",
                params={
                    # "verbosity": "debug",
                    "abnormal_event_handling": {
                        "encountered_events": {
                            "throw": False,
                        },
                    },
                },
                output_dir=tmp_path / "output",
            ),
            temperature=300.0,
            mol_composition=[0.799, 0.2, 0.001],
            seed=0,
            warn_if_encountered_event_is_abnormal=True,
            throw_if_encountered_event_is_abnormal=False,
            disallow_if_encountered_event_is_abnormal=False,
            n_write_if_encountered_event_is_abnormal=100,
            warn_if_selected_event_is_abnormal=True,
            throw_if_selected_event_is_abnormal=True,
            n_write_if_selected_event_is_abnormal=100,
            max_n_attempts=10,
        )


def test_abnormal_event_handling_1d(FCCBinaryVacancy_kmc_System_2, tmp_path):
    """Test no throw on abnormal events

    Seed = 699, mol_composition=[0.799, 0.2, 0.001],
    for an initial state that has an abnormal event.
    """

    with working_dir(tmp_path):
        run_test(
            runner=CalculatorTestRunner(
                system=FCCBinaryVacancy_kmc_System_2,
                method="kinetic",
                params={
                    # "verbosity": "quiet",
                    "abnormal_event_handling": {
                        "encountered_events": {
                            "throw": False,
                        },
                        "selected_events": {
                            "throw": False,
                        },
                    }
                },
                output_dir=tmp_path / "output",
            ),
            temperature=1200.0,
            mol_composition=[0.799, 0.2, 0.001],
            seed=0,
            warn_if_encountered_event_is_abnormal=True,
            throw_if_encountered_event_is_abnormal=False,
            disallow_if_encountered_event_is_abnormal=False,
            n_write_if_encountered_event_is_abnormal=100,
            warn_if_selected_event_is_abnormal=True,
            throw_if_selected_event_is_abnormal=False,
            n_write_if_selected_event_is_abnormal=100,
            max_n_attempts=10,
        )


def test_abnormal_event_handling_1e(FCCBinaryVacancy_kmc_System_2, tmp_path):
    """Test no throw on abnormal events

    Seed = 699, mol_composition=[0.799, 0.2, 0.001],
    for an initial state that has an abnormal event.
    """

    with working_dir(tmp_path):
        run_test(
            runner=CalculatorTestRunner(
                system=FCCBinaryVacancy_kmc_System_2,
                method="kinetic",
                params={
                    # "verbosity": "quiet",
                    "abnormal_event_handling": {
                        "encountered_events": {
                            "throw": False,
                            "disallow": True,
                        },
                        "selected_events": {
                            "throw": True,
                        },
                    }
                },
                output_dir=tmp_path / "output",
            ),
            temperature=1200.0,
            mol_composition=[0.799, 0.2, 0.001],
            seed=0,
            warn_if_encountered_event_is_abnormal=True,
            throw_if_encountered_event_is_abnormal=False,
            disallow_if_encountered_event_is_abnormal=True,
            n_write_if_encountered_event_is_abnormal=100,
            warn_if_selected_event_is_abnormal=True,
            throw_if_selected_event_is_abnormal=True,
            n_write_if_selected_event_is_abnormal=100,
            max_n_attempts=10,
        )
