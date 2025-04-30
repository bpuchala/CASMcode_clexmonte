import pytest

import libcasm.clexmonte as clexmonte

CalculatorTestRunner = pytest.helpers.CalculatorTestRunner


@pytest.fixture
def FCCBinaryVacancy_runner(FCCBinaryVacancy_kmc_System, tmp_path):
    runner = CalculatorTestRunner(
        system=FCCBinaryVacancy_kmc_System,
        method="kinetic",
        params={
            # "verbosity": "debug",
        },
        output_dir=tmp_path / "output",
    )

    return runner


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
        output_dir=str(runner.output_dir),
        write_observations=True,
    )
    kinetics_params.sample_by_step()
    kinetics_params.set_max_count(100)
    runner.kinetics_params = kinetics_params

    calculator = runner.calculator
    # calculator.collect("selected_event.by_equivalent_index")
    # calculator.collect("dE_activated.by_type", bin_width=0.01)

    # A custom event state calculation function which just increments a counter, prints
    # a message, and calls the default event state calculation
    runner.n_event_calculations = 0

    def event_state_calculation_f(
        state: clexmonte.EventState,  # <-- this is a mutable reference
        event_state_calculator: clexmonte.EventStateCalculator,
    ):
        runner.n_event_calculations += 1

        print(f"CUSTOM EVENT STATE CALCULATION {runner.n_event_calculations}")
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
        print(f"Event was selected {runner.n_selected_events}")

    calculator.add_generic_function(
        name="print_selected_event_f",
        description="Print information about the selected event state",
        requires_event_state=False,
        function=print_selected_event_f,
        order=0,
    )
    calculator.evaluate("print_selected_event_f")

    return runner


def initial_state_1(runner: CalculatorTestRunner):
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
            "temperature": 1200.0,
            # one of param/mol composition is needed
            # "param_composition": [0.0, 0.0],
            "mol_composition": [0.899, 0.1, 0.001],
        },
        min_volume=1000,
    )
    runner.initial_state = initial_state
    runner.state = initial_state.copy()
    runner.motif = motif
    return runner


def test_custom_event_state_calculation_f(FCCBinaryVacancy_runner):
    """Test using custom event state calculation functions

    Expected number of custom event state calculation calls: 1300

    - 12 possible 1NN Va hops in FCC
    - 1 selected event state calculation (for handling selected events without barriers)
    - 100 steps
    - 1 step calculated but not applied

    -> (12 + 1) * 100 + 12 = 1312 expected calls

    Expected number of selected events: 100

    """

    runner = setup_1(FCCBinaryVacancy_runner)
    runner = initial_state_1(runner)

    # Run
    kinetics = runner.calculator.run_fixture(
        state=runner.state,
        sampling_fixture_params=runner.kinetics_params,
    )
    assert isinstance(kinetics, clexmonte.SamplingFixture)

    # n_event_calculations is incremented for every allowed event rate update...
    # so at least 100
    assert runner.n_event_calculations == 1312
    assert runner.n_selected_events == 100
