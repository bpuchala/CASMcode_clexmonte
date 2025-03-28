import pytest

import libcasm.clexmonte as clexmonte
import libcasm.monte as monte
import libcasm.xtal as xtal


def test_run_series_1(Clex_ZrO_Occ_System, tmp_path):
    """A single run, using a fixture"""
    system = Clex_ZrO_Occ_System
    output_dir = tmp_path / "output"
    summary_file = output_dir / "summary.json"

    # construct a semi-grand canonical MonteCalculator
    mc_calculator = clexmonte.MonteCalculator(
        method="semigrand_canonical", system=system
    )

    # construct default sampling fixture parameters
    thermo = mc_calculator.make_default_sampling_fixture_params(
        label="thermo",
        output_dir=str(output_dir),
    )
    print(xtal.pretty_json(thermo.to_dict()))

    initial_conditions = {
        "temperature": 300.0,
        "param_chem_pot": [-4.0],
    }

    conditions_increment = {
        "temperature": 0.0,
        "param_chem_pot": [0.5],
    }

    n_states = 9

    # construct the initial state (default configuration)
    initial_state, motif, motif_id = clexmonte.make_initial_state(
        calculator=mc_calculator,
        conditions=initial_conditions,
        min_volume=1000,
    )

    output_params = clexmonte.RunDataOutputParams(
        do_save_all_initial_states=True,
        do_save_all_final_states=True,
        write_initial_states=True,
        write_final_states=True,
        output_dir=output_dir,
    )
    config_generator = clexmonte.FixedConfigGenerator(
        configuration=initial_state.configuration,
    )
    assert isinstance(config_generator, clexmonte.FixedConfigGenerator)

    state_generator = clexmonte.IncrementalConditionsStateGenerator(
        output_params=output_params,
        initial_conditions=initial_conditions,
        conditions_increment=conditions_increment,
        n_states=n_states,
        config_generator=config_generator,
        dependent_runs=True,
        modifiers=[],
    )

    engine = monte.RandomNumberEngine()
    global_cutoff = True
    # before_first_run = [] // TODO
    # before_each_run = [] // TODO
    sampling_fixture_params = [thermo]

    log = monte.MethodLog(
        logfile_path=str(output_dir / "log.txt"),
    )
    # log.reset_to_stdout()
    log.section("Begin: Monte Carlo calculation series")

    run_manager = clexmonte.RunManager(
        engine=engine,
        sampling_fixture_params=sampling_fixture_params,
        global_cutoff=global_cutoff,
    )

    log.print("Checking for completed runs...\n")
    state_generator.read_completed_runs()
    log.print(f"Found {state_generator.n_completed_runs}\n\n")

    while not state_generator.is_complete:
        run_manager.run_index = state_generator.n_completed_runs + 1

        log.print("Generating next state...\n")
        state = state_generator.next_state
        log.print(xtal.pretty_json(state.conditions.to_dict()))
        log.print("Done\n")

        run_data = clexmonte.RunData(
            initial_state=state,
        )

        log.print(f"Performing Run {run_manager.run_index}...\n")
        mc_calculator.run(
            state=state,
            run_manager=run_manager,
        )
        log.print(f"Run {run_manager.run_index} Done\n\n")

        run_data.final_state = state
        state_generator.append(run_data)
        state_generator.write_completed_runs()
    log.print("Monte Carlo calculation series complete\n")

    assert (output_dir / "completed_runs.json").exists()
    pytest.helpers.validate_summary_file(
        summary_file=summary_file, expected_size=n_states
    )
