#include "KMCCompleteEventListTestSystem.hh"
#include "casm/clexmonte/kinetic/io/stream/EventState_stream_io.hh"
#include "casm/clexmonte/kinetic/kinetic_events.hh"
#include "casm/clexmonte/state/Conditions.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "gtest/gtest.h"
#include "teststructures.hh"

using namespace CASM;

/// NOTE:
/// - This test is designed to copy data to the same directory each time, so
///   that the Clexulators do not need to be re-compiled.
/// - To clear existing data, remove the directory:
//    CASM_test_projects/FCCBinaryVacancy_default directory
class events_EventStateCalculator_Test
    : public test::KMCCompleteEventListTestSystem {
 protected:
  void run_checks() {
    using namespace clexmonte;

    // --- State setup ---

    // Create default state
    Eigen::Matrix3l T = test::fcc_conventional_transf_mat() * 10;
    state_type state(make_default_configuration(*system, T));
    // Set configuration - A, with single Va
    Eigen::VectorXi &occupation = get_occupation(state);
    occupation(0) = 2;

    // Set conditions
    state.conditions.scalar_values.emplace("temperature", 600.0);

    /// --- KMC implementation ---

    make_prim_event_list();
    // std::cout << "#prim events: " << prim_event_list.size() << std::endl;

    // Note: This calls occ_location->initialize. For correct atom tracking and
    // stochastic canonical / semi-grand canoncical event choosing,
    // occ_location->initialize must be called again if the configuration is
    // modified directly instead of via occ_location->apply. Event calculations
    // would be still be correct.
    make_complete_event_list(state);
    // std::cout << "#events: " << event_list.events.size() << std::endl;

    /// Make std::shared_ptr<clexmonte::Conditions> object from state.conditions
    auto conditions = make_conditions(*system, state);

    std::vector<kinetic::EventStateCalculator> prim_event_calculators =
        kinetic::make_prim_event_calculators(system, state, prim_event_list,
                                             conditions);
    EXPECT_EQ(prim_event_calculators.size(), 24);
    // std::cout << "#prim event calculators: " << prim_event_calculators.size()
    //           << std::endl;

    double expected_Ef_kra = 1.0;
    double expected_freq = 1e12;
    double expected_rate =
        expected_freq * exp(-conditions->beta * expected_Ef_kra);

    Index i = 0;
    Index n_allowed = 0;
    kinetic::EventState event_state;
    for (auto const &event : event_list.events) {
      auto const &event_id = event.first;
      auto const &event_data = event.second;
      auto const &prim_event_data = prim_event_list[event_id.prim_event_index];
      auto const &prim_event_calculator =
          prim_event_calculators[event_id.prim_event_index];
      prim_event_calculator.calculate_event_state(event_state, event_data,
                                                  prim_event_data);

      if (event_state.is_allowed) {
        // std::cout << "--- " << i << " ---" << std::endl;
        // print(std::cout, event_state, event_data, prim_event_data);
        // std::cout << std::endl;
        EXPECT_TRUE(CASM::almost_equal(event_state.dE_final, 0.0));
        EXPECT_TRUE(
            CASM::almost_equal(event_state.dE_activated, expected_Ef_kra));
        EXPECT_TRUE(CASM::almost_equal(event_state.Ef_kra, expected_Ef_kra));
        EXPECT_TRUE(CASM::almost_equal(event_state.freq, expected_freq));
        EXPECT_TRUE(CASM::almost_equal(event_state.rate, expected_rate));
        ++n_allowed;
      }
      ++i;
    }
    EXPECT_EQ(n_allowed, 12);
  }
};

/// \brief Test constructing event lists and calculating initial event states
///
/// Notes:
/// - FCC A-B-Va, 1NN interactions, A-Va and B-Va hops
/// - 10 x 10 x 10 (of the conventional 4-atom cell)
/// - expected runtime ~5s
TEST_F(events_EventStateCalculator_Test, Test1) {
  setup_input_files(false /*use_sparse_format_eci*/);
  run_checks();
}

/// \brief Test constructing event lists and calculating initial event states
///
/// Notes:
/// - FCC A-B-Va, 1NN interactions, A-Va and B-Va hops
/// - 10 x 10 x 10 (of the conventional 4-atom cell)
/// - expected runtime ~5s
TEST_F(events_EventStateCalculator_Test, Test2) {
  setup_input_files(true /*use_sparse_format_eci*/);
  run_checks();
}
