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
class events_CompleteEventCalculator_Test
    : public test::KMCCompleteEventListTestSystem {};

/// \brief Test constructing event lists and calculating initial event states
///
/// Notes:
/// - FCC A-B-Va, 1NN interactions, A-Va and B-Va hops
/// - 10 x 10 x 10 (of the conventional 4-atom cell)
/// - expected runtime ~8s
TEST_F(events_CompleteEventCalculator_Test, Test1) {
  using namespace clexmonte;
  // --- State setup ---
  setup_input_files(false /*use_sparse_format_eci*/);

  // Create default state
  Index dim = 10;
  Eigen::Matrix3l T = test::fcc_conventional_transf_mat() * dim;
  monte::State<clexmonte::Configuration> state(
      make_default_configuration(*system, T));

  // Set configuration - A, with single Va
  Eigen::VectorXi &occupation = get_occupation(state);
  occupation(0) = 2;

  // Set conditions
  state.conditions.scalar_values.emplace("temperature", 600.0);

  /// --- KMC implementation ---

  make_prim_event_list();

  // Note: For correct atom tracking and stochastic canonical / semi-grand
  // canoncical
  //  event choosing, after this, occ_location->initialize must be called again
  // if the configuration is modified directly instead of via
  // occ_location->apply. Event calculations would be still be correct.
  make_complete_event_list(state);

  /// Make std::shared_ptr<clexmonte::Conditions> object from state.conditions
  auto conditions = make_conditions(*system, state);

  std::vector<kinetic::EventStateCalculator> prim_event_calculators =
      clexmonte::kinetic::make_prim_event_calculators(
          system, state, prim_event_list, conditions);
  EXPECT_EQ(prim_event_calculators.size(), 24);

  // Construct CompleteEventListEventCalculator
  auto event_calculator =
      std::make_shared<clexmonte::kinetic::CompleteEventCalculator>(
          prim_event_list, prim_event_calculators, event_list.events);

  double expected_Ef_kra = 1.0;
  double expected_freq = 1e12;
  double expected_rate =
      expected_freq * exp(-conditions->beta * expected_Ef_kra);

  Index i = 0;
  Index n_not_allowed = 0;
  Index n_allowed = 0;
  for (auto const &event : event_list.events) {
    auto const &event_id = event.first;
    // auto const &event_data = event.second;
    // auto const &prim_event_data = prim_event_list[event_id.prim_event_index];
    double rate = event_calculator->calculate_rate(event_id);
    auto const &event_state = event_calculator->event_state;

    if (CASM::almost_equal(rate, 0.0)) {
      // std::cout << "--- " << i << " ---" << std::endl;
      // print(std::cout, event_state, event_data, prim_event_data);
      // std::cout << std::endl;
      ++n_not_allowed;
    } else {
      EXPECT_EQ(rate, expected_rate);

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
  EXPECT_EQ(n_not_allowed, occupation.size() * 24 - 12);
  EXPECT_EQ(n_allowed, 12);
}
