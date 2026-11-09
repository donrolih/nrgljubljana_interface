// Calls into the NRG Ljubljana libraries the interface links: the workflow API and numerics.
// No NRG run; checks that the libraries link and compute.

#include <gtest/gtest.h>

#include <nrgljubljana/numerics.hpp>
#include <nrgljubljana/solver_preparation.hpp>
#include <nrgljubljana/solver_workflow.hpp>

#include <cmath>
#include <map>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

static_assert(NRG::Workflow::workflow_api_version == 2, "nrgljubljana_interface requires NRG Ljubljana workflow API 2");

using NRG::Numerics::Table;

// A box spectrum S = -1 on [-D, D], on a symmetric mesh without w = 0
static Table box(double D, int n) {
  Table t;
  for (int i = n; i >= 1; --i) t.x.push_back(-D * i / n);
  for (int i = 1; i <= n; ++i) t.x.push_back(D * i / n);
  t.y.assign(t.x.size(), -1.0);
  return t;
}

TEST(NrgNumerics, Integrate) {
  auto t = box(2.0, 50);
  // The mesh starts and ends at -D and D, so the trapezoidal integral of -1 is exact
  EXPECT_NEAR(NRG::Numerics::integrate(t), -4.0, 1e-12);
}

TEST(NrgNumerics, KramersKronigBox) {
  double D = 2.0;
  auto t   = box(D, 50);
  auto k   = NRG::Numerics::kramers_kronig(t);
  ASSERT_EQ(k.x, t.x);
  // K(w) = PV int S(E)/(E-w) dE / pi = -log((D-w)/(D+w)) / pi inside the band
  for (std::size_t i = 1; i + 1 < k.x.size(); ++i) {
    double w = k.x[i];
    EXPECT_NEAR(k.y[i], -std::log((D - w) / (D + w)) / std::numbers::pi, 1e-10) << "w = " << w;
  }
}

TEST(NrgWorkflow, ResolveChainLength) {
  // Explicit Nmax is used as given
  EXPECT_EQ(NRG::Workflow::resolve_chain_length({{"Nmax", "7"}}), 7u);
  // Lambda = 2, Tmin = 1e-5 with the defaults z = 1, bandrescale = 1, discretization Z:
  // energy scales (1 - 1/Lambda)/log(Lambda) * Lambda^(-n/2) >= Tmin for n = 0..32
  EXPECT_EQ(NRG::Workflow::resolve_chain_length({{"Lambda", "2"}, {"Tmin", "1e-5"}}), 33u);
}

TEST(NrgWorkflow, ValidatePreparationParameters) {
  EXPECT_NO_THROW(NRG::Workflow::validate_preparation_parameters({{"Lambda", "2"}, {"Nmax", "33"}}));
  // Lambda must be larger than 1
  EXPECT_THROW(NRG::Workflow::validate_preparation_parameters({{"Lambda", "0.5"}, {"Nmax", "33"}}), std::runtime_error);
}
