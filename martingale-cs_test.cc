#include "martingale-cs.h"

#include <climits>
#include <cmath>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace {
using ::testing::DoubleNear;
using ::testing::Lt;

TEST(MartingaleCs, ConstantsOk)
{
	EXPECT_EQ(martingale_cs_check_constants(), 0);
}

// Darling and Robbins have an example with
// a = c = 2, m = 32, eps = 0.05.
//  -> A = 80/9,
//  -> 2 f_n(A) = 3 sqrt[n(log log n + 1.457) / 2]
TEST(MartingaleCs, Golden)
{
	for (size_t i = 32; i < 64; ++i) {
		const double expected = 3
		    * std::sqrt(0.5 * i * (std::log(std::log(i)) + 1.457));
		EXPECT_THAT(martingale_cs_threshold(
				i, 32, std::log(0.05) + martingale_cs_eq),
		    DoubleNear(expected, 1e-2));
	}
}

// n < min_count -> infinite threshold
TEST(MartingaleCs, TooEarly)
{
	EXPECT_EQ(martingale_cs_threshold(1, 10, -10), HUGE_VAL);
}

TEST(MartingaleCs, DefaultMinCount)
{
	EXPECT_EQ(martingale_cs_threshold(1, 1, -10), HUGE_VAL);

	EXPECT_EQ(martingale_cs_threshold(1000000, 1, -2),
	    martingale_cs_threshold(1000000, 2, -2));
}

// Higher n -> higher absolute threshold, lower relative to n.
TEST(MartingaleCs, MonotonicN)
{
	EXPECT_GT(martingale_cs_threshold(1001, 10, -10),
	    martingale_cs_threshold(1000, 10, -10));

	EXPECT_LT(martingale_cs_threshold(1001, 10, -10) / 1001,
	    martingale_cs_threshold(1000, 10, -10) / 1000);
}

// Higher min count -> lower threshold.
TEST(MartingaleCs, MonotonicMinCount)
{
	EXPECT_LT(martingale_cs_threshold(1000, 11, -10),
	    martingale_cs_threshold(1000, 10, -10));
}

// Lower eps -> higher threhsold.
TEST(MartingaleCs, MonotonicEps)
{
	EXPECT_GT(martingale_cs_threshold(1000, 10, -5),
	    martingale_cs_threshold(1000, 10, -4));
}

TEST(MartingaleCs, Quantile)
{
	// Compare against the expression given in the paper.
	EXPECT_EQ(martingale_cs_quantile_slop(0.5, 1000, 32, std::log(0.05)),
	    0.5 * martingale_cs_threshold(
		      1000, 32, std::log(0.05) + martingale_cs_eq));

	// And make sure we handle the symmetry correctly.
	EXPECT_EQ(martingale_cs_quantile_slop(0.1, 10000, 3, std::log(0.001)),
	    0.9 * martingale_cs_threshold(
		      10000, 3, std::log(0.001) + martingale_cs_eq));

	EXPECT_EQ(martingale_cs_quantile_slop(0.9, 10000, 3, std::log(0.01)),
	    0.9 * martingale_cs_threshold(
		      10000, 3, std::log(0.01) + martingale_cs_eq));
}
} // namespace
