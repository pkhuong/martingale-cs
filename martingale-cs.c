#include "martingale-cs.h"

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <string.h>

/* Pairwise <= test is the base case. */
const double martingale_cs_le = 0;

/* -log 2 rounded away from 0. */
const double martingale_cs_eq = -0.6931471805599454;

/* -1/2 log log 2, rounded up. */
static const double minus_half_log_log_2_up = 0.1832564602908322;

/*
 * Safe-rounding utilities.  Not because it makes a difference, but
 * because extreme p-values means we should be extra confidence.
 */
static inline uint64_t float_bits(double x)
{
	uint64_t bits;
	uint64_t mask;

	memcpy(&bits, &x, sizeof(bits));
	/* extract the sign bit. */
	mask = (int64_t)bits >> 63;
	/*
	 * If negative, flip the significand bits to convert from
	 * sign-magnitude to 2's complement.
	 */
	return bits ^ (mask >> 1);
}

static inline double bits_float(uint64_t bits)
{
	double ret;
	uint64_t mask;

	mask = (int64_t)bits >> 63;
	/* Undo the bit-flipping above. */
	bits ^= (mask >> 1);
	memcpy(&ret, &bits, sizeof(ret));
	return ret;
}

static inline double next_k(double x, uint64_t delta)
{
	return bits_float(float_bits(x) + delta);
}

static inline double next(double x) { return next_k(x, 1); }

static inline double prev_k(double x, uint64_t delta)
{
	return bits_float(float_bits(x) - delta);
}

static inline double prev(double x) { return prev_k(x, 1); }

/* Assume libm is off by < 4 ULPs. */
static const uint64_t libm_error_limit = 4;

static inline double log_up(double x)
{
	return next_k(log(x), libm_error_limit);
}

static inline double log2_down(double x)
{
	return prev_k(log2(x), libm_error_limit);
}

static inline double sqrt_up(double x)
{
	/* sqrt is supposed to be rounded correctly. */
	return next(sqrt(x));
}

int martingale_cs_check_constants(void)
{
	int ret = 0;

	/*
	 * Use memcpy instead of float_bits to directly compare bit
	 * patterns in sign-magnitude instead of float_bits's
	 * conversion to 2's complement.
	 */

	size_t index = 0;
#define CHECK(NAME, EXPECTED)                                                \
	do {                                                                 \
		assert(index < CHAR_BIT * sizeof(int) - 1);                  \
		uint64_t actual;                                             \
		memcpy(&actual, &NAME, sizeof(actual));                      \
		if (actual != (uint64_t)EXPECTED) {                          \
			ret |= 1 << index;                                   \
		}                                                            \
		++index;                                                     \
	} while (0)

	/* le is the default: no adjustment. */
	CHECK(martingale_cs_le, 0);
	CHECK(martingale_cs_eq, -4618953502541334032ULL);
	CHECK(minus_half_log_log_2_up, 4595770530100767648LL);
#undef CHECK

	return ret;
}

/* We let C and alpha = 2, like Darling and Robbins. */
static const double c = 2;

/*
 * Returns the log(A) term, the main factor in how far away
 * we expect the martingale to stray from the mean of 0.
 *
 * Scales linearly with log(eps), and inversely with log log
 * min_count.
 *
 *  Q_m = 1 / [lg_2 m - 1/2],
 * and
 *  Q_m / A <= eps
 *
 * <-> log(Q_m) - log(A) <= log(eps)
 * <-> log(A) >= log(Q_m) - log(eps)
 */
static const double log_a_up(uint64_t min_count, double log_eps)
{
	/*
	 * 1 / Q_m = lg_2 m - 1/2, and we want to round down, in
	 * order to round Q_m up, and thus also over-approximate
	 * log(A).
	 *
	 * We assume the conversion of `min_count` to double is
	 * exact. If it isn't, the values are so large that
	 * rounding hopefully doesn't matter.
	 */
	const double inv_q_m = prev(log2_down(min_count) - 0.5);

	return log_up(next(1.0 / inv_q_m)) - log_eps;
}

double martingale_cs_threshold(uint64_t n, uint64_t min_count, double log_eps)
{
	assert(log_eps <= 0 && "Positive log_eps means > 100% false positive "
			       "rate. Should it be negated?");

	if (min_count < c) {
		min_count = c;
	}

	if (n < min_count) {
		return HUGE_VAL;
	}

	if (log_eps >= 0) {
		/* >= 100% false positive rate: just always reject. */
		return -HUGE_VAL;
	}

	const double log_a = log_a_up(min_count, log_eps);

	/*
	 * n f_n(A)
	 *   = sqrt(n) (3 / 2sqrt(2)) sqrt(4 log log n - 4 log log2 + 2 log A)
	 *   = 3 sqrt(n) sqrt[(4 log log n - 4 log log 2 + 2 log A) / 8]
	 *   = 3 sqrt[n (1/2 log log n - 1/2 log log 2 + 1/4 log A)].
	 */

	const double inner
	    = next(next(.5 * log_up(log_up(n)) + minus_half_log_log_2_up)
		+ 0.25 * log_a);
	return next(3 * sqrt_up(next(n * inner)));
}

/*
 * Hoeffding's lemma guarantees that any zero-mean distribution with a
 * range of span 2 satisfies our constraint that `mgf <= exp(t^2 /
 * 2)`.
 *
 * Rescale the width returned by `martingale_cs_threshold` as if the
 * `width = 2`.
 */
double martingale_cs_threshold_span(
    uint64_t n, uint64_t min_count, double span, double log_eps)
{
	const double scale = span / 2; /* Division by 2 is exact. */
	return next(scale * martingale_cs_threshold(n, min_count, log_eps));
}

/*
 * The classic proof of Hoeffding's lemma eventually gets to a point
 * where we upper bound the expression
 *   t (1 - t),
 * where
 *   t = (rho exp(v)) / (1 - rho + rho exp(v)),
 *   rho = -lo / (hi - lo)
 *   v an arbitrary real >= 0.
 *
 * t is a monotonically decreasing function of v, and
 *   lim_{v -> \infty} t = 1 from below.
 *
 * t (1 - t) is maximised at t = 0.5.  When rho <= 0.5, the mean value
 * theorem tells us there exists a v such that t = 0.5 may be
 * achieved, and there's nothing to gain compared to
 * martingale_cs_threshold_span.  However, when rho > 0.5, the
 * expression is maximised at v = 0.  In that case, we may widen the
 * span that satisfies Darling and Robbins's condition on the mgf
 *
 * We have that mgf <= exp[1/2 (rho (1 - rho)) (hi - lo)^2 \lambda^t],
 * and thus only need
 *   hi - lo <= 1/sqrt[rho (1 - rho)]
 * to guarantee mgf(\lambda) <= exp(1/2 \lambda^2).
 */
double martingale_cs_threshold_range(
    uint64_t n, uint64_t min_count, double lo, double hi, double log_eps)
{
	/*
	 * With this kind of range, the random values must all be exactly 0
	 * to achieve a mean of zero.
	 */
	if (lo >= 0 || hi <= 0) {
		return 0;
	}

	const double span = next(hi - lo);
	const double rho = prev(-lo / span);
	double scale;
	if (rho <= 0.5) {
		scale = span / 2;
	} else {
		/*
		 * Ideal span is 1 / sqrt[rho (1 - rho)], so we must scale
		 * `span` by `span / ideal_span = sqrt[rho (1 - rho)] * span`
		 */
		scale = next(sqrt_up(rho * next(1 - rho)) * span);
	}

	return next(scale * martingale_cs_threshold(n, min_count, log_eps));
}

double martingale_cs_quantile_slop(
    double quantile, uint64_t n, uint64_t min_count, double log_eps)
{
	assert(quantile >= 0 && quantile <= 1.0
	    && "Quantile is a fraction in [0, 1]. Was a percentile passed in "
	       "without dividing by 100?");

	if (quantile <= 0.0 || quantile >= 1.0) {
		return 1;
	}

	/*
	 * Extend the range from Darling and Robbins to account for
	 * equality.
	 *
	 * We can't use f(x) = -0.5 if x <= median else 0.5 (for
	 * example): since `x = median` happens with non-zero
	 * probability, we must add a third case:
	 *
	 *  f(x) = -0.5 if x < median
	 *       |  0.0 if x = median
	 *       |  0.5 if x > median
	 *
	 * We must thus extend the range by one more observation,
	 * since that last observation might have incurred 0 "cost" in
	 * the martingale.
	 *
	 * Letting the range be 1 means that each unexpected value
	 * over or under the quantile "costs" 1 in terms of distance
	 * from the expected sum value, which is exactly what we want
	 * for this quantile slop.
	 */
	return 1 + martingale_cs_threshold_span(
		       n, min_count, 1.0, log_eps + martingale_cs_eq);
}

double martingale_cs_quantile_slop_hi(
    double quantile, uint64_t n, uint64_t min_count, double log_eps)
{
	assert(quantile >= 0 && quantile <= 1.0
	    && "Quantile is a fraction in [0, 1]. Was a percentile passed in "
	       "without dividing by 100?");

	if (quantile <= 0.0) {
		return 1;
	}

	if (quantile >= 1.0) {
		return HUGE_VAL;
	}

	/*
	 * If, e.g. quantile = 0.9, then we pay -0.1 for x < quantile,
	 * and .9 for x > quantile.
	 */
	return 1 + martingale_cs_threshold_range(n, min_count, quantile - 1,
		       quantile, log_eps + martingale_cs_eq);
}

double martingale_cs_quantile_slop_lo(
    double quantile, uint64_t n, uint64_t min_count, double log_eps)
{
	assert(quantile >= 0 && quantile <= 1.0
	    && "Quantile is a fraction in [0, 1]. Was a percentile passed in "
	       "without dividing by 100?");

	if (quantile <= 0.0) {
		return -HUGE_VAL;
	}

	if (quantile >= 1.0) {
		return -1;
	}

	return -1 - martingale_cs_threshold_range(n, min_count, -quantile,
			1 - quantile, log_eps + martingale_cs_eq);
}
