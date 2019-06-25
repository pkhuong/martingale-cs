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
	if (min_count < 2) {
		min_count = c;
	}

	if (n < min_count) {
		return HUGE_VAL;
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
