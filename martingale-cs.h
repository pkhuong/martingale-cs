#ifndef MARTINGALE_CS_H
#define MARTINGALE_CS_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

extern const double martingale_cs_le;
extern const double martingale_cs_eq;

/*
 * Returns 0 if the constants were definitely compiled correctly,
 * non-zero otherwise.
 *
 * When non-zero, the return value is a bitmask with ones for each
 * constant with an incorrect value, in order:
 *
 * bit 0: cs_le
 * bit 1: cs_eq
 * bit 2: an internal constant
 */
int martingale_cs_check_constants(void);

/*
 * Confidence interval sequence for simple martingales.
 *
 * `n` is the number of values observed so far, `min_count` is the
 * minimum number of values at which we output useful thresholds (when
 * `n < min_count`, the return value is always +infty), and `log_eps`
 * is the natural log of the allowed false positive (type I) rate.
 *
 * `min_count` should be at least 2, abd `log_eps` should be strictly
 * negative.
 *
 * Let X be a random variable with 0 mean and momemt generating
 * function mfg such that `mgf(t) = E[exp(tx)] <= exp(t^2 / 2)` for
 * all `t >= 0`.
 *
 * This function returns the width of a `1 - exp(log_eps)`-confidence
 * interval for the sum of `n` i.i.d. values sampled from `X`. This
 * interval is so conservative that we can compare our running sum
 * against it after every observation of a new value from `X`, and
 * only risk a false positive with probability at lmost
 * `exp(log_eps)`, regardless of how many comparisons we make, assuming
 * that we start comparing at `min_count >= 2`.
 *
 * By default, this function implements a one-sided test, i.e., the
 * confidence interval (CI) for the sum is (-infty, threshold).  Add
 * `martingale_cs_eq` to the `log_eps` for the half-interval of a
 * two-sided test: the CI becomes (-threshold, threshold).
 *
 * This confidence sequence can be used to test any hypothesis about a
 * statistic that is a sum of a value derived from each observation
 * independently, where the null hypothesis is that the statistic has
 * expected value 0.
 *
 * See [Confidence sequences for mean, variance, and
 * median](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC335597/) by
 * Darling and Robbins (1967) for the math.
 *
 * For example, compare the mean of two bounded random variables Y and
 * Z by renormalizing their range so that the difference of two
 * variables falls in [-1, 1], and accumulate the difference between
 * i.i.d. observations from Y and X.  This sum should be concentrated
 * around 0.  If strays far away to exceed `martingale_cs_threshold`,
 * we can reject the null: the means of the two random variables most
 * likely differ.
 *
 * We can also use this interval to estimate quantiles.  Let q be the
 * unknown 90th percentile for some random variable W.  We can derive
 * another random variable from W: W_i = -(0.1/0.9) if X_i <= q, and 1
 * if X_i > q.  W has zero mean and remains in the [-1, 1] range.  The
 * 2-sided confidence interval returned by
 * `martingale_cs_threshold(..., + martingale_cs_eq)` bounds how far
 * away we can expect |Sum w_i| to be from 0.  Set delta equal to
 * `martingale_cs_threshold(..., log(eps) + martingale_cs_le)`.
 *
 * We know that we can expect Sum w_i to fall in [-delta, delta] with
 * probability 1 - eps.  At the bottom edge, we have
 *
 *   -0.1/0.9 #{X_i <= q} + #{X_i > q} = -delta
 *            #{X_i <= q} + #{X_i > q} = n
 *
 * -> (1 + 0.1/0.9) #{X_i <= q} = 10/9 #{X_i <= q } = n + delta
 *
 * and thus at most 9/10 (n + delta) observations should fall at or
 * below the 90th percentile q.  Conversely, we can thus expect that q
 * is at most the `ceil[9/10 (n + delta)]` smallest observations.  In
 * the other direction, we solve
 *
 *   -0.1/0.9 #{X_i <= q} + #{X_i > q} = delta
 *            #{X_i <= q} + #{X_i > q} = n
 *
 * -> (1 + 0.1/0.9) #{X_i <= q} = 10/9 #{X_i <= q} = n - delta
 *
 * and thus at least 9/10 (n - delta) observations should fall at or
 * below the 90th percentile q. Conversely, we can expect that q is at
 * least the `floor[9/10 (n - delta)]` smallest observation.
 *
 * This pair of lower and upper bounds give us an eps-level confidence
 * interval for the 90th percentile, or for any quantile in general.
 */
double martingale_cs_threshold(
    uint64_t n, uint64_t min_count, double log_eps);

/*
 * Returns the width of a `1 - exp(log_eps)`-confidence interval for
 * the sum of `n` i.i.d. values sampled from `X`, where `X` has a
 * zero mean and a range of the form `[lo, lo + span]`.
 */
double martingale_cs_threshold_span(
    uint64_t n, uint64_t min_count, double span, double log_eps);

/*
 * Returns the width of a one-sided `1 - exp(log_eps)`-confidence
 * interval the sum of `n` i.i.d. values sampled from `X`, where `X`
 * has a zero mean and range `[lo, hi]` (`lo <= 0 <= hi`).
 *
 * The one-sided confidence interval is on [Sum X_i \leq Interval_n]
 * with `p = 1 - exp(log_eps)` for all n, where `Interval` is the
 * return value of `martingale_cs_threshold_range`.
 *
 * This function finds tighter intervals than the `_span` variant when
 * `|lo| > |hi|`: in that case, a sum > 0 means that we experienced a
 * lot of events with small individual impact on the sum.
 * Intuitively, that's less likely to happen than being unlucky once,
 * with that rare event having a large impact on the sum.
 *
 * This function differs from `martingale_cs_threshold_span` because
 * it returns stronger asymmetric intervals.  The other half-interval
 * for [Sum X_i >= Interval_n] may be obtained by flipping the sign of
 * the variate and swapping `-lo` and `-hi`.  At least one of the two
 * half-intervals will match the symmetric one returned by `_span`,
 * but the other one is tighter (unless the range is symmetric, in
 * which case both are equal).  The more negative `lo / (hi - lo)`,
 * the smaller the interval returned by this function.
 */
double martingale_cs_threshold_range(
    uint64_t n, uint64_t min_count, double lo, double hi, double log_eps);

/*
 * Uses the martingale confidence sequence to return the width of the
 * confidence interval on the index of a given `quantile` in `n`
 * observations, as described earlier.
 *
 * `quantile` must be in [0, 1].  The other parameters are as in
 * `martingale_cs_threshold`.
 *
 * Given a bag of `n` observations, and assuming that we only compute
 * quantile confidence intervals once we have `min_count` observations
 * or more, this function will return the slop for the index of the
 * `quantile`, for a `1 - exp(log_eps)` confidence interval.  The
 * quantile should be between the value at index `floor(quantile * n -
 * slop)` in the sorted list of observations and at at index
 * `ceiling(quantile * n + slop)`.  Either index might be out of
 * bounds (or negative); in that case we have too few observations to
 * provide a lower or upper bound on the quantile.
 *
 * The martingale confidence sequence guarantees that the actual
 * distribution quantile lies in that interval *for every n* with
 * probability at least `1 - exp(log_eps)`.
 */
double martingale_cs_quantile_slop(
    double quantile, uint64_t n, uint64_t min_count, double log_eps);

/*
 * Returns the upper end of an asymmetric confidence interval for
 * `quantile`: with probability `1 - exp(log_eps)`, the true quantile
 * is always less than or equal to `quantile * n + slop_hi`
 * observations.
 *
 * This half-interval has the same properties as that returned by
 * `martingale_cs_quantile_slop`, but is tighter in one direction when
 * `quantile != 0.5`, and otherwise equal.
 */
double martingale_cs_quantile_slop_hi(
    double quantile, uint64_t n, uint64_t min_count, double log_eps);

/*
 * Returns the lower end of an asymmetric confidence interval for
 * `quantile`: with probability `1 - exp(log_eps)`, the true quantile
 * is always greater than or equal to `quantile * n + slop_lo`
 * observations.
 *
 * This half-interval has the same properties as that returned by
 * `martingale_cs_quantile_slop`, but is tighter in one direction when
 * `quantile != 0.5`, and otherwise equal.
 */
double martingale_cs_quantile_slop_lo(
    double quantile, uint64_t n, uint64_t min_count, double log_eps);
#ifdef __cplusplus
} /* extern "C" */
#endif
#endif /* !MARTINGALE_CS_THRESHOLD */
