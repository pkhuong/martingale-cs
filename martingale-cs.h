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
 * bit 2: internal_constant
 */
int martingale_cs_check_constants(void);

/*
 * Confidence interval sequence for simple martingales.
 *
 * Let X be a random variable with 0 mean and range in [-1, 1].
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
 * confidence interval for the sum is (-infty, threshold).  Add
 * `martingale_cs_eq` to the `log_eps` for a two-sided test: the
 * confidence interval becomes (-threshold, threshold).
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
 */
double martingale_cs_threshold(
    uint64_t n, uint64_t min_count, double log_eps);

#ifdef __cplusplus
} /* extern "C" */
#endif
#endif /* !MARTINGALE_CS_THRESHOLD */
