Martingale Confidence Sequence
==============================

This is an implementation of the martingale confidence sequence
underlying Darling and Robbins's
[Confidence sequences for mean, variance, and median](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC335597/).

This confidence sequence method is appropriate to reject the null
hypothesis that a sum of i.i.d. random values in [-1, 1] has zero mean
(i.e., that the i.i.d. random values themselves have zero mean).  This
test is actually more widely applicable: it suffices for the moment
generating function `mgf(t) = E[exp(tx)] <= exp(t^2 / 2)` for all
`t >= 0`.  However, this property can be hard to prove, and
[Hoeffding's lemma](https://en.wikipedia.org/wiki/Hoeffding%27s_lemma)
says that a range in [-1, 1] suffices to satisfy the condition.

This library implements a confidence sequence on the rank of any
specific quantile in the observations on top of the martingale
confidence sequence, as demonstrated in the aforementioned paper of
Darling and Robbins.

We could also use this martingale confidence sequence to compare the
mean of two random variables X and Y in [0, 1] (e.g., runtimes).
Their difference Z = X - Y has range [-1, 1] and, if X and Y have the
same mean, Z has a mean of 0.  We can thus accumulate the sum `z_1 +
z_2 + ... + z_i`, and compare against the threshold returned by
`martingale_cs_threshold` at each iteration.  If we observe that the
sum exceeds the threshold even once, we may reject the hypothesis that
X has mean equal to or less than that of Y (a two-tailed test simply
needs a Bonferroni correction by adding `martingale_cs_eq` to
`log_eps`).

A call to `martingale_cs_threshold` generates a confidence sequence at
level `1 - exp(log_eps)`, for a sum of `n` values, assuming that the
first `min_count` values can be accumulated without expecting any
useful confidence interval (`martingale_cs_threshold` returns `+infty`
when `n < min_count`).  The confidence sequence guarantees that, if
the summands have zero mean and range [-1, 1] (or satisfy the
constraint on the mgf), the probability that sum exceeds the value
returned by `martingale_cs_threshold` at a single iteration is at most
`exp(log_eps)`, regardless of the total number of iterations (i.e.,
even if it's unbounded).  Moreover, the interval keeps shrinking, so
any mean that's strictly positive will eventually be detected.

For a two-tailed "equality" comparison (i.e., to determine when the
running sum is too positive *or too negative*), add `martingale_cs_eq`
to `log_eps` in the call to `martingale_cs_threshold`.  This will
simply ask for a confidence interval with half the false positive
rate, so that we can use the threshold symmetrically to check if the
sum is too high or too low, and still guarantee a total false positive
rate of at most `exp(log_eps)`.

See also
--------

The martingale-CS can work with a large number of point statistics,
but tends to need a lot of data to reject the null hypothesis.  If the
question can be cast as a Binomial test, [confidence sequence
method](https://github.com/pkhuong/csm) should terminate much more
quickly.  Otherwise, if we can compare full distributions,
[one-sided-KS](https://github.com/pkhuong/one-sided-ks)
may also be applicable: this Kolmogorov-Smirnov test tends to
require more data points than the binomial CSM, but still
less so than martingale-CS.
