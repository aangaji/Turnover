# Turnover
## Theory

> #### lineage turnover
> A mutant offspring lineage or clone is **orphaned** if it did not go extinct itselft <br>
> while its parental lineage (excluding the offspring lineage) went extinct. <br>
> W(q, N) = q / (2 log(N)) * (1-N^(-2))


> #### clone turnover
> A mutant offspring lineage or clone is **estranged** if it did not go extinct itself <br>
> while its parental clone went extinct. <br>
> W(b; a, μ, T) = μ (b + μ^2) / (1-μ)^2 / (a(1-μ/2)-b) (1-exp(-2*(a-b-aμ)T)) / (1-exp(-2aμT))

## Simulation

We simulated the exponential growth of mutant `clones` within a population of tumor cells with a Gillespie algorithm of cell division and death. A cell divides with a probability proportional to the birth rate `a` or dies with rate `b`. At birth any of the two cells can mutate with probability `μ` thus founding a new clone.

To measure the turnover from simulations a tumor population is first grown to a threshold size `N` or time `T`. The tumor then grows to a size that is sufficiently high for the parental lineages / clones to go extinct. In the case of lineage turnover the mutation rate can be set to zero for this growth phase.
The fraction of orphaned / estranged lineages that appeared up to the threshold gives the turnover. The results are averaged over a few hundred (between 200 and 1200) simulated tumors depending on mutation rate and threshold size.

For the clone turnover we additionally measured the expected fraction of estranged lineages. During simulation at birth of a mutation the current size `n` of its parental clone is recorded from which the probability of becoming estranged is computed as `q^n`, where `q` is the extinction probability for clones. Setting a time threshold instead of size slightly improves the fit to the analytical curve as expected. In this case no second growth phase is needed.

Determining whether a mutant lineage is orphaned or estranged is straight forward with knowledge of the complete genealogical tree of a all mutations. If the order of occurrence is unknown this is not possible in the case of subsequent mutations that become clonal to a lineage. Nevertheless the fraction of orphaned / estranged lineages and thus the turnover can still be computed exactly with an appropriate algorithm.

We also considered the number of mutations arising at division to be distributed as `~ Poisson(μ)`. If mutations are said to be single nucleotide variations (SNVs), this is a more realistic model. However, the results do not differ considerably at small mutation rates. The deviation in both lineage- and clone turnover from Poisson distributed mutation numbers is highest at death rate *d=0* and decays for increasing *d*.

## Results
<img src="turnover_plots/theory/orphaned_turnover.png" alt="lineage turnover" width="350"/>

![lineage turnover](turnover_plots/theory/orphaned_turnover_series.png)  
**Lineage turnover** - Fractions of orphaned lineages from simulations (blue markers) plotted together with the analytical turnover curve (red line) for different mutation rates. Lineage turnover is independent from mutation rates and the measured curves coincide. Simulations were run with 200 to 800 repetitions, birth rate *a=1.*, a threshold of 2000 cells and final size of 6000. Error bars indicate standard deviation.

<img src="turnover_plots/theory/estranged_turnover_expected.png" alt="Expected clone turnover" width="350"/>

![Expected clone turnover](turnover_plots/theory/estranged_turnover_expected_series.png)  
**Expected clone turnover** - Expected fraction of estranged lineages / clones from simulations (blue markers) plotted together with the analytical turnover curve (red line) for different mutation rates. A lineage's / clone's expected probability of becoming estranged to the parental clone is computed as `q^n`, with `q` being the clone extinction probability and `n` parental clone size at birth. Simulations were run with 200 to 500 repetitions, birth rate *a=1.* and time thresholds corresponding to a threshold size of 200 cells i.e. *T=log(200)/(a-b)*. Error bars indicate standard deviation.

<img src="turnover_plots/theory/estranged_turnover_measured.png" alt="Measured clone turnover" width="350">

![Measured clone turnover](turnover_plots/theory/estranged_turnover_measured_series.png)  
**Measured clone turnover** - Fractions of estranged lineages / clones from simulations (blue markers) plotted together with the analytical turnover curve (red line) for different mutation rates. Simulations were run with 400 to 1200 repetitions, birth rate *a=1.*, a threshold of 200 cells and final size of 2000. Error bars indicate standard deviation.

## Inference
Using the analytical results we can infer both death rate `b` and mutation rate `μ` for a given tumor and set frequency threshold `N`. Measuring its `lineage turnover` gives us an estimate of the lineage extinction probability. <br>
`q = b/a = 2log(N)(1-N^(-2))^-1*W_l ≈ 2log(N)*W_l`

As observed for simulations measured lineage turnover deviates from the theoretical result for high death rates `b` and low threshold size `N`. In cases where the estimated lineage extinction probability `q` does not fall below one the mutation rate cannot be inferred. However, using the results from simulations for the lineage turnover might allow a better estimation of `q`.

Using the inferred value of `b` there are two ways of obtaining the mutation rate `μ` by the clone turnover. The first and straight forward way is to determine the value of `μ` at which the theoretical clone turnover curve for given `b` assumes the measured turnover. Here the threshold time is to be estimated by `T=log(N)/(a-b)`. The other option is to (uniformly) subsample the set of sequenced mutations or consider smaller parts of the genome and thus effectively reduce the mutation rate. A mutation is selected with probability `L`, for each value of `L` the tumor can be sampled several times and the mean clone turnover computed. One obtains a pseudo clone turnover curve which can be fitted against the analytical result. Although the sets of mutations all share the same background and the curve keeps the tumors deviation from the expected clone turnover the fit generally results in a better match.  

Note: With the estimated lineage extinction probability `q` the (scaled) mutation rate `μ` can also be determined from the cumulative distribution of inverse frequencies (sottoriva plot). For an exponentially growing tumor under neutral evolution this curve is linear and the slope is given by `μ/(1-q)`.


![Inference - `b` dependence](turnover_plots/inference/vary_d_solve.png)  <br>
**Inference - `b` dependence** - Inference of death `b` and mutation rate `μ` for single tumors at fixed `μ`. Mutation rate is inferred using the estimated `b` by computing the tumors clone turnover `W_c` and comparing to the theoretical curve. For each death rate we simulate 10 tumors at `a = 1.0`, `μ = 0.2`, threshold size `N = 2000` and final size 100000. However, at `b = 0.7` two and at `b = 0.8` four of these show high lineage turnover such that `b` cannot be estimated and are excluded. Error bars indicate standard deviation.

![Inference with fit - `b` dependence](turnover_plots/inference/vary_d_fit.png)  <br>
**Inference with fit- `b` dependence** - Inference of death `b` and mutation rate `μ` for single tumors at fixed `μ`. Mutation rate is inferred using the estimated `b` by subsampling mutations thus effectively reducing the mutation rate, computing the clone turnover for each subset and fitting the resulting pseudo-turnover curve to the analytical curve. For each death rate we simulate 10 tumors at `a = 1.0`, `μ = 0.2`, threshold size `N = 2000` and final size 100000.  However, at `b = 0.7` two and at `b = 0.8` four of these show high lineage turnover such that `b` cannot be estimated and are excluded. The set of mutations is subsampled at probabilities `L = 0.1, 0.2..., 1.0` with 10 repetitions each. Error bars indicate standard deviation.

![Inference - `μ` dependence](turnover_plots/inference/vary_mu_solve.png)  <br>
**Inference - `μ` dependence** - Inference of death `b` and mutation rate `μ` for single tumors at fixed `b`. Mutation rate is inferred using the estimated `b` by computing the tumor's clone turnover `W_c` and comparing to the theoretical curve. For each death rate we simulate 10 tumors at `a = 1.0`, `b = 0.4`, threshold size `N = 2000` and final size 100000. Error bars indicate standard deviation.

![Inference with fit - `μ` dependence](turnover_plots/inference/vary_mu_fit.png)  <br>
**Inference with fit- `μ` dependence** - Inference of death `b` and mutation rate `μ` for single tumors at fixed `b`. Mutation rate is inferred using the estimated `b` by subsampling mutations thus effectively reducing the mutation rate, computing the clone turnover for each subset and fitting the resulting pseudo-turnover curve to the analytical curve. For each death rate we simulate 10 tumors at `a = 1.0`, `b = 0.4`, threshold size `N = 2000` and final size 100000. The set of mutations is subsampled at probabilities `L = 0.1, 0.2..., 1.0` with 10 repetitions each. Error bars indicate standard deviation.


## Treeless turnover

For the simulated populations the mutations' order of occurrence is known and the computation of both lineage and clone turnover is straight forward. Regarding clone turnover, given an extant mutation one only needs to consider the most recent ancestral clones. For the lineage turnover we shift perspective to ancestral clones. Given a mutant lineage we go over all mutant offspring lineages that appeared on this background and count them as orphaned if they are clonal.
However, for real tumors the genealogy is a priori unknown. A problem arises when a mutant `m_2` appears on and sweeps through its parental lineage `m_1` we cannot tell which one of them is estranged or orphaned. For instance `m_1` could be considered orphaned with respect to `m_2` since its lineage makes up all of the `m_2` lineage.

Nevertheless, even in such cases we know the number of estranged or orphaned lineages and knowledge of the complete genealogy is not required in order to compute turnover. We therefore defined `treeless` algorithms which compute turnover for unsorted haplotypes. The idea behind the method is the following. For a given mutant lineage `m` we first identify its truncal i.e. clonal set of mutations (the `trunk`) and then the private subset of truncal mutations which are unique to the lineage. The latter is a series of sweeps following the branching point of the lineage `m` (the point in the genealogical tree at which the lineage `m` branches off its most recent surviving ancestral lineage). In the case where the subset only contains `m` itself `m` is simply estranged/orphaned or not. But if it contains other clonal mutations we know that only one can be estranged, because the other clones are extinct. Also, all mutants (expect for the first one) in the subset are orphaned with respect to the ones situated before them within the subset. This is because a subset mutation is clonal to all previous subset mutations. In the example from before $m_2$ is orphaned wrt. `m_1` and if a third mutant `m_3` became clonal on the background of `m_2` it would be an orphan to both `m_2` and `m_1`. Then `m` as well as the other mutants of the subset are assigned a weight such that they eventually add up to the right number of estranged/orphaned mutants.  

We can illustrate this explicitly in the following simple example. A population gave rise to the lineages `1, 2, 3, 4, 5` within which the clones `{1}, {1,2}, {1,2,3,4,5} ` survived. Mutation `1` has trunk `{1}` and therefore no ancestor. Mutation `2` has trunk `{1,2}`, private subset `{2}` and is neither estranged nor orphaned. Mutations `3, 4, 5` on the other hand all have trunk `{1,2,3,4,5}` and private subset `{3,4,5}`. The trunk itself survived as a clone, thus one of the three mutations is estranged. The orphan count is `sum^t_{n=0} n = t*(t-1)/2 = 3` where `t` is the subset size.

<img src="drawing.png" alt="Treeless turnover" width="350"/> <br>
**Treeless turnover** - An exemplary genealogical tree illustrating the idea behind the `treeless turnover` algorithm. Mutations `3, 4, 5` (in red) have the same surviving lineage with the same truncal mutations `{1,2,3,4,5}` where `{1,2}` are not private but have their own distinct surviving clones. Even without knowing the order of occurence we find that one out of {3,4,5} is estranged (`5`) and three are orphaned (`4` wrt. `3` and `5` wrt. `3` and `4`).
