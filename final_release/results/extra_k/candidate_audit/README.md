# Circle rank study and probability model

The study combines the exact local rank identity with a rank-only sample-size experiment. It measures closure at connectivity and after extra integer-k steps on uniform, cosine and trough profiles. The probability model describes independent one-gap marks; it does not establish a universal law for dependent circle gaps.

## What the local identity proves

For a neighboring arc gap G_i < 1/2 on a clean unit circle, let L_i count points on the backward arc of length G_i from its first endpoint, and R_i count points on the forward arc of that length from its second endpoint. The intervals are open and exclude endpoints. With distinct chords, no distance tie needs treatment.

There are no points inside the gap. Strict chord monotonicity with minor-arc length therefore makes the endpoint ranks 1 + L_i and 1 + R_i. The OR key is 1 + min(L_i, R_i); the mutual key is 1 + max(L_i, R_i). Theorem 2 turns the two largest neighboring keys into k_conn and k_close. This identity computes the rank keys used by the sample-size study.

The implementation checks 18,980 neighboring keys against the released full chord-rank kernel on 240 samples, and checks 480 actual graph-connectivity events. All agree. A separate 30-sample check at n=128 and 512 also compares the vectorized sampler and its 19,200 keys with the sampler and full chord ranks; all keys agree. The largest coordinate discrepancy in that additional check is 8.71e-13.

## Experiment

The study uses the uniform, cosine and quarter-circle trough cumulative distributions. Sizes and samples per profile are:

- n=128, 512 and 4,096: 1,500 seeds each;
- n=32,768: 500 seeds;
- n=262,144: 100 seeds;
- n=1,000,000: 30 seeds.

The same uniforms and phase are used across profiles at a given size/seed. There are **15,390 profile–seed samples from 5,130 shared seed–size combinations**, yielding 30,780 method records. This is not 15,390 independent seed streams. Sampling seeds are `20330101 + 10000 + 100000 * size_index + replicate`, with the sizes in the order above. The 240 small validation samples use seeds starting at 20330101.

Sorted cyclic counts give the rank keys directly, avoiding a quadratic pair-distance matrix. Cosine inversion uses bracketed vectorized Newton steps with CDF residual at most 2e-13; the trough inverse is piecewise linear. Small/full-rank validation checks the resulting ranks, rather than coordinate agreement alone. All sweep samples have g < 1/3 (maximum 0.271062), so the positive-reference-loop and integer-k theorem applies. No new Ripser diagrams or diagram-error experiments are computed.

The data support these statements:

- Uniform n=128 through 32,768: closure at connectivity stays between 31.2% and 34.3% for OR, and between 28.7% and 29.4% for mutual. More samples do not make connectivity a reliable closure rule in these tested settings.
- At n=512, k_conn+3 closes 94.7% of uniform OR samples (pointwise 95% Wilson interval 93.5–95.8%) and 91.1% of uniform mutual samples (89.5–92.4%).
- At n=512 in the trough, the corresponding rates are 94.1% for OR (92.8–95.2%) and **38.1% for mutual** (35.6–40.6%). The largest extra-k requirements are 25 and 71, respectively.
- At n=128, the trough OR rate is only **71.1%** (68.8–73.4%). Thus “OR is unchanged by a density jump” and “add three works at any n” overstate even this numerical evidence.
- Cosine n=128 gives a mutual +3 closure rate of 78.7%, rising to about 91% at n=4,096. Smooth density does not remove finite-size effects.

The manuscript includes the identity, a practical paragraph and Table D.2 for all six sizes, with pointwise 95% Wilson intervals on the two exploratory large sizes. Profiles reuse uniform draws and phase: at n=1,000,000 uniform and cosine have identical closure indicators seed by seed for both methods, and trough OR matches them. The raw-record and interval checks are in results/rank_scaling_audit.json. The full extra-k distribution, including j=0 through 5 and tails beyond five, remains in the sample records. Intervals are pointwise, conditional on a profile and size; profiles are paired and should not be pooled as independent replications.

The defensible practical message is: **connectivity can leave a gap even at large sample sizes; a small extra increment often closes uniform OR graphs, while mutual needs substantially more steps across the tested density jump.** The data do not guarantee an increment, establish an optimal method across other densities, or supply a general noisy-data tuning rule.

## The Poisson calculation: correct marginal tails, incomplete global law

On an infinite homogeneous Poisson line, condition on a consecutive gap X ~ Exp(1). The counts on its two exterior sides are independent Poisson(X). After integrating out X,

P(L=l, R=r) = (l+r)! / [l! r! 3^(l+r+1)].

Hence P(L >= k)=2^-k. For k>=1 the OR tail is

T_OR(k) = integral_0^infinity exp(-x) P(Poisson(x)>=k)^2 dx
        = 2^(1-k) - 2 sum_(a=0)^(k-1) binom(k+a-1,a)/3^(k+a).

The mutual tail is T_mu(k)=2^(1-k)-T_OR(k), with both tails at zero equal to one. The script uses exact rational arithmetic for the finite sum before converting to floating point. Twelve independent numerical integral checks agree within 1.11e-16.

A positive version is T_OR(k)=2^(1-k) P(Binomial(2k-1, 2/3)<=k-1). The last term in that binomial tail dominates up to a factor tending to two: consecutive term ratios, moving downward from k-1, are bounded by 1/2 and approach 1/2 for a fixed number of steps. The central-binomial estimate then gives

T_OR(k) ~ 3/sqrt(pi*k) * (4/9)^k.

Consequently the OR tail ratio tends to 4/9 and the mutual ratio to 1/2. These are correct **one-gap Poisson statements**, not exact finite-circle probabilities.

For *independent* copies Z_1,...,Z_n of a mark with T(m)=P(Z>=m) and p_m=T(m)-T(m+1), the exact extra-k CDF is

P(max Z - second-largest Z <= j)
  = 1 - n sum_(m=j+1)^infinity p_m [1-T(m-j)]^(n-1).

Replacing the bracket by exp(-n*T(m-j)) gives the Poissonized approximation. At n=512 the exact independent-model CDFs for j=0..5 are:

- OR: 0.32877, 0.71449, 0.87804, 0.94773, 0.97753, 0.99032;
- mutual: 0.27492, 0.63368, 0.81524, 0.90694, 0.95318, 0.97647.

They are close to the uniform-circle simulation. Nevertheless, ranks of different gaps on one circle share surrounding points and are dependent. Marginal Poisson tails do not prove that the global two largest ranks behave like independent marks. An extreme-event dependence or clustering argument is still missing.

For the independent geometric heuristic, **(1-q)/log(1/q) is the phase-averaged unique-maximum probability**. The tie probability is its complement, **1-(1-q)/log(1/q)**, giving 0.314916 for q=4/9 and 0.278652 for q=1/2. The classical independent geometric model has lattice oscillations; these numbers are not an exact sample-size-independent distribution. The OR tail also has a k^(-1/2) factor, so substituting its limiting ratio into a pure geometric model loses finite-size information.

A density jump's finite-size deviations do not themselves prove a different limiting law as n tends to infinity. Neither direction is established here. The manuscript therefore does not promote the probability model to a theorem.

## Files and reproduction

Run `python3 code/check_knn_miss_candidate.py` from the release. To recompute all samples rather than resume, move the existing `results/extra_k/candidate_audit/samples.csv` and `summary.csv` aside together, then run the command. `--verify-only` checks all summaries and model values and replays one seed per method/profile/size group. It currently passes 949 checks.

- `samples.csv`: every paired sample's connecting k, closing k, extra k and largest gap;
- `summary.csv`: full CDFs, Wilson intervals, medians and maxima;
- `config.json`: seeds, sizes, sampling and validation scope;
- `iid_predictions.csv` and `poisson_tails.csv`: independent-model computations;
- `manuscript_table.tex`: the inline appendix table's reproducible numerical content;
- `verification.json`: the independent tail integral and stored-summary checks.
