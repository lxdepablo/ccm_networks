# Log: bootstrapping + multivariate CCM

Work log for the task in `claude.md`: (1) bootstrap-based false-positive pruning for the
pairwise CCM pipeline, and (2) a multivariate CCM implementation (MXMap, Zhang et al. 2025,
`zhang25a.pdf`) to compare against the pairwise strategy.

## Interpretation of "compare results from the two approaches"

"Approach A" = the existing pairwise CCM pipeline (`generate_xmaps.R` → `xmap_analysis.R`)
plus a new false-positive pruning step. "Approach B" = MXMap's two-phase design: the same
pairwise CCM graph, refined by pruning edges that are fully explained by an indirect
(mediated) path (multiPCM). `xmap_analysis.R` now produces all three edge lists so they can
be compared directly:

- `data/edge_lists.csv` - unchanged, raw pairwise (convergence-filtered only)
- `data/edge_lists_bootstrap.csv` - Approach A (+ bootstrap significance pruning)
- `data/edge_lists_multivariate.csv` - Approach B (+ multiPCM indirect-edge pruning, applied
  on top of the bootstrap-pruned graph)

`code/compare_networks.R` reports network-level stats for all three *and* which specific
edges each pruning step removed - bootstrap pruning and multiPCM pruning catch different
kinds of false positives (see below), so counts alone would be misleading.

## Important pre-existing finding: edge orientation

`par_calc_all_xmaps()`/`CCM()` name their output columns `"A:B"` for
`CCM(columns = A, target = B)` - i.e. "use A's manifold to predict B". `xmap_analysis.R`
splits that string on `:` and records the edge as `sp1 (= A) -> sp2 (= B)`.

I verified empirically (simulated `x` autonomous, `y` driven by `x`, i.e. true `x -> y`) that
a high `"A:B"` score means **B's dynamics leave a footprint in A** - CCM's classic result
that the *effect's* manifold reconstructs the *cause*. In the test system, `"y:x"` was high
(≈0.83-0.96) and `"x:y"` was low (≈0.2-0.4), and the true relationship was `x -> y`. So a high
`"A:B"` score actually supports the causal claim **B -> A**, the reverse of how
`xmap_analysis.R` draws the edge (`A -> B`).

This looks like a pre-existing convention issue in the pipeline (note `network_analysis.R`
already applies `reverse_edges()` to the metaweb, but not to the CCM-derived network, which
may or may not be intentional). **I did not change it** - fixing edge orientation across the
whole pipeline is a separate decision outside this task's scope, and it might already be
understood/compensated for downstream. All new code in `edm_utils.R` is written in terms of
`from_col` (= `sp1`, the manifold/`columns` variable) and `to_col` (= `sp2`, the
target/reconstructed variable) - i.e. matching `edge_lists.csv` literally - so it plugs into
the existing edge lists regardless of which way the true causal arrow points. Worth a
deliberate decision from whoever owns this repo before drawing causal conclusions from the
network's edge directions.

## Bootstrap pruning (`bootstrap_ccm_significance()`, `prune_nonsignificant_edges()`)

### What shipped

A block-**permutation** significance test per edge: `to_col`'s blocks (size `block_size`,
default 10) are reordered *without replacement* (`block_permute()`) while `from_col` stays in
its original order, breaking any real temporal coupling while preserving each series' own
within-block autocorrelation. Cross-map skill is recomputed on this surrogate `n_boot` times
(default 200) to build a null distribution; `p_value = (1 + #(null_rho >= rho_obs)) / (1 +
n_boot)`, significant if `p_value < alpha` (default 0.05).

### What I tried first, and why it was wrong

`edm_utils.R` already had an unused `block_bootstrap()` helper (resamples blocks *with*
replacement, preserving the pairing between whichever columns are passed in). My first
version reused it directly as the null generator. Empirically, on a genuinely independent
pair (`x`, `w`, both autonomous chaotic series, no coupling), this flagged the edge as
"significant" regardless of block size (5/10/20/40 all gave CIs excluding zero, and the
"significant" CI-lower-bound got *worse*, not better, as block size grew). Two compounding
problems, both confirmed:

1. Resampling the **joint** `(from_col, to_col)` pair keeps `from_col[i]` paired with
   `to_col[i]` in every row - it cannot test "not coupled" because it never breaks the
   coupling. It's a sampling-variability estimator, not a null-hypothesis surrogate (the
   original docstring's H0 claim was simply wrong for what the code did).
2. Sampling blocks *with replacement* creates duplicate blocks; a duplicated point is its own
   nearest neighbor in the simplex projection (distance 0, exact target match), which
   trivially inflates cross-map skill - worse with larger blocks, matching the observed trend.

Switching to a without-replacement block permutation of one series fixed both: verified on
the same independent pair (now correctly not significant, p≈0.10) plus a known-real direct
edge and a known-real-but-indirect (transitive) edge (both correctly significant, p<0.01) -
see "Validation on simulated data" below.

### Runtime

Reruns `Simplex` `n_boot` times per surviving edge per site - this is the part of the pipeline
most likely to be slow at full scale; `xmap_analysis.R` only runs it on edges that already
survived `filter_xmaps()`'s convergence screen, not every pair, but `n_boot` and `block_size`
are exposed as parameters for HPC vs. local tuning.

## Multivariate CCM / MXMap (`multi_pcm()`, `prune_indirect_edges()`)

### What it does

Implements Zhang et al. 2025's two-phase MXMap framework: phase 1 is the existing pairwise
CCM graph; phase 2 (`prune_indirect_edges()`) finds, for each edge `from_col -> to_col`,
other node(s) `k` already sitting on a 2-hop path (`from_col -> k` and `k -> to_col` both
present in the current graph) and runs `multi_pcm()` to decide whether the edge is direct or
fully explained by that indirect path.

`multi_pcm()` is a genuine **two-hop** composition, matching Eq. 7 of the paper exactly:

1. apparent: `to_col` reconstructed from `from_col`'s manifold
2. reconstruct each intermediate *from `from_col`'s manifold* (the X2 → Conds hop)
3. conditioned: `to_col` reconstructed from that *reconstructed* conds block (the Conds → X1
   hop) - **not** from the true conds values directly, which would over-prune real direct
   links whenever conds also happens to correlate with `to_col`
4. `rho_direct` = partial correlation of (`to_col`, apparent-reconstruction) controlling for
   (conditioned-reconstruction); `ratio = rho_direct / rho_all`; prune if `ratio <
   gamma_threshold` (default 0.45, the paper's own value)

`knn` defaults to 10, matching the paper's own choice - chained (two-hop) reconstruction is
noticeably more noise-sensitive than a single cross-map, and rEDM's own default (`knn = E+1`)
wasn't enough to get clean separation in testing (see below).

### Adaptation for short series (deliberate, documented)

The paper's multivariate embedding (multiSSR, Eq. 6) stacks a full E-dimensional lag
embedding *per conditioning variable*. With ~26 timepoints per site that blows the
conditioned manifold's dimensionality past what the library can support, so here each
conditioning variable contributes a single (unlagged) reconstructed value to the conditioned
block instead of a full E-dim lag stack - the conditioned embedding dimension is
`length(conds)` regardless of E. `max_conds` (default 3) caps how many intermediates get
conditioned on, for the same reason.

### Validation on simulated data

Built a simulated chain `x -> y -> z` (x autonomous, y driven by x, z driven by y; logistic
maps, matching the style of the paper's own Lotka-Volterra-derived test systems) plus an
independent series `w`, following the paper's own multiPCM validation recipe (Section 4.2:
noise-free systems, and note their multiPCM-specific runs use L=3500, much longer than the
L≤6000 used for their main causal-discovery benchmarks - I matched L=3500 for this check
after finding L≈250-700 gave too little reconstruction fidelity to separate direct vs.
indirect at all).

Ratio results (`E=3, tau=1, knn=10`, at L=3500):

| Edge | Relationship | ratio |
|---|---|---|
| `y -> x` (recovers true `x -> y`) | direct | ~1.000 |
| `z -> x` (recovers true, transitive `x =>...=> z` via y) | **indirect** | ~0.975-0.985 |

This holds directionally (indirect ratio consistently below direct ratio) across a small
grid of `(E, tau)` values (E ∈ {3,5,8}, tau ∈ {1,2}), i.e. the method correctly *ranks*
indirect below direct, but the indirect ratio never dropped below the paper's own
`gamma_threshold = 0.45` in this coupling regime (`beta = 0.35-0.5`). The paper's own
Appendix C makes exactly this point about the threshold: it was tuned on noise-free
Lotka-Volterra-style systems at their own coupling strength and "the appropriate threshold
may vary depending on the system." I did not chase a threshold that would "work" on my
specific simulated system, since that would just be curve-fitting one more simulation rather
than validating the method - the directional check (does it correctly rank indirect below
direct, and does the bootstrap test correctly separate independent/real/indirect-but-real) is
the more meaningful validation.

**Where I initially went wrong, for anyone extending this**: a first attempt conditioned on
the *true* (not reconstructed) intermediate values, contemporaneously, and saw *zero*
discrimination between direct and indirect edges (ratio ≈ 1.0 either way). That result is a
broken test, not evidence the method doesn't work: contemporaneous linear correlation between
chaotic variables is near zero even under strong coupling (the "mirage correlation"
phenomenon CCM exists to route around), so partialling on it can't move a partial correlation
computed from a mirage-correlated conditioning variable. The paper never conditions on raw
values for this reason - always on a cross-map reconstruction, as implemented above.

### Real-data limitation (documented, not "fixed")

Per-site series here are **~26 timepoints** (confirmed: `prep_site_ccm_data()` on real site 1
gives a 26 x 35 wide table), vs. the paper's own L=3500 validation length. The two-hop
composition is noise-sensitive even at L=3500-4000 in testing above; at L=26 it will be far
less reliable, and MXMap's own Algorithm 1 conditioning on "all intermediate nodes along any
path" isn't feasible at this scale regardless of embedding tricks. This is a genuine data
limitation, not an implementation gap - `multi_pcm()`/`prune_indirect_edges()` are faithful
to the paper and validated on simulated data; on the real per-site data they should be treated
as a secondary, lower-confidence view, not the primary result.

A small-subset preliminary run against real data (site 1, top-5-variance species + sst,
N=26) completed end to end with no errors and produced one interpretable pruning decision:
`brown_cushion -> brown_thin_blade_d` (ratio 0.148, well under 0.45) was flagged as indirect,
mediated by `acsp` (which had direct edges to both). This confirms the pipeline runs on real
data and produces sensible output; it is not a claim that this specific edge is correctly
classified given N=26's low statistical power - see above.

### Bottom line for "compare the two approaches"

- **Bootstrap-pruned pairwise network** (`edge_lists_bootstrap.csv`) is the primary, more
  reliable output at this series length.
- **Multivariate/MXMap network** (`edge_lists_multivariate.csv`) is implemented faithfully to
  the paper, validated on simulated data (correct ranking of indirect vs. direct; threshold
  needs system-specific calibration per the paper's own caveat), and available as a secondary
  comparison - but its two-hop reconstruction is data-hungry, and ~26 points/site is short
  enough that its pruning decisions should be treated cautiously until run on more data (e.g.
  pooling across a wider date range, or lower-frequency validation on simulated systems at
  comparable N).
- `compare_networks.R` reports which edges differ between all three networks per site, so this
  comparison can be inspected directly once run on the full HPC-scale data.

## What was NOT run locally (per claude.md's instruction)

The full pipeline (`generate_xmaps.R` over all 8 sites x ~30-46 species/site, then bootstrap
pruning with `n_boot=200` and multiPCM pruning over every surviving edge) was **not** run
locally - it's the same CCM computation the task says is HPC-scale, now with bootstrap and
multiPCM adding more compute on top. What was run locally:

- Package installation + empirical verification of rEDM's CCM/Simplex conventions (small
  simulated data, informed the edge-orientation finding above)
- Full validation of bootstrap pruning and multiPCM on simulated systems (documented above)
- One small-subset preliminary pass on real data (1 site, 6 variables, reduced `n_boot`/`maxE`)
  to confirm the wiring works end-to-end without crashing

`generate_xmaps.R`/`xmap_analysis.R` are otherwise structured to run as before (same
`setwd()` HPC path, same inputs) - the new bootstrap/multivariate steps are additive stages
at the end of `xmap_analysis.R`, and `generate_xmaps.R` now also writes each site's prepped
wide time series to `data/ccm_data/site_<n>.csv` so those new steps don't need to redo the
data wrangling or re-run CCM.

## Files touched

- `code/edm_utils.R` - new: `prep_site_ccm_data()` (refactored out of `generate_xmaps.R`,
  fixes a latent bug where a later `select(-survey_group)` referenced an already-dropped
  column), `get_optimal_E()`, `xmap_reconstruct()`, `multi_xmap_reconstruct()`,
  `partial_cor()`, `multi_pcm()`, `prune_indirect_edges()`, `block_permute()`,
  `bootstrap_ccm_significance()`, `prune_nonsignificant_edges()`
- `code/generate_xmaps.R` - uses `prep_site_ccm_data()`; saves per-site wide `ccm_data` to
  `data/ccm_data/site_<n>.csv`
- `code/xmap_analysis.R` - adds the bootstrap and multivariate pruning stages, writing
  `data/edge_lists_bootstrap.csv` and `data/edge_lists_multivariate.csv`
- `code/compare_networks.R` - new: compares the three networks, reporting which edges each
  pruning step removed and per-site network stats
- `code/network_analysis.R`, `code/sst_scripts/pull_sst.R` - pre-existing uncommitted changes
  from the user's working copy, carried into this worktree unmodified so this branch reflects
  the actual current state (not otherwise touched by this task)
