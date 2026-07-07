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

## Edge orientation fix

`par_calc_all_xmaps()`/`CCM()` name their output columns `"A:B"` for
`CCM(columns = A, target = B)` - i.e. "use A's manifold to predict B". `xmap_analysis.R` used
to split that string on `:` and record the edge as `sp1 (= A) -> sp2 (= B)`.

I verified empirically (simulated `x` autonomous, `y` driven by `x`, i.e. true `x -> y`) that
a high `"A:B"` score means **B's dynamics leave a footprint in A** - CCM's classic result
that the *effect's* manifold reconstructs the *cause*. In the test system, `"y:x"` was high
(≈0.83-0.96) and `"x:y"` was low (≈0.2-0.4), and the true relationship was `x -> y`. So a high
`"A:B"` score actually supports the causal claim **B -> A**, the reverse of how the edge used
to get drawn (`A -> B`).

**Fixed**: `xmap_analysis.R` now records the edge as `sp1 (= B = cause) -> sp2 (= A = effect)`
- i.e. `edge_lists.csv`'s `sp1 -> sp2` now points in the direction of causality, not the raw
`"columns:target"` order. Re-validated on the same simulated chain (`x -> y -> z` + an
independent `w`) after the fix: edges now come out as `x -> y`, `y -> z` (both direct, correct
orientation) and `x -> z` (the correctly-identified transitive/indirect edge) - previously
these came out reversed (`y -> x`, `z -> y`, `z -> x`).

This also resolves a related inconsistency: `do_smap.R` already treats `sp1` as the "driver"
of `sp2` (`these_drivers <- filter(this_site, sp2 == this_target)$sp1`) when fitting its
multivariate S-map - i.e. it already assumed `sp1 -> sp2` meant "`sp1` causes `sp2`". Before
this fix that assumption was actually backwards relative to what `edge_lists.csv` encoded; now
it's correct, with no change needed to `do_smap.R` itself. `network_analysis.R` similarly
needed no changes - it just reads `edge_lists.csv` and builds a directed graph from `sp1`/
`sp2`, so it now produces a correctly-oriented causal network automatically. (Note
`network_analysis.R` separately applies `reverse_edges()` to the *metaweb* - that's an
unrelated, pre-existing convention for the trophic metaweb's own edge direction, not touched
here.)

Since `sp1`/`sp2` now mean cause/effect rather than columns/target, the pruning functions in
`edm_utils.R` (`prune_indirect_edges()`, `prune_nonsignificant_edges()`) needed a
corresponding fix: they use `sp1`/`sp2` as literal cause/effect for graph-topology work
(finding 2-hop mediators), but must still call `multi_pcm()`/`bootstrap_ccm_significance()`
with the manifold variable and reconstructed variable in the right (reversed) order -
`from_col = effect`, `to_col = cause` - since that's the actual direction the underlying
`Simplex` call needs to run in. Getting this backwards would silently retest the wrong
direction. Both were updated and re-validated (see below).

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

Both `prune_nonsignificant_edges()` and `prune_indirect_edges()` (below) parallelize across a
site's edges via `future_map()`/`plan(multisession)`, matching the convention already used by
`par_calc_all_xmaps()` and `do_smap.R` - parallel within a site, sequential across sites
(unchanged in `xmap_analysis.R`), so there's no nested parallelism to manage. Re-validated
after parallelizing: identical results to the serial version on the simulated chain.

While parallelizing, found and fixed an unrelated inefficiency: `EmbedDimension()` defaults to
`showPlot = TRUE`, so every parallel worker was opening a graphics device (visible as
"MultisessionFuture ... opened the default graphics device" warnings, and would otherwise
leave orphaned `Rplots.pdf` files scattered across worker processes on the cluster). Fixed in
`get_optimal_E()` and in the pre-existing `EmbedDimension()` calls in
`par_calc_all_xmaps()`/`calc_all_xmaps()`.

## Multivariate CCM / MXMap (`multi_pcm()`, `prune_indirect_edges()`)

### What it does

Implements Zhang et al. 2025's two-phase MXMap framework: phase 1 is the existing pairwise
CCM graph (`sp1 = cause -> sp2 = effect`, per the edge orientation fix above); phase 2
(`prune_indirect_edges()`) finds, for each edge `cause -> effect`, other node(s) `k` already
sitting on a 2-hop causal path (`cause -> k` and `k -> effect` both present in the current
graph) and runs `multi_pcm()` to decide whether the edge is direct or fully explained by that
indirect path. `multi_pcm()` itself works in *Simplex* terms, not causal-graph terms - it
takes `from_col` (the manifold/columns variable) and `to_col` (the variable being
reconstructed), so it's called as `multi_pcm(..., from_col = effect, to_col = cause, conds,
...)` - reconstructing the cause from the effect's manifold, matching how the edge's apparent
CCM score was computed in the first place.

`multi_pcm()` is a genuine **two-hop** composition, matching Eq. 7 of the paper exactly:

1. apparent: `to_col` (cause) reconstructed from `from_col`'s (effect's) manifold
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

**Post-orientation-fix re-check**: the table above was produced by calling `multi_pcm()`/
`xmap_reconstruct()` directly with explicit `from_col`/`to_col`, independent of the
`sp1`/`sp2` edge-list convention, so it wasn't affected by the edge-orientation bug or its
fix. After fixing `prune_indirect_edges()` to map `sp1`/`sp2` (cause/effect) onto
`from_col`/`to_col` in the (reversed) order `multi_pcm()` actually needs, I re-ran the full
chain simulation end-to-end through the *fixed* pipeline functions (edges built the same way
`xmap_analysis.R` now builds them, then pruned via `prune_indirect_edges()`/
`prune_nonsignificant_edges()`): edges came out correctly oriented (`x -> y`, `y -> z` direct;
`x -> z` transitive) and reproduced the same ratios as above (`x -> y` ratio ≈0.999, `x -> z`
ratio ≈0.991) and the same bootstrap significance results (all four edges significant,
p<0.01) - confirming the fix didn't change the underlying math, only which physical variables
get passed to it.

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
  `data/edge_lists_bootstrap.csv` and `data/edge_lists_multivariate.csv`; fixes edge
  orientation so `sp1 -> sp2` points in the direction of causality (see above)
- `code/compare_networks.R` - new: compares the three networks, reporting which edges each
  pruning step removed and per-site network stats
- `code/network_analysis.R`, `code/sst_scripts/pull_sst.R` - pre-existing uncommitted changes
  from the user's working copy, carried into this worktree unmodified so this branch reflects
  the actual current state (not otherwise touched by this task)
