Okay, this peer review provides excellent, actionable feedback. Below is the complete revised proposal, integrating all suggestions to create a robust and transparent plan for advancing eye-movement similarity analysis.

## Revised Proposal: Advancing Eye-Movement Similarity Analysis for Encoding-Retrieval Studies

**1. Introduction & Motivation**

Current methods for quantifying eye-movement similarity between encoding and retrieval (e.g., cosine similarity on single-scale Fixation Density Maps - FDMs) offer valuable insights but may lack the sensitivity to capture the full richness of gaze reinstatement. This proposal outlines a staged approach to implement and rigorously evaluate increasingly powerful techniques for measuring eye-movement similarity. Our goal is to enhance sensitivity, improve interpretability where possible, and establish a data-driven framework for selecting the optimal metric(s) for understanding gaze-based memory reinstatement.

We will pay close attention to **stimulus similarity foils**, evaluating metric performance not just on identical image pairs but also on same-category, different-instance pictures to ensure our measures capture specific memory reinstatement rather than generic saliency overlap.

**1.1. Gaze-Data Quality Control**
Prior to any similarity analysis, rigorous data quality control will be implemented:
*   **Track-Loss Handling:** Frames with >[X]% track loss within a [Y]ms window (or equivalent pupil-size validity metrics) will be identified. Short periods of track loss (<[Z]ms) may be interpolated if appropriate for the analysis stage; longer periods will lead to fixation exclusion or trial rejection.
*   **Exclusion Thresholds:** Trials with <[N] valid fixations or total gaze data below [P]% of trial duration will be excluded. Participants with >[Q]% trial exclusion will be considered for removal from the analysis.
*   **Calibration:** Re-calibration procedures will be performed between encoding and retrieval blocks (or at fixed intervals) to maintain high spatial accuracy. Drift correction will be applied if necessary.

**2. Overarching Strategy & Guiding Principles**

Our approach is built on the following principles:

*   **Incremental Complexity:** We will progress from simpler, drop-in upgrades to more sophisticated, novel methods, ensuring each step justifies its added complexity.
*   **Quantitative Gating:** Pre-defined statistical thresholds for improvement will guide decisions to proceed to more complex stages.
*   **Interpretability Balance:** While seeking maximum sensitivity, we will prioritize interpretable methods and include strategies to "peek inside" more black-box approaches.
*   **Rigorous Validation:** A comprehensive validation framework, including noise ceiling estimation and qualitative error analysis, will be applied consistently.
*   **Participant-Level Sensitivity:** We will explore and account for individual differences in optimal metric parameters.
*   **Task Variant Consideration:** This proposal assumes [e.g., full-image re-presentation during retrieval]. If multiple task variants are used (e.g., full-image vs. partial-cue delay), metrics will be evaluated separately for each. Metrics will only be pooled across tasks if a pilot analysis demonstrates that metric × task interactions are negligible for the primary outcomes.

**3. Methodological Stages & Evaluation**

*(A graphical concept figure will be developed to visually illustrate scanpaths and their representation at each stage: FDM, multi-scale FDM, DTW path, graph representation. This will appear as Figure 1 in the final methods.)*

**Stage 0: Establish Current Baseline**
*   *(Effort ≈ 1 week; RAM ≈ 4 GB)*
*   **Method:** Re-implement the current standard: Fixation Density Map (FDM) smoothed by a single σ, with cosine similarity.
*   **Evaluation:** Reproduce published/internal effect sizes for key behavioural correlates (e.g., recognition accuracy, vividness). This serves as the primary benchmark.

**Stage 1: Upgrade the Heat-Map Metric (Same Representation, Better Distance)**
*   *(Effort ≈ 2 weeks; RAM ≈ 8 GB)*
*   **Improvements:**
    1.  **Earth-Mover's Distance (EMD/Wasserstein-1):** Replace cosine similarity with EMD between FDMs. Optionally, sign EMD by subtracting a per-image saliency map (e.g., using a pre-trained model like DeepGaze III, ITTI, or similar, specified *a priori*). Scene-dependent failures of the saliency model will be monitored during the qualitative error-triage.
        *   *Tool:* `transport::wasserstein()` in R.
    2.  **Multi-Scale Similarity:** Compute FDMs at multiple spatial scales (e.g., σ = 0.25°, 0.5°, 1.0°) and average the (1 – EMD) similarity across scales.
    3.  **Duration-Weighted Density:** Weight fixations by their duration (or inverse saccade velocity) when creating density maps.
*   **Evaluation:**
    *   Compare multi-scale, duration-weighted EMD against Stage 0 using permutation z-scores and ΔAUC (macro) for same-image classification.
    *   **Quantitative Gate:** Proceed if ΔAUC ≥ `max(0.02, 0.05 * baseline_AUC_standard_error)` and the 95% CI of ΔAUC > 0.
    *   A power analysis (via simulation) will be pre-registered to ensure adequate power (e.g., 80%) to detect this ΔAUC gate with the expected number of participant-item pairs.
*   **Interpretability:** EMD remains relatively interpretable.

**Stage 2: Bring Back Time – Soft Sequence Alignment**
*   *(Effort ≈ 2-3 weeks; RAM ≈ 8-16 GB)*
*   **Improvements:**
    1.  **Soft Dynamic Time Warping (soft-DTW):** Apply soft-DTW to raw (x,y) fixation sequences. The DTW cost function (e.g., Euclidean distance) will be fixed *a priori*.
        *   *Tool:* `dtwclust` or `tsclust` in R.
    2.  **(Alternative) Temporal Cross-Recurrence Plot (CRP):** Summarise diagonality from a binary matrix of pairwise distances < ε between encoding and retrieval fixations. The bandwidth ε will be fixed *a priori*.
    3.  **Combined Metric:** Fuse the spatial metric from Stage 1 with the temporal metric (soft-DTW distance) using a simple weighted sum.
*   **Evaluation:**
    *   Compare the Stage 1 + soft-DTW fusion against Stage 1 alone. Assess impact on behavioural correlations using mixed-effects models.
    *   **Quantitative Gate:** Proceed if Δβ (standardized, for vividness prediction) ≥ +0.10 or p < 0.01 (Holm-corrected) for the interaction term or main effect of the improved metric.
*   **Interpretability:** DTW alignment paths and CRPs offer visual insights.

**Stage 3: Graph-Theoretic Representation (Novel & Interpretable Features)**
*   *(Effort ≈ 1 month; RAM ≈ 30-40 GB for kernels on ~10k pairs; HPC/GPU access planned if needed)*
*   **Representation:** Convert each scanpath into an undirected weighted graph (Nodes = Fixation centroids; Edge weight = Saccade probability).
*   **Similarity Metrics/Features:**
    1.  Global Similarity (Graph Edit Distance, Delta-Con).
    2.  Graph Kernels (e.g., Weisfeiler-Lehman) fed to an SVM.
    3.  Extracted Interpretable Features (Degree centrality Gini, Edge-betweenness of AOIs, Triad census).
        *   *Tool:* `igraph` in R, `graphkernels` via reticulate.
*   **Evaluation:**
    *   Compare a model using EMD + DTW + extracted graph features against Stage 2.
    *   Test using 10-fold nested CV on same-image classification.
    *   Demonstrate unique predictive variance using hierarchical variance-partitioning (LRT), cross-validated ΔAUC decomposition, collinearity checks (VIF < 3), and feature importance.
    *   **Quantitative Gate:** Proceed to Stage 4 if Δχ² (LRT for graph features) shows p < 0.005.
*   **Interpretability:** Graph features are directly interpretable; kernels can be probed.

**Stage 4: Contrastive Metric Learning with a Siamese ScanPath Encoder (Powerful & New)**
*   *(Effort ≈ 1-2 months; RAM ≈ 30-40 GB for training on ~10k pairs; HPC/GPU access planned)*
*   **Representation:** Sequence of fixations (x, y, duration, velocity) → positional embeddings → small Transformer or 1D CNN.
*   **Loss Function:** InfoNCE / Triplet loss.
*   **Output:** Learned embedding space where cosine similarity reflects scanpath similarity.
    *   *Tool:* `torch` or `keras` in R.
*   **Decision Rule & Data Requirements:** "Siamese encoder proceeds if Stage 3 < 85 % of picture-ID ceiling (see Section 4.1) and sample size ≥ 3,000 positive pairs." Heavy data augmentation will be used.
*   **Evaluation:** Compare learned similarity against Stage 3.
*   **Interpretability Probes:** t-SNE plots, Gradient × Input, frozen backbone + GLM.

**Stage 5: Hybrid “Gaze-Conditioned Optimal Transport” (Blue-Sky Idea)**
*   *(Effort ≈ 1 month (if pursued); RAM dependent on implementation, likely >16GB)*
*   **Representation:** Fixation as 5-D point: (x, y, t, Δt, v).
*   **Similarity:** 5-D Wasserstein distance with optimized λ_t, λ_v per participant.
*   **Evaluation:** Ablation tests (λ_v=0 vs. optimized), variance-partitioning.

**4. Cross-Cutting Validation, Power Boost, and Decision Framework**

**4.1. Estimating the "Noise Ceiling" / Classification Head-Room**
To inform the 85% decision rule for Stage 4:

1.  **Split-Half Scanpath Reliability:** As previously described.
2.  **Oracle Classifiers:**
    *   **Picture-ID Oracle:** (Leave-picture-out CV).
    *   **Subject-ID Oracle:** (Leave-subject-out CV).
    *   **Joint Oracle:** (Random 80/20 split).
    *   *All three will be reported. The "effective" ceiling for cross-subject generalization will be considered the lower of the picture-ID oracle and the GBT proxy to avoid inflation by subject fingerprints.*
3.  **Bayes Error Proxy (GBT):** XGBoost on flattened fixation sequences (first 30 enc + 30 ret, NA padding, plus fixation counts).

**Rule Update for Stage 4:** The multi-scale-EMD + soft-DTW (Stage 2) or graph-based model (Stage 3) performance will be compared against 85% of the **Picture-ID Oracle ceiling** (or the GBT proxy if lower and deemed more appropriate for generalization).

**4.2. General Validation Tools**

*   **Permutation Framework:** As previously described (z-scores).
*   **Multi-Metric Fusion:** "Multi-metric fusion will be explored only after single-metric benchmarking to avoid curse-of-dimensionality." If pursued, logistic regression or random forest will be used.
*   **Reliability Check:** As previously described (split-half).
*   **Net Re-classification Improvement (NRI):** Pre-registered as a secondary criterion for model comparison, guarding against cases where rank order improves (AUC) but calibration worsens.

**4.3. Systematic Participant-Level Optimization & Heterogeneity**
"Participant-level optimisation: subject-specific σ/ω (and other relevant parameters like λ for Stage 5) become random slopes in mixed-effects models; inference reported at group level after partial pooling (e.g., via hierarchical Bayesian models if tailoring shows significant benefit via LRT comparison of random slope vs. fixed slope models)."
1.  Per-subject hyper-parameter grid optimization.
2.  Random-slopes mixed models.
3.  Report heterogeneity of optimal parameters and correlate with traits.
4.  Shrinkage estimators (hierarchical Bayesian models) if tailoring helps.

**4.4. Qualitative Error Triage (Embedded at Stage Ends)**
A systematic qualitative audit of "difficult" trials:
1.  Identify discrepant pairs (top false positives/negatives; high-vividness/low-similarity quartiles).
2.  Visual diagnostics (dual heat-maps, animated replays, AOI strings).
3.  Annotation worksheet for two coders (tagging factors like partial cues, central bias, idiosyncratic strategy).
4.  Feedback loop to inform metric tweaks.
A tiny **Shiny dashboard** will be developed early on (heat-map overlay, GIF replay, coding form) to streamline this process and reduce coder burden; this tool will be mentioned in the Methods. A "Qualitative Error Analysis" subsection will feature representative vignettes.

**5. Computational Logistics & Reproducibility**

*   **Tool Chain:** R (`transport`, `dtwclust`, `igraph`, `targets` for workflow management), Python via `reticulate` (`graphkernels`), `torch`/`keras` (for R or Python).
*   **Hardware Budget:** Explicit RAM and potential HPC/GPU needs noted per stage.
*   **Data Storage & Pre-computation:** Large O(n²) similarity matrices will be managed, potentially with HDF5 or `fst` for efficient caching.
*   **Data Versioning:** A data-versioning plan (e.g., DVC - Data Version Control, or Git-LFS for key matrices/CSVs) will be implemented to ensure pre-computed similarity matrices and intermediate datasets are traceable and shareable.
*   **Reproducibility:** A minimal Docker/Singularity recipe will be provided in an Appendix to facilitate turnkey reproducibility of the computational environment.

**6. Implementation Roadmap & Quantitative Gates**

| Month(s)      | Key Milestones                                                                    | Quantitative Gate / Action Item                                                                                                                                                                                                 |
| :------------ | :-------------------------------------------------------------------------------- | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| **1**         | Stage 0 (Baseline). Implement Stage 1 (EMD). Gaze QC.                            | ΔAUC (S1 vs S0) ≥ `max(0.02, 0.05*baseline_SE)` & 95% CI > 0. Compute oracles; update 85% rule. Power analysis for S1 gate.                                                                                                       |
| **1.5-2**     | Stage 1 Eval & Error Triage (Shiny app ready). Implement Stage 2 (soft-DTW).     | Δβ (vividness, S2 vs S1) ≥ +0.10 or p < 0.01 (Holm). Run GBT proxy.                                                                                                                                                                  |
| **2.5-3**     | Stage 2 Eval & Error Triage. Implement Stage 3 (Graph features). Define DTW cost/ε. | Δχ² (LRT for graph feat. vs S2) p < 0.005. VIF < 3. NRI check.                                                                                                                                                                    |
| **3.5-4**     | Stage 3 Eval & Error Triage. Assess participant-level optimization.               | LRT p < 0.01 favoring random slopes. **Decision Point:** Stage 4 if <85% of (Picture-ID oracle OR GBT proxy, whichever is lower & relevant for generalization) ceiling met AND sample size ≥ 3k positive pairs. NRI check. |
| **4+ (If needed)** | Prototype Stage 4 (Siamese). Stage 5 (5D OT) consideration.                      | Compare to S3. Assess remaining noise-ceiling gap. NRI check.                                                                                                                                                                   |
| **Ongoing**   | Qualitative Error Triage packets produced at each stage end.                       | Analyst notes and visual diagnostics via Shiny app.                                                                                                                                                                             |

**7. Conclusion & Expected Outcomes**

This systematic, multi-stage proposal provides a clear path to significantly enhance the analysis of eye-movement similarity. By combining incremental advancements with rigorous, quantitatively-gated evaluations, qualitative insights, and a commitment to open science, we aim to:

*   Develop more sensitive and nuanced measures of gaze reinstatement.
*   Understand the contributions of spatial, temporal, and dynamic gaze features.
*   Identify when complex models offer benefits over simpler ones.
*   Provide a robust framework adaptable for future research.

**Dissemination Plan:**
We commit to open science practices. This will include:
*   Pre-registration of the study plan (e.g., on OSF).
*   Making analysis code open-source (e.g., via GitHub).
*   Sharing anonymized data and pre-computed similarity matrices where ethically permissible and feasible.
*   Publishing findings in peer-reviewed journals and presenting at conferences, with a preprint made available (e.g., on bioRxiv/PsyArXiv).

Okay, I will add this explanation as an addendum to the main proposal.

## Revised Proposal: Advancing Eye-Movement Similarity Analysis for Encoding-Retrieval Studies

**(Main proposal content as previously detailed remains here – from Section 1. Introduction & Motivation through Section 7. Conclusion & Expected Outcomes)**

---

**Addendum: Explanation of "Signing" the Earth-Mover's Distance (EMD)**

This addendum clarifies the concept and implementation of "signing" the Earth-Mover's Distance, a technique optionally employed in Stage 1 of the proposed methodology to refine eye-movement similarity measures by accounting for bottom-up visual saliency.

**A.1. Plain-English Idea**

Eye-tracking heat maps are naturally all-positive probability distributions: every pixel holds non-negative fixation density, and the whole map sums to 1. The vanilla Earth-Mover Distance (EMD) between two such maps is always non-negative—it only tells you how much work it takes to morph one distribution into the other.

However, suppose both observers (e.g., during encoding and retrieval of the same image) spent a lot of time on regions that were already highly bottom-up salient (e.g., bright areas, high-contrast corners, faces). A significant portion of the measured EMD might simply reflect this shared "saliency gravity," rather than the idiosyncratic, memory-driven component of gaze that is often of primary interest in cognitive studies.

By subtracting a saliency baseline map first, each fixation heat map is converted into a *signed residual map*. In these residual maps:
*   Positive values indicate regions where the viewer looked *more* than predicted by bottom-up saliency.
*   Negative values indicate regions where the viewer looked *less* than predicted by bottom-up saliency.

Computing EMD on these residual maps, therefore, measures how similarly the two scanpaths *depart from* the baseline saliency predictions. The resulting score can be positive or negative (hence "signed EMD"), reflecting the nature and similarity of these departures.

**A.2. Step-by-Step Recipe for Signed EMD**

1.  **Compute a Saliency Map for Each Stimulus:**
    *   `S(x,y)` ← Prediction from a pre-specified saliency model (e.g., DeepGaze III, ITTI & Koch, or another chosen *a priori* model), normalized to integrate to 1. The choice of saliency model will be documented, and its potential scene-dependent failures noted as an area for monitoring during the qualitative error-triage.

2.  **Generate Standard Fixation-Density Maps (FDMs):**
    *   `F_enc(x,y)` (for encoding) and `F_ret(x,y)` (for retrieval), each smoothed with the same kernel (as defined in Stage 1) and normalized to integrate to 1.

3.  **Form Signed Residual Maps:**
    *   `R_enc = F_enc - S`
    *   `R_ret = F_ret - S`
    *   These maps now represent the deviation of observed gaze from saliency predictions. The sum of pixel values in a residual map is zero.

4.  **Split Residuals into Positive and Negative "Piles of Earth":**
    For each residual map (e.g., `R_enc`):
    *   `R_enc_positive(x,y) = max(R_enc(x,y), 0)`
    *   `R_enc_negative(x,y) = max(-R_enc(x,y), 0)`
    And similarly for `R_ret`.
    *   Each of these `R_positive` and `R_negative` maps now forms a proper, non-negative distribution. The total mass (sum of pixel values) in `R_positive` will equal the total mass in `R_negative` for a given original residual map (because the positive and negative parts of the residual map `R` must cancel each other out as `R` sums to zero).

5.  **Compute Two EMDs:**
    *   `EMD_positive` = EMD between `R_enc_positive` and `R_ret_positive`.
    *   `EMD_negative` = EMD between `R_enc_negative` and `R_ret_negative`.
    *   These measure the "work" required to transform the "looked more than saliency" patterns into each other, and similarly for the "looked less than saliency" patterns.

6.  **Derive a Single Signed Similarity Score:**
    *   `SignedEMD = -(EMD_positive + EMD_negative)`
    *   **Interpretation:**
        *   Small (more negative) values indicate a large total EMD cost, meaning the two observers departed from saliency in dissimilar ways.
        *   Large (less negative or even positive, approaching zero) values indicate a low total EMD cost, meaning the two observers over-shot and under-shot saliency in similar locations, suggesting genuine reinstatement of gaze patterns beyond bottom-up pull.
    *   *(Note: The sign convention is arbitrary as long as it's consistent. Some researchers might report the negative of the total cost so that "higher = more similar." The chosen convention will be clearly stated.)*

**A.3. Why Bother with Signed EMD?**

*   **Controls for Bottom-Up Bias:** If both scanpaths align with highly salient features simply because those features attract everyone's gaze, subtracting the saliency map neutralizes this shared variance, allowing for a cleaner measure of top-down influences.
*   **Highlights Mnemonic or Task-Driven Guidance:** The residual gaze patterns (regions looked at more or less than predicted by saliency) are stronger candidates for reflecting top-down, memory-driven, or task-relevant fixation choices.
*   **Produces an Effect Direction (Potentially):** While the sum of two non-negative EMDs is always non-negative, the specific formulation and interpretation focus on the *similarity of deviations*. The term "signed" primarily refers to the residual maps themselves. The final similarity score, as defined above, will range from highly negative (very dissimilar deviations) towards zero (very similar deviations).

**A.4. Quick R Pseudo-Code Example**

```R
# Assume F_enc, F_ret, and S are matrices (e.g., from imager::as.cimg())
# representing fixation density maps and the saliency map,
# all appropriately smoothed and normalized to sum to 1.

# library(imager) # for image handling, pmax
# library(transport) # for wasserstein()

# 1. Saliency map S is pre-computed
# 2. F_enc and F_ret are generated

# 3. Form signed residual maps
res_enc <- F_enc - S
res_ret <- F_ret - S

# 4. Split residuals into positive and negative "piles"
pos_enc <- pmax(res_enc, 0)
neg_enc <- pmax(-res_enc, 0) # Note: -res_enc, so positive values here represent where original was negative

pos_ret <- pmax(res_ret, 0)
neg_ret <- pmax(-res_ret, 0)

# Ensure these are suitable for transport::wasserstein (e.g., as 2D histograms or point patterns)
# This might require converting matrix representations to the format expected by the function.
# For 2D histograms, transport::wasserstein() expects inputs that are themselves
# representations of distributions (e.g., objects from hist() or similar structures).
# If F_enc etc. are simple matrices, they might need to be converted.
# For simplicity, let's assume they can be directly used or easily converted.

# 5. Compute two EMDs (using p=1 for Wasserstein-1 distance)
# The exact call to wasserstein() will depend on how pos_enc etc. are structured.
# Assuming they are compatible 2D histograms / distributions:
emd_pos <- transport::wasserstein(pp(pos_enc), pp(pos_ret), p = 1) # pp() might be a placeholder for conversion
emd_neg <- transport::wasserstein(pp(neg_enc), pp(neg_ret), p = 1) # to point pattern or suitable histogram format

# 6. Derive a single signed similarity score
signed_similarity_score <- -(emd_pos + emd_neg) # Higher values (closer to 0) mean more similar deviations

# This signed_similarity_score is then used in subsequent analyses
# (permutation framework, mixed models, etc.) in place of the vanilla EMD.
```

*(The `pp()` in the pseudo-code is a placeholder; actual conversion to the input format required by `transport::wasserstein` for 2D histograms would be needed, e.g., ensuring they are structured as outputs from `hist` or as two-column matrices of coordinates and masses for point patterns.)*

**A.5. Take-Home Message**

"Signing the EMD" is a methodological refinement that aims to isolate the cognitive component of gaze similarity. By first correcting each fixation map for baseline, bottom-up visual saliency, the subsequent EMD calculation reflects how similarly two gaze patterns deviate from these predictions. This approach provides a more targeted measure of memory-driven or task-relevant gaze reinstatement, rather than simply quantifying the overlap in raw saliency-influenced gaze distributions. The resulting `signed_similarity_score` is then used within the same statistical framework as other similarity metrics.


Okay, I will add a second addendum detailing how to wrap DeepGaze in R via `reticulate`, including a helper function for installation. This assumes a Python environment with DeepGaze III (or a similar version) is accessible or can be set up.

**(Main proposal content as previously detailed remains here – from Section 1. Introduction & Motivation through Section 7. Conclusion & Expected Outcomes, followed by Addendum A: Explanation of "Signing" the Earth-Mover's Distance (EMD))**

---

**Addendum B: Integrating DeepGaze Saliency Model in R via Reticulate**

This addendum provides a practical guide for using a Python-based saliency model like DeepGaze III within an R workflow using the `reticulate` package. This is relevant for Stage 1 when optionally subtracting a per-image saliency map to compute Signed EMD.

**B.1. Overview**

DeepGaze models are typically implemented in Python using deep learning frameworks like PyTorch or TensorFlow. The `reticulate` package in R allows for seamless interoperability between R and Python, enabling us to call Python functions and use Python objects directly from R.

**B.2. Prerequisites**

1.  **Python Installation:** A working Python installation (e.g., Python 3.7+) is required.
2.  **R and Reticulate:** R and the `reticulate` package must be installed (`install.packages("reticulate")`).
3.  **Python Environment for DeepGaze:** It's highly recommended to create a dedicated Python virtual environment (e.g., using `conda` or `venv`) to manage DeepGaze and its dependencies without conflicting with other Python projects.
4.  **DeepGaze Installation:** The chosen DeepGaze implementation (e.g., DeepGaze III from a specific GitHub repository or package) needs to be installed into this Python environment along with its dependencies (PyTorch, NumPy, PIL/Pillow, etc.).

**B.3. Helper Function: `install_deepgaze_env()`**

This R function aims to simplify the setup of a Conda environment with DeepGaze. *Users should adapt paths and package names based on the specific DeepGaze version and their system.* This example assumes a hypothetical `deepgaze_pytorch` package available via pip and its dependencies.

```R
#' Install a Conda Environment for DeepGaze
#'
#' This function attempts to create a new Conda environment named 'r_deepgaze_env'
#' and install necessary Python packages for running a PyTorch-based DeepGaze model.
#'
#' @param envname Character string, the name of the Conda environment to create.
#'                Defaults to "r_deepgaze_env".
#' @param python_version Character string, the Python version for the environment.
#'                       Defaults to "3.9".
#' @param pytorch_cpu_only Logical, if TRUE, installs the CPU-only version of PyTorch.
#'                         Defaults to TRUE for wider compatibility.
#'                         Set to FALSE if a GPU version is desired and CUDA is configured.
#' @param deepgaze_package_name Character string, the pip installable name for DeepGaze.
#'                              This is HYPOTHETICAL and needs to be replaced with the actual
#'                              package name or installation command (e.g., from a GitHub repo).
#'                              Defaults to "deepgaze_pytorch".
#'
#' @details
#' Users might need to adjust 'deepgaze_package_name' and potentially other
#' dependencies based on the specific DeepGaze implementation they intend to use.
#' This function requires Conda to be installed and accessible in the system PATH.
#'
#' It's often more robust to manually create the Python environment and install
#' DeepGaze following its official documentation, then point reticulate to it.
#' This function is provided as a convenience but might require troubleshooting.
#'
#' @examples
#' \dontrun{
#'   # Before running, ensure 'deepgaze_package_name' is correct for your DeepGaze source
#'   # install_deepgaze_env(deepgaze_package_name = "pip_package_for_my_deepgaze_version")
#'   # OR, if installing from GitHub:
#'   # install_deepgaze_env(deepgaze_package_name = "git+https://github.com/user/deepgaze-repo.git")
#' }
#'
install_deepgaze_env <- function(envname = "r_deepgaze_env",
                                 python_version = "3.9",
                                 pytorch_cpu_only = TRUE,
                                 deepgaze_package_name = "deepgaze_pytorch") { # HYPOTHETICAL

  if (!reticulate::conda_is_installed()) {
    message("Conda is not installed or not found by reticulate.")
    message("Please install Miniconda/Anaconda and ensure it's in your PATH,")
    message("or use reticulate::install_miniconda().")
    return(invisible(NULL))
  }

  existing_envs <- reticulate::conda_list()$name
  if (envname %in% existing_envs) {
    message(paste("Conda environment '", envname, "' already exists.", sep = ""))
    message("To use it: reticulate::use_condaenv('", envname, "', required = TRUE)", sep = "")
    return(invisible(NULL))
  }

  message(paste("Attempting to create Conda environment: '", envname, "'...", sep = ""))
  tryCatch({
    reticulate::conda_create(envname = envname, python_version = python_version)

    packages_to_install <- c("numpy", "Pillow") # Basic dependencies

    # PyTorch installation command varies
    if (pytorch_cpu_only) {
      # Example for CPU-only PyTorch from PyTorch.org (check for current command)
      # This command can be highly specific to OS and PyTorch version.
      # Using pip install within conda env for simplicity here.
      # A more robust way is to use conda install pytorch torchvision torchaudio cpuonly -c pytorch
      packages_to_install <- c(packages_to_install, "torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu")
    } else {
      # GPU PyTorch installation is more complex and system-dependent
      # (e.g., "conda install pytorch torchvision torchaudio cudatoolkit=11.3 -c pytorch")
      # For this example, we'll stick to a generic pip install command
      # that the user would need to ensure installs the GPU version.
      packages_to_install <- c(packages_to_install, "torch torchvision torchaudio")
      message("Note: For GPU PyTorch, ensure your CUDA drivers and toolkit are correctly set up.")
    }

    # Add the DeepGaze package itself (THIS IS HYPOTHETICAL)
    # The user MUST replace 'deepgaze_package_name' with the actual pip install command
    # or instructions for their chosen DeepGaze version.
    # It might be 'pip install deepgaze3', 'pip install git+https://github.com/user/repo.git', etc.
    if (startsWith(deepgaze_package_name, "git+")) {
         packages_to_install <- c(packages_to_install, deepgaze_package_name)
    } else {
         packages_to_install <- c(packages_to_install, deepgaze_package_name)
    }


    message("Installing packages into '", envname, "': ", paste(packages_to_install, collapse=", "), "...")
    # Install using pip within the conda environment
    # Note: reticulate::conda_install can be tricky with complex packages like PyTorch.
    # Using pip directly via conda_run is often more reliable for these.
    for (pkg in strsplit(packages_to_install, " ")[[1]]) { # simplistic split
        if (grepl("torch", pkg) && grepl("http", pkg)) { # handle complex torch install command
             reticulate::conda_run(paste("pip install", paste(strsplit(packages_to_install, " ")[[1]][which(grepl("torch", strsplit(packages_to_install, " ")[[1]])):length(strsplit(packages_to_install, " ")[[1]])], collapse=" ")), envname = envname)
             break # Assuming torch command is last
        } else if (!grepl("http", pkg) && !grepl("git+", pkg)) {
             reticulate::conda_install(envname = envname, packages = pkg, pip = TRUE)
        } else if (grepl("git+", pkg)) {
             reticulate::conda_run(paste("pip install", pkg), envname = envname)
        }
    }


    message(paste("Conda environment '", envname, "' created and packages installed (hopefully!).", sep = ""))
    message("To use this environment in your R session, run:")
    message("reticulate::use_condaenv('", envname, "', required = TRUE)", sep = "")
    message("Then import the DeepGaze module using reticulate::import('deepgaze_module_name').")

  }, error = function(e) {
    message(paste("Error creating Conda environment '", envname, "':", sep = ""))
    message(e)
    message("Please try creating the environment and installing packages manually.")
  })
  return(invisible(NULL))
}
```
**Important Considerations for `install_deepgaze_env()`:**
*   **Hypothetical Package:** `deepgaze_pytorch` is a placeholder. The actual installation command for DeepGaze (e.g., from PyPI or GitHub) must be used.
*   **PyTorch Installation:** PyTorch installation can be tricky and platform-dependent, especially for GPU versions. It's often best to follow official PyTorch installation instructions for Conda directly. The function provides a basic CPU-only pip attempt.
*   **Manual Setup Recommended:** For complex dependencies, manually setting up the Conda environment as per the DeepGaze model's documentation and then pointing `reticulate` to it (using `reticulate::use_condaenv()`) is often more reliable.

**B.4. Wrapping DeepGaze for Use in R**

Once the Python environment is set up and DeepGaze is installed, you can use `reticulate` to call it from R.

```R
# 0. Ensure the correct Conda environment is activated for reticulate
# This should be done ONCE per R session, ideally at the beginning.
# Replace 'r_deepgaze_env' with the actual name of your Conda environment.
tryCatch({
  reticulate::use_condaenv("r_deepgaze_env", required = TRUE)
  message("Successfully using Conda environment: r_deepgaze_env")
}, error = function(e) {
  message("Could not activate Conda environment 'r_deepgaze_env'.")
  message("Ensure it's created, DeepGaze is installed in it, and the name is correct.")
  message("You might need to run install_deepgaze_env() or set it up manually.")
  stop(e)
})


# 1. Import the necessary Python modules
# The exact module name for DeepGaze will depend on its Python implementation.
# This is HYPOTHETICAL.
# For example, if DeepGaze III is in a module called 'deepgaze.deepgaze_pytorch':
deepgaze_module <- NULL
tryCatch({
  # This name 'deepgaze_models' or similar is HYPOTHETICAL
  # You need to find the actual importable module name from your DeepGaze installation
  deepgaze_module <- reticulate::import("deepgaze.pytorch_model_for_example", delay_load = TRUE) # Replace with actual module
  # Or, if it's a specific class within a module:
  # DeepGazeModel <- reticulate::import_main()$MyDeepGazeClass # If in __main__ after a script run
  # Or:
  # dg_models <- reticulate::import("some_deepgaze_library.models")
  # DeepGazeModel <- dg_models$DeepGazeIII
  message("DeepGaze Python module imported (or load delayed).")
}, error = function(e) {
  message("Failed to import the DeepGaze Python module.")
  message("Ensure the module name is correct and it's installed in the active Conda environment.")
  stop(e)
})

# 2. Load the pre-trained DeepGaze model (Python side)
# This step is highly dependent on the specific DeepGaze API.
# It might involve specifying model paths, device (CPU/GPU), etc.
deepgaze_model_instance <- NULL # Global variable to hold the model instance

# Helper function to initialize the model (Python side) ONCE
initialize_deepgaze_model <- function(model_weights_path = "path/to/deepgaze_weights.pth",
                                      device = "cpu") { # "cuda" for GPU
  if (!is.null(deepgaze_model_instance) && inherits(deepgaze_model_instance, "python.builtin.object")) {
    message("DeepGaze model already initialized.")
    return(deepgaze_model_instance)
  }

  message("Initializing DeepGaze model instance...")
  # This is HYPOTHETICAL - adapt to your DeepGaze model's API
  tryCatch({
    # Example: Assuming the imported module has a function like `load_model`
    # or a class constructor.
    # deepgaze_model_instance <<- deepgaze_module$load_pretrained_model(weights_path = model_weights_path, device = device)
    # OR for a class:
    # deepgaze_model_instance <<- deepgaze_module$DeepGazeIII(pretrained=TRUE, device=device) # Hypothetical
    
    # --- REPLACE THIS BLOCK WITH ACTUAL DEEPGAZE MODEL LOADING CODE ---
    # For demonstration, let's assume a placeholder function exists in the Python module
    # You would replace 'PlaceholderDeepGazeModel' and its arguments.
    # Ensure your Python 'deepgaze.pytorch_model_for_example' has this.
    if("PlaceholderDeepGazeModel" %in% names(deepgaze_module)) {
        deepgaze_model_instance <<- deepgaze_module$PlaceholderDeepGazeModel(device = device)
        deepgaze_model_instance$eval() # Typically set model to evaluation mode
        message("DeepGaze model instance initialized and set to eval mode.")
    } else {
        warning("PlaceholderDeepGazeModel not found in the imported module. Model not loaded.")
        message("You MUST replace the model loading logic with your actual DeepGaze API.")
    }
    # --- END OF REPLACEABLE BLOCK ---

  }, error = function(e) {
    message("Error initializing DeepGaze model instance in Python:")
    stop(e)
  })
  return(deepgaze_model_instance)
}


# 3. R function to get saliency map for an image
#' Get Saliency Map using DeepGaze
#'
#' Loads an image, preprocesses it, gets saliency prediction from DeepGaze,
#' and returns the saliency map as an R matrix or imager::cimg object.
#'
#' @param image_path Character string, path to the input image file.
#' @param model_instance A Python object representing the loaded DeepGaze model.
#' @param target_size Numeric vector c(width, height), e.g., c(1024, 768),
#'                    to resize the saliency map to match FDMs.
#'                    If NULL, original saliency map size is returned.
#' @return A matrix or imager::cimg object representing the saliency map,
#'         normalized to sum to 1. Returns NULL on error.
#'
get_saliency_deepgaze <- function(image_path, model_instance, target_size = NULL) {
  if (is.null(model_instance) || !inherits(model_instance, "python.builtin.object")) {
    message("DeepGaze model instance is not valid. Initialize it first.")
    return(NULL)
  }
  if (!file.exists(image_path)) {
    message("Image file not found: ", image_path)
    return(NULL)
  }

  PIL <- reticulate::import("PIL.Image", delay_load = TRUE)
  transforms <- reticulate::import("torchvision.transforms", delay_load = TRUE) # If needed for preprocessing
  torch <- reticulate::import("torch", delay_load = TRUE) # If manual tensor ops needed

  tryCatch({
    # Load image using PIL (Python Imaging Library)
    img_pil <- PIL$open(image_path)$convert('RGB')

    # --- REPLACE THIS BLOCK WITH ACTUAL DEEPGAZE PREPROCESSING AND PREDICTION ---
    # Preprocessing: This is highly model-specific.
    # Example: typical PyTorch preprocessing
    # preprocess <- transforms$Compose(list(
    #   transforms$Resize(list(224L, 224L)), # Example size
    #   transforms$ToTensor(),
    #   transforms$Normalize(mean=c(0.485, 0.456, 0.406), std=c(0.229, 0.224, 0.225))
    # ))
    # img_tensor <- preprocess(img_pil)$unsqueeze(0L) # Add batch dimension

    # If your model takes a raw PIL image or numpy array, adapt accordingly.
    # For our placeholder, let's assume it takes a PIL image and returns a tensor.
    # This part is CRITICALLY dependent on your specific DeepGaze model's API.
    
    # Get prediction (output is usually a PyTorch tensor)
    # saliency_tensor_pred <- model_instance$predict(img_tensor) # Hypothetical predict method
    
    # Using our placeholder method that might exist in the Python module
    # Ensure your Python 'deepgaze.pytorch_model_for_example' has this.
    if("predict_saliency_from_pil" %in% names(model_instance)) {
        saliency_tensor_pred <- model_instance$predict_saliency_from_pil(img_pil)
    } else {
        warning("predict_saliency_from_pil method not found on model instance. Cannot predict.")
        message("You MUST replace prediction logic with your actual DeepGaze API.")
        return(NULL)
    }
    # --- END OF REPLACEABLE BLOCK ---


    # Postprocess: Convert tensor to R matrix, resize, normalize
    # Squeeze batch and channel dimensions if necessary (common for saliency maps)
    saliency_map_py <- saliency_tensor_pred$squeeze()$cpu()$detach()$numpy() # To NumPy array

    # Convert to R matrix
    saliency_map_r <- as.matrix(saliency_map_py)

    # Resize if target_size is provided (using imager for convenience)
    if (!is.null(target_size) && requireNamespace("imager", quietly = TRUE)) {
      saliency_cimg <- imager::as.cimg(saliency_map_r)
      # imager's resize uses x,y,z,c dimensions. Saliency is usually 2D.
      # Ensure aspect ratio is handled as desired (e.g., stretch vs. pad)
      saliency_cimg_resized <- imager::resize(saliency_cimg,
                                              size_x = target_size[1],
                                              size_y = target_size[2],
                                              interpolation_type = 6) # Lanczos for quality
      saliency_map_r <- as.matrix(saliency_cimg_resized)
    }

    # Normalize to sum to 1
    saliency_map_r <- saliency_map_r - min(saliency_map_r) # Ensure non-negative
    if (sum(saliency_map_r) > 1e-6) {
      saliency_map_r <- saliency_map_r / sum(saliency_map_r)
    } else {
      # Handle case of all-zero map (e.g., return uniform distribution)
      saliency_map_r[,] <- 1 / (nrow(saliency_map_r) * ncol(saliency_map_r))
      warning("Saliency map was near zero; returning uniform distribution.")
    }

    return(saliency_map_r)

  }, error = function(e) {
    message(paste("Error getting saliency map for", image_path, ":"))
    message(e)
    return(NULL)
  })
}

# --- Example Usage (Illustrative - requires actual DeepGaze setup) ---
#
# # A. One-time setup (or at start of script/project):
# # 1. Ensure you have a Conda env named "r_deepgaze_env" with Python, PyTorch,
# #    Pillow, NumPy, and your specific DeepGaze package installed.
# #    You might use install_deepgaze_env() CAREFULLY or set it up manually.
# #
# # install_deepgaze_env(deepgaze_package_name = "actual_pip_name_or_git_url_for_deepgaze")
#
# # 2. Activate the environment (already done at the top of this script section)
# # reticulate::use_condaenv("r_deepgaze_env", required = TRUE)
#
# # 3. Initialize the DeepGaze model (Python side)
# #    Replace "path/to/deepgaze_weights.pth" with the actual path if your model
# #    requires loading weights explicitly. Some models load them automatically if pretrained.
# #    The Python module `deepgaze.pytorch_model_for_example` and the class/function
# #    `PlaceholderDeepGazeModel` and its method `predict_saliency_from_pil`
# #    ARE HYPOTHETICAL and need to be replaced with your actual DeepGaze model's API.
#
# # dg_model <- initialize_deepgaze_model(device = "cpu") # or "cuda"
#
# # B. For each image:
# # if (!is.null(dg_model)) {
# #   # Define target dimensions for the saliency map, e.g., matching your screen/FDM resolution
# #   screen_width <- 1024
# #   screen_height <- 768
# #   saliency_map <- get_saliency_deepgaze(
# #     image_path = "path/to/your/stimulus_image.jpg",
# #     model_instance = dg_model,
# #     target_size = c(screen_width, screen_height)
# #   )
# #
# #   if (!is.null(saliency_map)) {
# #     # Now 'saliency_map' is an R matrix, normalized, and ready for use
# #     # in the Signed EMD calculation (see Addendum A)
# #     print(paste("Saliency map dimensions:", paste(dim(saliency_map), collapse="x")))
# #     print(paste("Sum of saliency map:", sum(saliency_map))) # Should be 1
# #
# #     if (requireNamespace("imager", quietly = TRUE) && interactive()) {
# #        # Quick plot if imager is available
# #        # plot(imager::as.cimg(saliency_map), main="DeepGaze Saliency Map")
# #     }
# #   }
# # } else {
# #   message("DeepGaze model not initialized. Cannot generate saliency map.")
# # }
```

**B.5. Crucial Considerations and Customization:**

*   **DeepGaze API Specifics:** The core of `initialize_deepgaze_model` and `get_saliency_deepgaze` (especially image preprocessing, model prediction call, and output handling) **MUST be adapted** to the precise API of the DeepGaze version being used. The provided code uses placeholders (e.g., `deepgaze.pytorch_model_for_example`, `PlaceholderDeepGazeModel`, `predict_saliency_from_pil`). These will not work out-of-the-box and require replacement with actual, functional Python code corresponding to the chosen DeepGaze library.
*   **Python Module Name:** The `reticulate::import("module_name")` call needs the correct Python module name where the DeepGaze model is defined.
*   **Model Loading:** How pre-trained weights are loaded is model-specific. Some models load them automatically with `pretrained=True`, others require a path to a weights file.
*   **Preprocessing:** Image preprocessing steps (resizing, normalization, tensor conversion) must match what the DeepGaze model expects during its training.
*   **Output Format:** Saliency models output data in various formats (e.g., PyTorch tensors, NumPy arrays). The R code needs to handle this and convert it to a suitable R matrix.
*   **Device Management (CPU/GPU):** If using a GPU, ensure PyTorch (or TensorFlow) is installed with GPU support in the Conda environment and that `device="cuda"` (or similar) is correctly passed to the model and tensors.
*   **Error Handling:** Robust error handling is important, especially around Python interactions.
*   **Efficiency:** For processing many images, initializing the Python model once and reusing the instance (`deepgaze_model_instance`) is crucial for performance, as model loading can be slow.

This addendum provides a structural template. Successful integration will require careful attention to the documentation and API of the specific DeepGaze implementation chosen for the project. It's advisable to first ensure the DeepGaze model runs correctly in a pure Python script within the target Conda environment before attempting to wrap it with `reticulate`.