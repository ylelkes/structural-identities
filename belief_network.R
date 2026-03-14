# ============================================================
# belief_network.R
#
# Constructs and visualises a weighted belief network from the
# pairwise association task (Block 2 only — Block 1 removed).
#
# RT-weighted edge strength: WITHIN-PERSON normalized.
# For each participant, rt_ratio = pair_median_rt / participant_median_rt
#   (ratio=1 = at that person's average; <1 = faster; >1 = slower)
# Then across participants: mean_ratio = mean(rt_ratio per pair)
#   rt_factor          = max(0.4, 1 − 0.3 × mean_ratio)
#   rt_weighted_strength = assoc_strength × rt_factor
#
# Scaling:   ratio=0.5 → 0.85  |  ratio=1.0 → 0.70  |  ratio≥2 → 0.40
#
# INPUT (from simulate.html or real experiment — place in simulated_data/):
#   block1_evaluation_summary.csv  — B1 good/bad responses, used for node valence colour
#   block2_pair_summary.csv        — B2 pair associations (must contain rt_ratio column)
#   (or all_trials.csv as fallback for B2 if the pair summary is missing)
#
# OUTPUTS (written to plots/):
#   01_full_network.pdf         main belief network, full sample
#   02_rt_scatter.pdf           RT vs association strength scatter
#   03_association_heatmap.pdf  pair-level association matrix
#   04_group_networks.pdf       liberal vs centrist vs conservative
#   05_me_connections.pdf       self-concept associations
# ============================================================

# ── 0. Packages ──────────────────────────────────────────────────────────────

required <- c("tidyverse", "igraph", "ggraph", "tidygraph",
              "graphlayouts", "scales", "patchwork", "ggrepel")

for (pkg in required) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  library(pkg, character.only = TRUE)
}

# ── 1. Paths & directories ───────────────────────────────────────────────────

script_dir <- if (interactive()) {
  dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  normalizePath(".")
}

PLOT_DIR <- file.path(script_dir, "plots")
dir.create(PLOT_DIR, showWarnings = FALSE)

find_file <- function(filename) {
  candidates <- c(
    file.path(script_dir, "simulated_data", filename),
    file.path(script_dir, filename)
  )
  found <- candidates[file.exists(candidates)]
  if (length(found)) found[1] else stop("Cannot find ", filename,
    "\nLooked in:\n", paste(candidates, collapse = "\n"))
}

# ── 2. Load data ─────────────────────────────────────────────────────────────

# Block 1 — evaluation summary (used for node valence coloring)
b1_path <- tryCatch(find_file("block1_evaluation_summary.csv"), error = function(e) NULL)
b1_raw  <- if (!is.null(b1_path)) {
  message("✓ Loaded block1_evaluation_summary.csv")
  read_csv(b1_path, show_col_types = FALSE)
} else {
  message("⚠ block1_evaluation_summary.csv not found — nodes will use category colours")
  NULL
}

b2_path <- tryCatch(find_file("block2_pair_summary.csv"), error = function(e) NULL)

if (!is.null(b2_path)) {
  b2_raw <- read_csv(b2_path, show_col_types = FALSE)
  message("✓ Loaded block2_pair_summary.csv")

} else {
  message("⚙ Summary file not found — computing from all_trials.csv")
  trials_path <- find_file("all_trials.csv")
  trials      <- read_csv(trials_path, show_col_types = FALSE)

  b2_trials <- trials %>%
    filter(task_type == "pairwise") %>%
    mutate(
      n1   = pmin(left_node_label,  right_node_label),
      n2   = pmax(left_node_label,  right_node_label),
      cat1 = if_else(left_node_label  == n1, left_node_category,  right_node_category),
      cat2 = if_else(left_node_label  == n1, right_node_category, left_node_category),
      pid  = paste(n1, n2, sep = " | ")
    )

  # Within-person median RT baseline (across all non-timeout B2 trials)
  participant_median_rt <- b2_trials %>%
    filter(timeout == 0) %>%
    group_by(participant_id) %>%
    summarise(participant_median_rt = median(rt_ms, na.rm = TRUE), .groups = "drop")

  b2_raw <- b2_trials %>%
    group_by(participant_id, participant_type,
             pair_id = pid,
             node1_label = n1, node2_label = n2,
             node1_category = cat1, node2_category = cat2,
             mapping_condition) %>%
    summarise(
      n_trials           = n(),
      n_together         = sum(response_label == "together", na.rm = TRUE),
      n_apart            = sum(response_label == "apart",    na.rm = TRUE),
      prop_together      = n_together / n_trials,
      prop_apart         = n_apart    / n_trials,
      median_rt_all      = median(rt_ms[timeout == 0], na.rm = TRUE),
      median_rt_together = median(rt_ms[response_label == "together" & timeout == 0],
                                  na.rm = TRUE),
      median_rt_apart    = median(rt_ms[response_label == "apart"    & timeout == 0],
                                  na.rm = TRUE),
      timeout_rate       = mean(timeout, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(participant_median_rt, by = "participant_id") %>%
    mutate(rt_ratio = if_else(
      !is.na(median_rt_all) & participant_median_rt > 0,
      median_rt_all / participant_median_rt,
      NA_real_
    ))

  message("✓ B2 summary computed from all_trials.csv")
}

N_participants <- n_distinct(b2_raw$participant_id)

# ── 3. Configuration ─────────────────────────────────────────────────────────

# Edge threshold: only edges where |mean_assoc − 0.5| exceeds this
EDGE_THRESHOLD <- 0.05

# Minimum participants who must have seen a pair for it to be included
MIN_PARTICIPANTS_PER_PAIR <- 3

FIRST_NAME_DISPLAY <- "me"

# Category colours (matches simulate.html CAT_COLORS)
CATEGORY_COLORS <- c(
  identity = "#5b8dd9",
  policy   = "#6ac46a",
  values   = "#c97de8"
)

# ── 4. Node list (derived from B2 pairs — no Block 1 needed) ─────────────────

node_summary <- bind_rows(
  b2_raw %>% transmute(node = node1_label, category = node1_category),
  b2_raw %>% transmute(node = node2_label, category = node2_category)
) %>%
  mutate(node = if_else(node == "first_name", FIRST_NAME_DISPLAY, node)) %>%
  distinct(node, category)

# ── 4b. Node valence from Block 1 (RT-weighted) ──────────────────────────────
#
# For each node, compute:
#   mean_ratio_b1 = mean(rt_ratio_b1) across participants
#   rt_factor_b1  = max(0.4, 1 − 0.3 × mean_ratio_b1)  [within-person normalized]
#   valence_centered = (mean_prop_good − 0.5) × rt_factor_b1
#   valence_score    = 0.5 + valence_centered   [0=bad, 0.5=neutral, 1=good]

if (!is.null(b1_raw)) {
  # Compute within-person B1 RT baseline if not already in the file
  if (!"rt_ratio" %in% names(b1_raw)) {
    b1_raw <- b1_raw %>%
      group_by(participant_id) %>%
      mutate(participant_median_rt_b1 = median(median_rt_all, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(rt_ratio = if_else(
        !is.na(median_rt_all) & participant_median_rt_b1 > 0,
        median_rt_all / participant_median_rt_b1,
        NA_real_
      ))
  }

  node_valence <- b1_raw %>%
    mutate(
      node = if_else(target_label == "first_name", FIRST_NAME_DISPLAY, target_label)
    ) %>%
    group_by(node) %>%
    summarise(
      mean_prop_good  = mean(prop_good,  na.rm = TRUE),
      mean_ratio_b1   = mean(rt_ratio,   na.rm = TRUE),
      n_participants  = n_distinct(participant_id),
      .groups = "drop"
    ) %>%
    mutate(
      rt_factor_b1      = pmax(0.4, 1 - 0.3 * coalesce(mean_ratio_b1, 1.0)),
      valence_centered  = (mean_prop_good - 0.5) * rt_factor_b1,
      valence_score     = 0.5 + valence_centered   # [~0.2, ~0.8] in practice
    )

  node_summary <- node_summary %>%
    left_join(node_valence %>% select(node, mean_prop_good, valence_score),
              by = "node")

  message(sprintf("✓ Valence computed for %d nodes", sum(!is.na(node_summary$valence_score))))
} else {
  node_summary <- node_summary %>%
    mutate(mean_prop_good = NA_real_, valence_score = NA_real_)
}

# Colour ramp: dark red (bad) → slate neutral → dark blue (good)
valence_to_color <- function(score) {
  neutral  <- col2rgb("#3a4455") / 255
  bad_col  <- col2rgb("#8b0000") / 255
  good_col <- col2rgb("#003580") / 255
  ifelse(
    is.na(score), "#3a4455",
    ifelse(
      score < 0.5,
      rgb(neutral + (bad_col  - neutral) * (1 - score * 2)),
      rgb(neutral + (good_col - neutral) * ((score - 0.5) * 2))
    )
  )
}

node_summary <- node_summary %>%
  mutate(node_fill = valence_to_color(valence_score))

# ── 5. Edge-level aggregation with RT weighting ───────────────────────────────
#
# Within-person RT normalization (mirrors simulate.html aggB2):
#   mean_ratio = mean(rt_ratio) across participants  [ratio = pairRT / personMeanRT]
#   rt_factor  = max(0.4, 1 − 0.3 × mean_ratio)
#
# rt_weighted_strength = assoc_strength × rt_factor
# Both "go together" and "do not go together" fast responses count —
# strength = absolute deviation from 0.5, rt_factor scales it by within-person RT speed.

# Verify rt_ratio column exists (present in CSV from updated simulate.html)
if (!"rt_ratio" %in% names(b2_raw)) {
  warning("rt_ratio column not found in block2_pair_summary.csv — ",
          "falling back to participant-level calculation. ",
          "Re-run simulate.html to generate updated CSV.")
  b2_raw <- b2_raw %>%
    group_by(participant_id) %>%
    mutate(participant_median_rt = median(median_rt_all, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(rt_ratio = if_else(
      !is.na(median_rt_all) & participant_median_rt > 0,
      median_rt_all / participant_median_rt,
      NA_real_
    ))
}

edge_summary <- b2_raw %>%
  mutate(
    node1_label = if_else(node1_label == "first_name", FIRST_NAME_DISPLAY, node1_label),
    node2_label = if_else(node2_label == "first_name", FIRST_NAME_DISPLAY, node2_label)
  ) %>%
  group_by(node1_label, node2_label, node1_category, node2_category) %>%
  summarise(
    n_participants     = n_distinct(participant_id),
    mean_assoc         = mean(prop_together,   na.rm = TRUE),
    sd_assoc           = sd(prop_together,     na.rm = TRUE),
    mean_rt            = mean(median_rt_all,   na.rm = TRUE),
    sd_rt              = sd(median_rt_all,     na.rm = TRUE),
    # Within-person RT ratio: 1.0 = at person's average; <1 = faster; >1 = slower
    mean_ratio         = mean(rt_ratio,        na.rm = TRUE),
    n_obs              = n(),
    .groups = "drop"
  ) %>%
  filter(n_participants >= MIN_PARTICIPANTS_PER_PAIR) %>%
  mutate(
    centered_assoc       = mean_assoc - 0.5,
    assoc_strength       = abs(centered_assoc),
    assoc_direction      = if_else(centered_assoc >= 0, "together", "apart"),
    # Within-person RT factor: ratio=1 → 0.70, ratio=0.5 → 0.85, ratio≥2 → 0.40
    rt_factor            = pmax(0.4, 1 - 0.3 * coalesce(mean_ratio, 1.0)),
    rt_weighted_strength = assoc_strength * rt_factor,
    pair_id = map2_chr(node1_label, node2_label,
                       ~ paste(sort(c(.x, .y)), collapse = " | "))
  ) %>%
  filter(assoc_strength >= EDGE_THRESHOLD)

# ── 6. Build igraph / tidygraph object ────────────────────────────────────────

edge_df <- edge_summary %>%
  transmute(
    from                 = node1_label,
    to                   = node2_label,
    weight               = rt_weighted_strength,   # layout attraction
    centered_assoc       = centered_assoc,
    assoc_strength       = assoc_strength,
    assoc_direction      = assoc_direction,
    rt_factor            = rt_factor,
    mean_rt              = mean_rt,
    rt_weighted_strength = rt_weighted_strength,
    n_participants       = n_participants,
    n_obs                = n_obs
  )

active_nodes <- union(edge_df$from, edge_df$to)
node_df <- node_summary %>%
  filter(node %in% active_nodes) %>%
  rename(name = node)

# Normalize valence scores to full [0,1] range within this node set,
# then recompute fill colors so the full red–blue gradient is always visible.
node_df <- node_df %>%
  mutate(
    valence_norm = if_else(
      !is.na(valence_score),
      (valence_score - min(valence_score, na.rm = TRUE)) /
        (max(valence_score, na.rm = TRUE) - min(valence_score, na.rm = TRUE)),
      NA_real_
    ),
    node_fill = valence_to_color(valence_norm),
    # "me" node gets a star shape (shape 8 in ggplot2 = asterisk; use 23=diamond rotated as star proxy,
    # or shape 11 = filled star). We use 11 (open star) with fill.
    node_shape = if_else(name == FIRST_NAME_DISPLAY, 11L, if_else(category == "identity", 21L,
                   if_else(category == "policy", 22L, 23L)))
  )

# Precompute named fill/colour vectors for ggraph
node_fill_vec  <- setNames(node_df$node_fill,                 node_df$name)
node_color_vec <- setNames(CATEGORY_COLORS[node_df$category], node_df$name)
node_shape_vec <- setNames(node_df$node_shape,                node_df$name)

g <- tbl_graph(nodes = node_df, edges = edge_df, directed = FALSE)

# ── 6b. Node centrality (RT-weighted) ────────────────────────────────────────
#
# Weights = rt_weighted_strength so that fast, strong associations drive
# centrality more than slow or weak ones.

ig <- as.igraph(g)

node_centrality <- node_df %>%
  mutate(
    strength    = strength(ig,    weights = E(ig)$rt_weighted_strength),
    betweenness = betweenness(ig, weights = 1 / E(ig)$rt_weighted_strength,
                               normalized = TRUE),
    eigenvector = eigen_centrality(ig, weights = E(ig)$rt_weighted_strength)$vector
  ) %>%
  mutate(
    strength_norm    = (strength    - min(strength))    / diff(range(strength)),
    eigenvector_norm = (eigenvector - min(eigenvector)) / diff(range(eigenvector))
  ) %>%
  arrange(desc(strength))

write_csv(node_centrality, file.path(PLOT_DIR, "node_centrality.csv"))
message("✓ Saved node_centrality.csv")

# ── 7. PLOT 1 — Full Belief Network ──────────────────────────────────────────
#
# Node colour  = category
# Node size    = fixed
# Edge width   = RT-weighted strength (fast strong assoc → thicker)
# Edge colour  = direction (blue = together, orange = apart)
# Edge alpha   = assoc_strength (weak edges more transparent)
# Layout       = stress, attracted by RT-weighted strength

set.seed(42)

p1 <- ggraph(g, layout = "stress", weights = E(g)$rt_weighted_strength) +

  geom_edge_link(
    aes(
      width  = rt_weighted_strength,
      colour = assoc_direction,
      alpha  = assoc_strength
    ),
    lineend = "round"
  ) +
  scale_edge_width(range = c(0.4, 3.5), guide = "none") +
  scale_edge_alpha(range = c(0.3, 0.95), guide = "none") +
  scale_edge_colour_manual(
    values = c(together = "#3a82c4", apart = "#d45f3c"),
    labels = c(together = "Go together", apart = "Do not go together"),
    name   = "Perceived association"
  ) +

  geom_node_point(
    aes(fill = name, colour = name, shape = name),
    size = 11, stroke = 2.5
  ) +
  scale_fill_manual(values   = node_fill_vec,  guide = "none") +
  scale_colour_manual(values = node_color_vec, guide = "none") +
  scale_shape_manual(values  = node_shape_vec, guide = "none") +

  geom_node_text(
    aes(label = name),
    size = 3.1, fontface = "bold", colour = "white", repel = FALSE
  ) +

  theme_graph(base_family = "sans", background = "#0d1117") +
  theme(
    legend.position   = "right",
    legend.background = element_rect(fill = "#1a1f2e", colour = NA),
    legend.text       = element_text(colour = "#cccccc", size = 9),
    legend.title      = element_text(colour = "#aaaaaa", size = 9),
    plot.background   = element_rect(fill = "#0d1117", colour = NA),
    plot.title        = element_text(colour = "#eeeeee", size = 14,
                                     face = "bold", margin = margin(b = 6)),
    plot.subtitle     = element_text(colour = "#888888", size = 10,
                                     margin = margin(b = 14)),
    plot.caption      = element_text(colour = "#555555", size = 8)
  ) +
  labs(
    title    = "Belief Network — Political Concepts",
    subtitle = paste0(
      "Node fill = RT-weighted valence (red=bad · gray=neutral · blue=good)  ·  ",
      "Stroke = category  ·  ",
      "Edge width = RT-weighted strength  ·  Edge colour = association direction"
    ),
    caption  = paste0(
      "N = ", N_participants, " participants  ·  ",
      "Edge threshold: |prop_together − 0.5| ≥ ", EDGE_THRESHOLD, "  ·  ",
      "rt_factor = max(0.4, 1 − 0.3 × mean_ratio)  [within-person normalized]"
    )
  )

ggsave(file.path(PLOT_DIR, "01_full_network.pdf"),
       p1, width = 13, height = 9, device = cairo_pdf)
message("✓ Saved 01_full_network.pdf")

# ── 8. PLOT 2 — RT vs Association Strength ───────────────────────────────────
#
# Faster responses to strongly-associated pairs indicate more implicit/automatic
# retrieval. This scatter shows where speed and strength converge or diverge.

p2 <- edge_summary %>%
  mutate(
    pair_label = paste(node1_label, "×", node2_label),
    rt_ms_label = sprintf("%.0fms", mean_rt)
  ) %>%
  ggplot(aes(x = mean_rt, y = assoc_strength,
             colour = assoc_direction, label = pair_label)) +
  geom_point(aes(size = rt_weighted_strength), alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, colour = "white",
              linewidth = 0.7, alpha = 0.12, inherit.aes = FALSE,
              aes(x = mean_rt, y = assoc_strength)) +
  geom_text_repel(size = 2.5, max.overlaps = 12,
                  colour = "#cccccc", segment.colour = "#444444") +
  scale_colour_manual(
    values = c(together = "#3a82c4", apart = "#d45f3c"),
    labels = c(together = "Go together", apart = "Do not go together"),
    name   = "Association direction"
  ) +
  scale_size(range = c(2, 8), name = "RT-weighted\nstrength") +
  scale_x_continuous(labels = function(x) paste0(x, "ms")) +
  theme_minimal(base_family = "sans") +
  theme(
    panel.background  = element_rect(fill = "#1a1f2e", colour = NA),
    plot.background   = element_rect(fill = "#0d1117", colour = NA),
    axis.text         = element_text(colour = "#cccccc", size = 9),
    axis.title        = element_text(colour = "#aaaaaa", size = 10),
    legend.background = element_rect(fill = "#1a1f2e", colour = NA),
    legend.text       = element_text(colour = "#cccccc", size = 9),
    legend.title      = element_text(colour = "#aaaaaa", size = 9),
    plot.title        = element_text(colour = "#eeeeee", size = 13, face = "bold"),
    plot.subtitle     = element_text(colour = "#888888", size = 9),
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(colour = "#2a2f3e", linewidth = 0.4)
  ) +
  labs(
    x        = "Mean RT (ms)  [lower = faster / more automatic]",
    y        = "Association strength  |prop_together − 0.5|",
    title    = "Response Latency vs Association Strength",
    subtitle = paste0(
      "Point size = RT-weighted strength  ·  ",
      "Fast + strong = upper-left quadrant (most implicit associations)"
    ),
    caption  = paste0("N = ", N_participants, " participants")
  )

ggsave(file.path(PLOT_DIR, "02_rt_scatter.pdf"),
       p2, width = 11, height = 8, device = cairo_pdf)
message("✓ Saved 02_rt_scatter.pdf")

# ── 9. PLOT 3 — Association Heatmap ──────────────────────────────────────────

all_nodes_ordered <- node_summary %>%
  arrange(category, node) %>%
  pull(node)

heatmap_df <- edge_summary %>%
  transmute(
    node1 = factor(node1_label, levels = all_nodes_ordered),
    node2 = factor(node2_label, levels = all_nodes_ordered),
    val   = mean_assoc,
    rt    = mean_rt
  ) %>%
  bind_rows(transmute(., node1 = node2, node2 = node1, val = val, rt = rt)) %>%
  distinct()

p3 <- heatmap_df %>%
  ggplot(aes(x = node1, y = fct_rev(node2), fill = val)) +
  geom_tile(colour = "#0d1117", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", val)), size = 2.6, colour = "white") +
  scale_fill_gradient2(
    low = "#c0392b", mid = "#2a2f3e", high = "#3a82c4",
    midpoint = 0.5, limits = c(0, 1),
    name = "Prop.\ntogether"
  ) +
  scale_x_discrete(position = "top") +
  theme_minimal(base_family = "sans") +
  theme(
    panel.background  = element_rect(fill = "#0d1117", colour = NA),
    plot.background   = element_rect(fill = "#0d1117", colour = NA),
    axis.text.x       = element_text(colour = "#cccccc", size = 9,
                                     angle = 40, hjust = 0),
    axis.text.y       = element_text(colour = "#cccccc", size = 9),
    axis.title        = element_blank(),
    legend.background = element_rect(fill = "#1a1f2e", colour = NA),
    legend.text       = element_text(colour = "#cccccc", size = 8),
    legend.title      = element_text(colour = "#aaaaaa", size = 8),
    plot.title        = element_text(colour = "#eeeeee", size = 13, face = "bold"),
    plot.subtitle     = element_text(colour = "#888888", size = 9),
    panel.grid        = element_blank()
  ) +
  labs(
    title    = "Pairwise Association Heatmap",
    subtitle = paste0(
      "Cell values = mean proportion GO TOGETHER  ·  N = ", N_participants, " participants\n",
      "Blue = strongly associated  ·  Red = strongly disassociated"
    )
  )

ggsave(file.path(PLOT_DIR, "03_association_heatmap.pdf"),
       p3, width = 10, height = 9, device = cairo_pdf)
message("✓ Saved 03_association_heatmap.pdf")

# ── 10. PLOT 4 — Group Comparison Networks ────────────────────────────────────

group_map <- tribble(
  ~participant_type,          ~political_camp,
  "strong_progressive",       "Liberal",
  "moderate_progressive",     "Liberal",
  "moderate_democrat",        "Liberal",
  "union_democrat",           "Liberal",
  "liberal",                  "Liberal",
  "centrist_independent",     "Centrist",
  "libertarian",              "Centrist",
  "economic_populist",        "Centrist",
  "older_moderate",           "Centrist",
  "moderate_republican",      "Conservative",
  "suburban_republican",      "Conservative",
  "business_republican",      "Conservative",
  "never_trump_republican",   "Conservative",
  "conservative_republican",  "Conservative",
  "evangelical_republican",   "Conservative",
  "gun_rights_voter",         "Conservative",
  "maga_republican",          "Conservative",
  "conservative",             "Conservative",
  "maga",                     "Conservative"
)

make_group_network <- function(camp, b2_data, b1_data = NULL) {

  pids <- b2_data %>%
    left_join(group_map, by = "participant_type") %>%
    filter(political_camp == camp) %>%
    pull(participant_id) %>% unique()

  if (length(pids) < 2) return(NULL)

  grp_edges <- b2_data %>%
    filter(participant_id %in% pids) %>%
    mutate(
      node1_label = if_else(node1_label == "first_name", FIRST_NAME_DISPLAY, node1_label),
      node2_label = if_else(node2_label == "first_name", FIRST_NAME_DISPLAY, node2_label)
    ) %>%
    group_by(node1_label, node2_label, node1_category, node2_category) %>%
    summarise(
      mean_assoc = mean(prop_together, na.rm = TRUE),
      mean_ratio = mean(rt_ratio,      na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      centered_assoc       = mean_assoc - 0.5,
      assoc_strength       = abs(centered_assoc),
      assoc_direction      = if_else(centered_assoc >= 0, "together", "apart"),
      rt_factor            = pmax(0.4, 1 - 0.3 * coalesce(mean_ratio, 1.0)),
      rt_weighted_strength = assoc_strength * rt_factor
    ) %>%
    filter(assoc_strength >= EDGE_THRESHOLD) %>%
    transmute(
      from = node1_label, to = node2_label,
      weight = rt_weighted_strength,
      centered_assoc, assoc_strength, assoc_direction,
      rt_factor, rt_weighted_strength
    )

  active <- union(grp_edges$from, grp_edges$to)
  node_sub <- node_summary %>%
    filter(node %in% active) %>%
    rename(name = node)

  if (nrow(node_sub) == 0 || nrow(grp_edges) == 0) return(NULL)

  g_grp <- tbl_graph(nodes = node_sub, edges = grp_edges, directed = FALSE)

  # Compute valence from this camp's own B1 rows, then normalize within the set
  if (!is.null(b1_data) && nrow(b1_data) > 0) {
    camp_valence <- b1_data %>%
      filter(participant_id %in% pids) %>%
      mutate(node = if_else(target_label == "first_name", FIRST_NAME_DISPLAY, target_label)) %>%
      group_by(node) %>%
      summarise(
        mean_prop_good = mean(prop_good,  na.rm = TRUE),
        mean_ratio_b1  = mean(rt_ratio,   na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        rt_factor_b1     = pmax(0.4, 1 - 0.3 * coalesce(mean_ratio_b1, 1.0)),
        valence_score     = 0.5 + (mean_prop_good - 0.5) * rt_factor_b1
      )
    node_sub <- node_sub %>%
      select(-any_of(c("mean_prop_good", "valence_score"))) %>%
      left_join(camp_valence %>% select(node = node, mean_prop_good, valence_score),
                by = c("name" = "node"))
  }

  node_sub <- node_sub %>%
    mutate(
      valence_norm = if_else(
        !is.na(valence_score),
        (valence_score - min(valence_score, na.rm = TRUE)) /
          (max(valence_score, na.rm = TRUE) - min(valence_score, na.rm = TRUE)),
        NA_real_
      ),
      node_fill  = valence_to_color(valence_norm),
      node_shape = if_else(name == FIRST_NAME_DISPLAY, 11L, if_else(category == "identity", 21L,
                     if_else(category == "policy", 22L, 23L)))
    )

  sub_fill_vec  <- setNames(node_sub$node_fill,                 node_sub$name)
  sub_color_vec <- setNames(CATEGORY_COLORS[node_sub$category], node_sub$name)
  sub_shape_vec <- setNames(node_sub$node_shape,                node_sub$name)

  set.seed(42)
  ggraph(g_grp, layout = "stress", weights = E(g_grp)$rt_weighted_strength) +
    geom_edge_link(
      aes(width = rt_weighted_strength, colour = assoc_direction, alpha = assoc_strength),
      lineend = "round"
    ) +
    scale_edge_width(range = c(0.4, 3.0), guide = "none") +
    scale_edge_alpha(range = c(0.3, 0.9),  guide = "none") +
    scale_edge_colour_manual(
      values = c(together = "#3a82c4", apart = "#d45f3c"),
      guide  = "none"
    ) +
    geom_node_point(
      aes(fill = name, colour = name, shape = name),
      size = 10, stroke = 2.0
    ) +
    scale_fill_manual(values   = sub_fill_vec,  guide = "none") +
    scale_colour_manual(values = sub_color_vec, guide = "none") +
    scale_shape_manual(values  = sub_shape_vec, guide = "none") +
    geom_node_text(
      aes(label = name), size = 2.8, fontface = "bold", colour = "white"
    ) +
    theme_graph(base_family = "sans", background = "#0d1117") +
    theme(
      plot.title    = element_text(colour = "#eeeeee", size = 12, face = "bold"),
      plot.subtitle = element_text(colour = "#888888", size = 8.5)
    ) +
    labs(
      title    = camp,
      subtitle = paste0("n = ", length(pids), " participants")
    )
}

camps       <- c("Liberal", "Centrist", "Conservative")
camp_plots  <- map(camps, make_group_network, b2_data = b2_raw, b1_data = b1_raw) %>% compact()

if (length(camp_plots) >= 2) {
  p4 <- wrap_plots(camp_plots, nrow = 1) +
    plot_annotation(
      title    = "Belief Networks by Political Camp",
      subtitle = paste0(
        "Node colour = category  ·  Edge colour: blue = GO TOGETHER, red = DO NOT  ·  ",
        "Edge width = RT-weighted strength"
      ),
      caption  = paste0("Edge threshold: |prop_together − 0.5| ≥ ", EDGE_THRESHOLD),
      theme    = theme(
        plot.background = element_rect(fill = "#0d1117", colour = NA),
        plot.title      = element_text(colour = "#eeeeee", size = 14, face = "bold"),
        plot.subtitle   = element_text(colour = "#888888", size = 9),
        plot.caption    = element_text(colour = "#555555", size = 8)
      )
    )

  ggsave(file.path(PLOT_DIR, "04_group_networks.pdf"),
         p4, width = 18, height = 7, device = cairo_pdf)
  message("✓ Saved 04_group_networks.pdf")
} else {
  message("⚠ Not enough groups to produce group comparison plot.")
}

# ── 11. PLOT 5 — Me Connections ───────────────────────────────────────────────

me_edges <- edge_summary %>%
  filter(node1_label == FIRST_NAME_DISPLAY | node2_label == FIRST_NAME_DISPLAY) %>%
  mutate(
    other_node     = if_else(node1_label == FIRST_NAME_DISPLAY, node2_label, node1_label),
    other_category = if_else(node1_label == FIRST_NAME_DISPLAY, node2_category, node1_category)
  ) %>%
  mutate(other_node = fct_reorder(other_node, centered_assoc))

p5 <- me_edges %>%
  ggplot(aes(x = centered_assoc, y = other_node, fill = assoc_direction)) +
  geom_col(aes(x = if_else(assoc_direction == "together",
                             rt_weighted_strength, -rt_weighted_strength)),
           width = 0.65, alpha = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "#888888", linewidth = 0.6) +
  geom_text(
    aes(
      x     = centered_assoc + if_else(centered_assoc >= 0, 0.02, -0.02),
      label = sprintf("%.2f (%.0fms)", mean_assoc, mean_rt),
      hjust = if_else(centered_assoc >= 0, 0, 1)
    ),
    size = 2.8, colour = "white"
  ) +
  facet_grid(other_category ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(
    values = c(together = "#3a82c4", apart = "#d45f3c"),
    labels = c(together = "Go together with me", apart = "Do not go together with me"),
    name   = NULL
  ) +
  scale_x_continuous(
    limits = c(-0.55, 0.7),
    breaks = c(-0.5, -0.25, 0, 0.25, 0.5),
    labels = c("−0.5", "−0.25", "0", "+0.25", "+0.5")
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    panel.background   = element_rect(fill = "#1a1f2e", colour = NA),
    plot.background    = element_rect(fill = "#0d1117", colour = NA),
    strip.background   = element_rect(fill = "#252b3b", colour = NA),
    strip.text         = element_text(colour = "#aaaaaa", size = 9),
    axis.text          = element_text(colour = "#cccccc", size = 9.5),
    axis.title         = element_text(colour = "#aaaaaa", size = 10),
    legend.position    = "bottom",
    legend.background  = element_rect(fill = "#1a1f2e", colour = NA),
    legend.text        = element_text(colour = "#cccccc", size = 9),
    plot.title         = element_text(colour = "#eeeeee", size = 13, face = "bold"),
    plot.subtitle      = element_text(colour = "#888888", size = 9),
    plot.caption       = element_text(colour = "#555555", size = 8),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(colour = "#2a2f3e", linewidth = 0.4)
  ) +
  labs(
    x       = "RT-weighted association with 'me'  (bar width = rt_weighted_strength, direction = centered_assoc)",
    y       = NULL,
    title   = "Self Connections — Association strength & speed with 'me'",
    subtitle = paste0(
      "Bar width = RT-weighted strength  ·  Label = prop_together (mean RT)  ·  ",
      "N = ", N_participants, " participants"
    ),
    caption = paste0(
      "Only pairs with ≥ ", MIN_PARTICIPANTS_PER_PAIR,
      " participants and |assoc − 0.5| ≥ ", EDGE_THRESHOLD
    )
  )

ggsave(file.path(PLOT_DIR, "05_me_connections.pdf"),
       p5, width = 9, height = 8, device = cairo_pdf)
message("✓ Saved 05_me_connections.pdf")

# ── 12. Console summary ───────────────────────────────────────────────────────

cat("\n══════════════════════════════════════════════════════\n")
cat("  Belief Network Summary  (RT-weighted)\n")
cat("══════════════════════════════════════════════════════\n\n")

cat(sprintf("Participants: %d\n", N_participants))
cat(sprintf("Nodes in network: %d\n", nrow(node_df)))
cat(sprintf("Edges above threshold: %d / %d possible pairs\n",
            nrow(edge_summary), choose(nrow(node_summary), 2)))
cat(sprintf("Mean RT across pairs: %.0fms  (mean_ratio: %.2f – %.2f  |  rt_factor: %.2f – %.2f)\n\n",
            mean(edge_summary$mean_rt,    na.rm = TRUE),
            min(edge_summary$mean_ratio,  na.rm = TRUE),
            max(edge_summary$mean_ratio,  na.rm = TRUE),
            min(edge_summary$rt_factor),  max(edge_summary$rt_factor)))

cat("Top 5 FASTEST + STRONGEST pairs (highest rt_weighted_strength):\n")
edge_summary %>%
  arrange(desc(rt_weighted_strength)) %>%
  slice_head(n = 5) %>%
  mutate(row = sprintf("  %-16s × %-24s  assoc=%.2f  rt=%.0fms  rtwt=%.3f",
                       node1_label, node2_label, mean_assoc, mean_rt, rt_weighted_strength)) %>%
  pull(row) %>% cat(sep = "\n")

cat("\n\nTop 5 most ASSOCIATED pairs (prop_together):\n")
edge_summary %>%
  arrange(desc(mean_assoc)) %>%
  slice_head(n = 5) %>%
  mutate(row = sprintf("  %-16s × %-24s  %.2f  (%.0fms)",
                       node1_label, node2_label, mean_assoc, mean_rt)) %>%
  pull(row) %>% cat(sep = "\n")

cat("\n\nNode centrality — RT-weighted strength (top 10):\n")
cat(sprintf("  %-26s  %8s  %11s  %11s\n",
            "Node", "Strength", "Betweenness", "Eigenvector"))
cat("  ", strrep("─", 62), "\n", sep = "")
node_centrality %>%
  slice_head(n = 10) %>%
  mutate(row = sprintf("  %-26s  %8.3f  %11.3f  %11.3f",
                       name, strength_norm, betweenness, eigenvector_norm)) %>%
  pull(row) %>% cat(sep = "\n")

cat(sprintf("\n  (Full table saved to plots/node_centrality.csv)\n"))
cat("\n✓ All plots saved to:", normalizePath(PLOT_DIR), "\n")
cat("══════════════════════════════════════════════════════\n\n")
