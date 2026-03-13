# ============================================================
# belief_network.R
#
# Constructs and visualizes a weighted belief network from the
# two-block hybrid task output.
#
# Block 1 (node valence)  → prop_good per concept
# Block 2 (edge strength) → prop_together per pair
#
# The network represents how a sample of people mentally represent
# a set of political concepts: where each node sits on a
# good/bad dimension, and how strongly each pair is associated.
#
# INPUTS (from simulate.html or real experiment):
#   simulated_data/block1_target_summary.csv
#   simulated_data/block2_pair_summary.csv
#
# OUTPUTS (written to plots/):
#   01_full_network.pdf          main belief network, full sample
#   02_valence_bar.pdf           node valence sorted bar chart
#   03_association_heatmap.pdf   pair-level association matrix
#   04_group_networks.pdf        liberal vs. conservative networks
#   05_valence_x_association.pdf scatter: edge assoc vs. node valence difference
# ============================================================

# ── 0. Packages ──────────────────────────────────────────────────────────────

required <- c("tidyverse", "igraph", "ggraph", "tidygraph",
              "graphlayouts", "scales", "patchwork", "ggrepel", "viridis")

for (pkg in required) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  library(pkg, character.only = TRUE)
}

# ── 1. Paths & directories ───────────────────────────────────────────────────

# Auto-detect data location: look for summary CSVs in simulated_data/ first,
# then fall back to the script's own directory; if only all_trials.csv is
# present, compute summaries on the fly.

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
# Strategy: use pre-computed summary CSVs if available.
# If not, derive them from all_trials.csv (output of simulate.html or the
# real experiment).

b1_path <- tryCatch(find_file("block1_target_summary.csv"), error = function(e) NULL)
b2_path <- tryCatch(find_file("block2_pair_summary.csv"),   error = function(e) NULL)

if (!is.null(b1_path) && !is.null(b2_path)) {
  # Pre-computed summaries exist — load directly
  b1_raw <- read_csv(b1_path, show_col_types = FALSE)
  b2_raw <- read_csv(b2_path, show_col_types = FALSE)
  message("✓ Loaded pre-computed summary files")

} else {
  # Fall back: derive summaries from all_trials.csv
  message("⚙ Summary files not found — computing from all_trials.csv")
  trials_path <- find_file("all_trials.csv")
  trials      <- read_csv(trials_path, show_col_types = FALSE)

  b1_raw <- trials %>%
    filter(task_type == "single_target") %>%
    group_by(participant_id, participant_type, target_label, target_category,
             mapping_condition) %>%
    summarise(
      n_trials   = n(),
      n_good     = sum(response_label == "good", na.rm = TRUE),
      n_bad      = sum(response_label == "bad",  na.rm = TRUE),
      prop_good  = n_good / n_trials,
      prop_bad   = n_bad  / n_trials,
      median_rt_all  = median(rt_ms[timeout == 0], na.rm = TRUE),
      median_rt_good = median(rt_ms[response_label == "good" & timeout == 0], na.rm = TRUE),
      median_rt_bad  = median(rt_ms[response_label == "bad"  & timeout == 0], na.rm = TRUE),
      timeout_rate   = mean(timeout, na.rm = TRUE),
      .groups = "drop"
    )

  b2_raw <- trials %>%
    filter(task_type == "pairwise") %>%
    # Canonicalise pair order (alphabetically) so left/right randomisation
    # doesn't create duplicate pair keys. Column names are from the unified
    # all_trials.csv schema produced by simulate.html / hybrid_task.html.
    mutate(
      n1     = pmin(left_node_label,  right_node_label),
      n2     = pmax(left_node_label,  right_node_label),
      cat1   = if_else(left_node_label  == n1, left_node_category,  right_node_category),
      cat2   = if_else(left_node_label  == n1, right_node_category, left_node_category),
      pid    = paste(n1, n2, sep = " | ")
    ) %>%
    group_by(participant_id, participant_type,
             pair_id = pid,
             node1_label = n1, node2_label = n2,
             node1_category = cat1, node2_category = cat2,
             mapping_condition) %>%
    summarise(
      n_trials       = n(),
      n_together     = sum(response_label == "together", na.rm = TRUE),
      n_apart        = sum(response_label == "apart",    na.rm = TRUE),
      prop_together  = n_together / n_trials,
      prop_apart     = n_apart    / n_trials,
      median_rt_all      = median(rt_ms[timeout == 0], na.rm = TRUE),
      median_rt_together = median(rt_ms[response_label == "together" & timeout == 0],
                                  na.rm = TRUE),
      median_rt_apart    = median(rt_ms[response_label == "apart"    & timeout == 0],
                                  na.rm = TRUE),
      timeout_rate       = mean(timeout, na.rm = TRUE),
      .groups = "drop"
    )

  message("✓ Summaries computed from all_trials.csv")
}

# Participant count — derived from data, used in plot labels throughout
N_participants <- n_distinct(b1_raw$participant_id)

# ── 3. Configuration ─────────────────────────────────────────────────────────

# Edge threshold: only show edges where |mean_assoc - 0.5| exceeds this.
# Raise to reduce clutter; lower to show all edges.
EDGE_THRESHOLD <- 0.05

# Node label to display for the personalised self node
FIRST_NAME_DISPLAY <- "Name"

# Category colour palette (used for node borders and group labels)
CATEGORY_COLORS <- c(
  identity = "#5b8dd9",
  policy   = "#6ac46a",
  values   = "#c97de8"
)

# ── 4. Node-level aggregation ─────────────────────────────────────────────────
# One row per concept; average across participants.

node_summary <- b1_raw %>%
  mutate(
    # Rename the personalised self node to a fixed display label
    target_label = if_else(target_label == "first_name",
                           FIRST_NAME_DISPLAY, target_label)
  ) %>%
  group_by(target_label, target_category) %>%
  summarise(
    mean_valence = mean(prop_good,    na.rm = TRUE),
    sd_valence   = sd(prop_good,      na.rm = TRUE),
    n_obs        = n(),
    .groups = "drop"
  ) %>%
  rename(node = target_label, category = target_category)

# ── 5. Edge-level aggregation ─────────────────────────────────────────────────
# One row per pair; pairwise complete observations (average prop_together only
# for participants who actually saw that pair — pairs unseen by a participant
# contribute NA and are excluded via na.rm = TRUE).
#
# With 30 possible nodes, 10 per participant (me always included), each
# non-me pair is expected to appear together ~(9/29)^2 × N times.
# MIN_PARTICIPANTS_PER_PAIR sets a quality floor; pairs seen by fewer
# participants are excluded as unreliable estimates.
MIN_PARTICIPANTS_PER_PAIR <- 3

edge_summary <- b2_raw %>%
  mutate(
    node1_label = if_else(node1_label == "first_name", FIRST_NAME_DISPLAY, node1_label),
    node2_label = if_else(node2_label == "first_name", FIRST_NAME_DISPLAY, node2_label)
  ) %>%
  group_by(node1_label, node2_label, node1_category, node2_category) %>%
  summarise(
    # n_participants = number of distinct participants who saw this pair
    n_participants = n_distinct(participant_id),
    mean_assoc     = mean(prop_together, na.rm = TRUE),
    sd_assoc       = sd(prop_together,   na.rm = TRUE),
    n_obs          = n(),
    .groups = "drop"
  ) %>%
  # Drop pairs seen by too few participants (unreliable estimates)
  filter(n_participants >= MIN_PARTICIPANTS_PER_PAIR) %>%
  mutate(
    # Centre at 0.5: positive = perceived as going together,
    #                negative = perceived as NOT going together
    centered_assoc    = mean_assoc - 0.5,
    assoc_strength    = abs(centered_assoc),
    assoc_direction   = if_else(centered_assoc >= 0, "together", "apart"),
    # Unique pair ID (sorted so order doesn't matter)
    pair_id = map2_chr(node1_label, node2_label,
                       ~ paste(sort(c(.x, .y)), collapse = " | "))
  ) %>%
  # Keep only edges above association threshold
  filter(assoc_strength >= EDGE_THRESHOLD)

# ── 6. Build igraph / tidygraph object ────────────────────────────────────────

edge_df <- edge_summary %>%
  transmute(
    from            = node1_label,
    to              = node2_label,
    weight          = mean_assoc,          # raw 0-1 (used for layout)
    centered_assoc  = centered_assoc,
    assoc_strength  = assoc_strength,
    assoc_direction = assoc_direction,
    n_participants  = n_participants,      # number of participants who saw this pair
    n_obs           = n_obs
  )

# Only include nodes that actually appear in the edge list
active_nodes <- union(edge_df$from, edge_df$to)
node_df <- node_summary %>%
  filter(node %in% active_nodes) %>%
  rename(name = node)

g <- tbl_graph(
  nodes    = node_df,
  edges    = edge_df,
  directed = FALSE
)

# ── 6b. Node centrality ───────────────────────────────────────────────────────
# All three measures use assoc_strength (absolute association) as edge weights,
# so both "go together" and "do not go together" strong associations count.
#
#   strength    – weighted degree: sum of incident edge strengths.
#                 Captures total cognitive connectedness.
#   betweenness – how often a node lies on shortest paths between other nodes.
#                 Identifies conceptual bridges between clusters.
#   eigenvector – importance weighted by the importance of neighbours.
#                 Captures whether a node is connected to other hub nodes.

ig <- as.igraph(g)   # plain igraph object for centrality functions

node_centrality <- node_df %>%
  mutate(
    strength    = strength(ig,    weights = E(ig)$assoc_strength),
    betweenness = betweenness(ig, weights = 1 / E(ig)$assoc_strength,
                               normalized = TRUE),
    eigenvector = eigen_centrality(ig, weights = E(ig)$assoc_strength)$vector
  ) %>%
  # Normalise strength and eigenvector to [0, 1] for comparability
  mutate(
    strength_norm    = (strength    - min(strength))    / diff(range(strength)),
    eigenvector_norm = (eigenvector - min(eigenvector)) / diff(range(eigenvector))
  ) %>%
  arrange(desc(strength))

write_csv(node_centrality,
          file.path(PLOT_DIR, "node_centrality.csv"))
message("✓ Saved node_centrality.csv")

# ── 7. PLOT 1 — Full Belief Network ──────────────────────────────────────────
#
# Node colour  = valence (red → green gradient)
# Node size    = fixed
# Edge width   = association strength
# Edge colour  = direction (blue = together, orange = apart)
# Layout       = weighted Fruchterman-Reingold (high weight → closer)

# Stress layout (graphlayouts): minimises a graph-distance stress function,
# which reliably places the most central (hub) nodes at the geometric centre.
# Weights = positive association strength; "go together" pairs are pulled close.
# The layout is deterministic so set.seed is not strictly needed but kept for
# reproducibility of any small tie-breaking steps.
set.seed(42)

p1 <- ggraph(g, layout = "stress",
             weights = pmax(E(g)$centered_assoc, 0.01)) +

  # Edges
  geom_edge_link(
    aes(
      width  = assoc_strength,
      colour = assoc_direction,
      alpha  = assoc_strength
    ),
    lineend = "round"
  ) +
  scale_edge_width(range = c(0.4, 3.2), guide = "none") +
  scale_edge_alpha(range = c(0.35, 0.95), guide = "none") +
  scale_edge_colour_manual(
    values = c(together = "#3a82c4", apart = "#d45f3c"),
    labels = c(together = "Go together", apart = "Do not go together"),
    name   = "Perceived association"
  ) +

  # Nodes — filled by valence, bordered by category
  geom_node_point(
    aes(fill = mean_valence, colour = category),
    shape = 21, size = 11, stroke = 2.2
  ) +
  scale_fill_gradient2(
    low      = "#c0392b",
    mid      = "#f5e642",
    high     = "#27ae60",
    midpoint = 0.5,
    limits   = c(0, 1),
    name     = "Mean valence\n(prop. rated Good)"
  ) +
  scale_colour_manual(
    values = CATEGORY_COLORS,
    name   = "Category"
  ) +

  # Node labels
  geom_node_text(
    aes(label = name),
    size = 3.1, fontface = "bold", colour = "white",
    repel = FALSE
  ) +

  # Theme
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
      "Node colour = valence (Block 1: proportion rated GOOD)  ·  ",
      "Edge width & colour = association (Block 2: proportion GO TOGETHER)"
    ),
    caption  = paste0("N = ", N_participants, " participants  ·  Edge threshold: |prop_together − 0.5| ≥ ", EDGE_THRESHOLD)
  )

ggsave(file.path(PLOT_DIR, "01_full_network.pdf"),
       p1, width = 13, height = 9, device = cairo_pdf)
message("✓ Saved 01_full_network.pdf")

# ── 8. PLOT 2 — Node Valence Bar Chart ───────────────────────────────────────

p2 <- node_summary %>%
  mutate(
    node         = fct_reorder(node, mean_valence),
    valence_label = sprintf("%.2f", mean_valence),
    # Confidence interval (approximate SE from binomial)
    se           = sqrt(mean_valence * (1 - mean_valence) / n_obs),
    ci_lo        = pmax(0, mean_valence - 1.96 * se),
    ci_hi        = pmin(1, mean_valence + 1.96 * se)
  ) %>%
  ggplot(aes(x = mean_valence, y = node, fill = mean_valence)) +
  geom_col(width = 0.65, alpha = 0.9) +
  geom_errorbarh(
    aes(xmin = ci_lo, xmax = ci_hi),
    height = 0.25, colour = "white", linewidth = 0.5
  ) +
  geom_vline(xintercept = 0.5, linetype = "dashed",
             colour = "#888888", linewidth = 0.6) +
  geom_text(
    aes(label = valence_label),
    hjust = -0.2, size = 3.2, colour = "white"
  ) +
  facet_grid(category ~ ., scales = "free_y", space = "free_y",
             labeller = label_value) +
  scale_fill_gradient2(
    low = "#c0392b", mid = "#f5e642", high = "#27ae60",
    midpoint = 0.5, limits = c(0, 1), guide = "none"
  ) +
  scale_x_continuous(
    limits = c(0, 1.12), breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = c("0", ".25", ".50", ".75", "1.0")
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    panel.background  = element_rect(fill = "#1a1f2e", colour = NA),
    plot.background   = element_rect(fill = "#0d1117", colour = NA),
    strip.background  = element_rect(fill = "#252b3b", colour = NA),
    strip.text        = element_text(colour = "#aaaaaa", size = 9),
    axis.text         = element_text(colour = "#cccccc", size = 9.5),
    axis.title        = element_text(colour = "#aaaaaa", size = 10),
    plot.title        = element_text(colour = "#eeeeee", size = 13, face = "bold"),
    plot.subtitle     = element_text(colour = "#888888", size = 9),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(colour = "#2a2f3e", linewidth = 0.4)
  ) +
  labs(
    x        = "Mean proportion rated GOOD  (Block 1)",
    y        = NULL,
    title    = "Node Valence — Mean Evaluation per Concept",
    subtitle = "Error bars = 95% CI  ·  Dashed line = neutral (0.50)"
  )

ggsave(file.path(PLOT_DIR, "02_valence_bar.pdf"),
       p2, width = 9, height = 7, device = cairo_pdf)
message("✓ Saved 02_valence_bar.pdf")

# ── 9. PLOT 3 — Association Heatmap ──────────────────────────────────────────

# Build a full symmetric matrix of all pairwise associations
all_nodes_ordered <- node_summary %>%
  arrange(category, desc(mean_valence)) %>%
  pull(node)

# Add diagonal (self-loops, not in data, set to NA)
heatmap_df <- edge_summary %>%
  transmute(
    node1   = factor(node1_label, levels = all_nodes_ordered),
    node2   = factor(node2_label, levels = all_nodes_ordered),
    val     = mean_assoc
  ) %>%
  # Mirror: add the reversed direction too
  bind_rows(transmute(., node1 = node2, node2 = node1, val = val)) %>%
  distinct()

p3 <- heatmap_df %>%
  ggplot(aes(x = node1, y = fct_rev(node2), fill = val)) +
  geom_tile(colour = "#0d1117", linewidth = 0.6) +
  geom_text(
    aes(label = sprintf("%.2f", val)),
    size = 2.8, colour = "white"
  ) +
  scale_fill_gradient2(
    low      = "#c0392b",
    mid      = "#2a2f3e",
    high     = "#3a82c4",
    midpoint = 0.5,
    limits   = c(0, 1),
    name     = "Prop.\ntogether"
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
    title    = "Pairwise Association Heatmap — Block 2",
    subtitle = paste0("Cell values = mean proportion of GO TOGETHER responses  ·  N = ", N_participants, " participants\nBlue = strongly associated  ·  Red = strongly disassociated")
  )

ggsave(file.path(PLOT_DIR, "03_association_heatmap.pdf"),
       p3, width = 10, height = 9, device = cairo_pdf)
message("✓ Saved 03_association_heatmap.pdf")

# ── 10. PLOT 4 — Group Comparison Networks ───────────────────────────────────
# Side-by-side networks for liberal vs. conservative participants.

# Classify each participant as liberal or conservative based on participant_type
group_map <- tribble(
  ~participant_type,        ~political_camp,
  "strong_progressive",     "Liberal",
  "moderate_democrat",      "Liberal",
  "union_democrat",         "Liberal",
  "centrist_independent",   "Centrist",
  "libertarian",            "Centrist",
  "moderate_republican",    "Conservative",
  "suburban_republican",    "Conservative",
  "business_republican",    "Conservative",
  "never_trump_republican", "Conservative",
  "conservative_republican","Conservative",
  "evangelical_republican", "Conservative",
  "maga_republican",        "Conservative",
  "gun_rights_voter",       "Conservative",
  "economic_populist",      "Centrist",
  "older_moderate",         "Centrist"
)

make_group_network <- function(camp, b1_data, b2_data, node_summary_all) {

  pids <- b1_data %>%
    left_join(group_map, by = "participant_type") %>%
    filter(political_camp == camp) %>%
    pull(participant_id) %>% unique()

  # Node valences for this group
  grp_nodes <- b1_data %>%
    filter(participant_id %in% pids) %>%
    mutate(target_label = if_else(target_label == "first_name",
                                  FIRST_NAME_DISPLAY, target_label)) %>%
    group_by(target_label, target_category) %>%
    summarise(mean_valence = mean(prop_good, na.rm = TRUE), .groups = "drop") %>%
    rename(name = target_label, category = target_category)

  # Edge associations for this group
  grp_edges <- b2_data %>%
    filter(participant_id %in% pids) %>%
    mutate(
      node1_label = if_else(node1_label == "first_name", FIRST_NAME_DISPLAY, node1_label),
      node2_label = if_else(node2_label == "first_name", FIRST_NAME_DISPLAY, node2_label)
    ) %>%
    group_by(node1_label, node2_label) %>%
    summarise(mean_assoc = mean(prop_together, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      centered_assoc  = mean_assoc - 0.5,
      assoc_strength  = abs(centered_assoc),
      assoc_direction = if_else(centered_assoc >= 0, "together", "apart")
    ) %>%
    filter(assoc_strength >= EDGE_THRESHOLD) %>%
    transmute(
      from = node1_label, to = node2_label,
      weight = mean_assoc, centered_assoc, assoc_strength, assoc_direction
    )

  active <- union(grp_edges$from, grp_edges$to)
  node_sub <- grp_nodes %>% filter(name %in% active)
  if (nrow(node_sub) == 0 || nrow(grp_edges) == 0) return(NULL)

  g_grp <- tbl_graph(nodes = node_sub, edges = grp_edges, directed = FALSE)

  set.seed(42)
  ggraph(g_grp, layout = "stress",
         weights = pmax(E(g_grp)$centered_assoc, 0.01)) +
    geom_edge_link(
      aes(width = assoc_strength, colour = assoc_direction, alpha = assoc_strength),
      lineend = "round"
    ) +
    scale_edge_width(range = c(0.4, 3.0), guide = "none") +
    scale_edge_alpha(range = c(0.3, 0.9),  guide = "none") +
    scale_edge_colour_manual(
      values = c(together = "#3a82c4", apart = "#d45f3c"),
      guide  = "none"
    ) +
    geom_node_point(
      aes(fill = mean_valence, colour = category),
      shape = 21, size = 10, stroke = 2.0
    ) +
    scale_fill_gradient2(
      low = "#c0392b", mid = "#f5e642", high = "#27ae60",
      midpoint = 0.5, limits = c(0, 1), guide = "none"
    ) +
    scale_colour_manual(values = CATEGORY_COLORS, guide = "none") +
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

camps <- c("Liberal", "Centrist", "Conservative")
camp_plots <- map(camps, make_group_network,
                  b1_data = b1_raw, b2_data = b2_raw,
                  node_summary_all = node_summary) %>%
  compact()                                   # drop NULLs

if (length(camp_plots) >= 2) {
  p4 <- wrap_plots(camp_plots, nrow = 1) +
    plot_annotation(
      title    = "Belief Networks by Political Camp",
      subtitle = "Node colour = valence  ·  Edge colour: blue = GO TOGETHER, red = DO NOT GO TOGETHER  ·  Edge width = strength",
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

# ── 11. PLOT 5 — Valence Gap vs. Association Strength ────────────────────────
# Do pairs with larger valence differences get rated as NOT going together?
# This tests whether the belief network has a consistency/balance structure.

valence_lookup <- node_summary %>%
  select(node, mean_valence)

valence_edge_df <- edge_summary %>%
  left_join(valence_lookup, by = c("node1_label" = "node")) %>%
  rename(val1 = mean_valence) %>%
  left_join(valence_lookup, by = c("node2_label" = "node")) %>%
  rename(val2 = mean_valence) %>%
  mutate(
    valence_gap      = abs(val1 - val2),
    valence_product  = val1 * val2,              # both good → high
    pair_label       = paste(node1_label, "×", node2_label),
    # Colour by whether both nodes are positively or negatively valenced
    both_valence_dir = case_when(
      val1 > 0.5 & val2 > 0.5 ~ "Both rated GOOD",
      val1 < 0.5 & val2 < 0.5 ~ "Both rated BAD",
      TRUE                    ~ "Mixed valence"
    )
  )

# Correlation label
cor_r <- with(valence_edge_df, cor(valence_gap, centered_assoc, use = "pairwise.complete.obs"))
cor_label <- sprintf("r = %.2f", cor_r)

p5 <- valence_edge_df %>%
  ggplot(aes(x = valence_gap, y = centered_assoc,
             colour = both_valence_dir, label = pair_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "#555555") +
  geom_point(aes(size = assoc_strength), alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, colour = "white",
              linewidth = 0.8, alpha = 0.15, inherit.aes = FALSE,
              aes(x = valence_gap, y = centered_assoc)) +
  geom_text_repel(
    size = 2.7, max.overlaps = 15,
    colour = "#cccccc", segment.colour = "#444444"
  ) +
  annotate("text", x = 0.05, y = max(valence_edge_df$centered_assoc) - 0.03,
           label = cor_label, colour = "#aaaaaa", size = 3.8) +
  scale_colour_manual(
    values = c("Both rated GOOD" = "#27ae60",
               "Both rated BAD"  = "#c0392b",
               "Mixed valence"   = "#e8883a"),
    name = "Pair valence profile"
  ) +
  scale_size(range = c(2, 7), guide = "none") +
  scale_x_continuous(limits = c(0, NA)) +
  scale_y_continuous(
    breaks = seq(-0.5, 0.5, 0.1),
    labels = function(x) sprintf("%+.1f", x)
  ) +
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
    x        = "Valence gap between the two nodes  |prop_good₁ − prop_good₂|",
    y        = "Centered association  (prop_together − 0.5)",
    title    = "Structural Balance: Valence Gap vs. Perceived Association",
    subtitle = paste0(
      "Positive y = pairs rated GO TOGETHER  ·  Negative y = DO NOT GO TOGETHER\n",
      "Cognitive balance predicts: similar-valence pairs → associated, different-valence → not associated"
    )
  )

ggsave(file.path(PLOT_DIR, "05_valence_x_association.pdf"),
       p5, width = 11, height = 8, device = cairo_pdf)
message("✓ Saved 05_valence_x_association.pdf")

# ── 12. PLOT 6 — Me Connections ───────────────────────────────────────────────
# Which concepts are most (and least) associated with the self node?
# Bar chart of centered association for every me-pair that cleared the coverage
# threshold, split by concept category and sorted by association strength.

self_label <- "me"   # matches the node label set in CONFIG

me_edges <- edge_summary %>%
  filter(node1_label == self_label | node2_label == self_label) %>%
  mutate(
    other_node     = if_else(node1_label == self_label, node2_label, node1_label),
    other_category = if_else(node1_label == self_label, node2_category, node1_category)
  ) %>%
  left_join(node_summary %>% select(node, mean_valence),
            by = c("other_node" = "node")) %>%
  mutate(other_node = fct_reorder(other_node, centered_assoc))

p6 <- me_edges %>%
  ggplot(aes(x = centered_assoc, y = other_node, fill = assoc_direction)) +
  geom_col(width = 0.65, alpha = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "#888888", linewidth = 0.6) +
  geom_text(
    aes(
      x     = centered_assoc + if_else(centered_assoc >= 0,  0.02, -0.02),
      label = sprintf("%.2f", mean_assoc),
      hjust = if_else(centered_assoc >= 0, 0, 1)
    ),
    size = 3, colour = "white"
  ) +
  facet_grid(other_category ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(
    values = c(together = "#3a82c4", apart = "#d45f3c"),
    labels = c(together = "Go together with me", apart = "Do not go together with me"),
    name   = NULL
  ) +
  scale_x_continuous(
    limits = c(-0.55, 0.55),
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
    x       = "Centered association with 'me'  (prop_together − 0.5)",
    y       = NULL,
    title   = "Self Connections — How strongly is each concept linked to 'me'?",
    subtitle = paste0(
      "Positive = rated GO TOGETHER with me  ·  Negative = DO NOT GO TOGETHER  ·  ",
      "N = ", N_participants, " participants"
    ),
    caption = paste0(
      "Only pairs with ≥ ", MIN_PARTICIPANTS_PER_PAIR,
      " participants and |assoc − 0.5| ≥ ", EDGE_THRESHOLD
    )
  )

ggsave(file.path(PLOT_DIR, "06_me_connections.pdf"),
       p6, width = 9, height = 8, device = cairo_pdf)
message("✓ Saved 06_me_connections.pdf")

# ── 14. Console summary ───────────────────────────────────────────────────────

cat("\n══════════════════════════════════════════════════════\n")
cat("  Belief Network Summary\n")
cat("══════════════════════════════════════════════════════\n\n")

cat("Nodes (concepts):\n")
node_summary %>%
  arrange(category, desc(mean_valence)) %>%
  mutate(bar = strrep("█", round(mean_valence * 20)),
         row = sprintf("  %-26s [%-8s]  %.2f  %s",
                       node, category, mean_valence, bar)) %>%
  pull(row) %>% cat(sep = "\n")

cat("\n\nTop 5 most ASSOCIATED pairs (prop_together):\n")
edge_summary %>%
  arrange(desc(mean_assoc)) %>%
  slice_head(n = 5) %>%
  mutate(row = sprintf("  %-14s ×  %-24s  %.2f",
                       node1_label, node2_label, mean_assoc)) %>%
  pull(row) %>% cat(sep = "\n")

cat("\n\nTop 5 most DISASSOCIATED pairs:\n")
edge_summary %>%
  arrange(mean_assoc) %>%
  slice_head(n = 5) %>%
  mutate(row = sprintf("  %-14s ×  %-24s  %.2f",
                       node1_label, node2_label, mean_assoc)) %>%
  pull(row) %>% cat(sep = "\n")

cat(sprintf("\n\nStructural balance correlation: r = %.3f\n", cor_r))
cat("(Negative → pairs with larger valence gaps are perceived as less associated)\n")

cat("\nNetwork density (edges above threshold):",
    nrow(edge_summary), "/",
    choose(nrow(node_summary), 2), "possible pairs\n")

cat("\n\nNode centrality (sorted by weighted strength):\n")
cat(sprintf("  %-26s  %8s  %11s  %11s\n",
            "Node", "Strength", "Betweenness", "Eigenvector"))
cat("  ", strrep("─", 62), "\n", sep = "")
node_centrality %>%
  mutate(row = sprintf("  %-26s  %8.3f  %11.3f  %11.3f",
                       name, strength_norm, betweenness, eigenvector_norm)) %>%
  pull(row) %>% cat(sep = "\n")
cat(sprintf("\n  (Saved full table to plots/node_centrality.csv)\n"))

cat("\n✓ All plots saved to:", normalizePath(PLOT_DIR), "\n")
cat("══════════════════════════════════════════════════════\n\n")
