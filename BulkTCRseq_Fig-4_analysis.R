###############################################################################
# bulk_tcrseq_analysis.R
#
# Reads TRUST4 TRB output for each bulk TCRseq well, assigns Global Clone IDs
# by matching CDR3 AA sequences against the scTCRseq-derived reference CSVs,
# then builds per-well frequency and UMI count summary tables.
###############################################################################

#setwd("~/Lab files/CMV project/TCRseq")

# ── 0. PACKAGES ───────────────────────────────────────────────────────────────
library(tidyverse)
library(compositions)


# ── 1. CONFIGURATION ──────────────────────────────────────────────────────────

# Folder for D2-D5 data (all wells imported)
main_folder_d2_d5 <- "20260320_FuncE_run_3_D2-5/TRUST4/TRUST4 output/reports"

# Folder for D1 data (only specific wells imported — see d1_wells_to_import)
main_folder_d1 <- "20251121_FuncE_run_2_D1_repeat/TRUST4/TRUST4 output/reports"

# Only these D1 wells are imported (Full- prefix only, no totals or mid)
d1_wells_to_import <- c(
  "Full-11-Control",
  "Full-10-Prolif",
  "Full-7-TCR-low",
  "Full-9-TCR-high",
  "Full-2-107a-pos",
  "Full-3-107a-neg",
  "Full-5-IFNg-pos",
  "Full-6-IFNg-neg"
)

# TCR sequence reference CSV (NT sequences, from D2-5 experiment folder)
tcr_ref_file <- "20260320_FuncE_run_3_D2-5/TRUST4/5-D_clonotypes_TCR_nt_sequences_filtered.csv"

# Minimum UMI count to retain a row (set to 1 to keep all)
min_umi <- 1

# Output directory
output_dir <- "Final_analysis_results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# ── 2. BUILD GLOBAL CLONE ID LOOKUP FROM REFERENCE CSV ───────────────────────
#
# The reference CSV has columns: Global_Clone_ID, CDR3_aa_sequence_1,
# CDR3_aa_sequence_2, ... (one column per unique CDR3 AA sequence seen
# for that clone in the scTCRseq data).
#
# We pivot these into a flat lookup: CDR3_aa_sequence -> Global_Clone_ID

tcr_ref <- read_csv(tcr_ref_file, show_col_types = FALSE)

# Gather all CDR3_nt_sequence_X columns into a single lookup table.
# Reference sequences are stored as "TRB:ATCGXX...;TRA:ATCGXX..." (with chain
# prefixes) but bulk TCRseq CDR3s are bare nucleotide sequences.
# We extract just the TRB portion from each reference sequence for matching.
clone_lookup <- tcr_ref %>%
  select(Global_Clone_ID, starts_with("CDR3_nt_sequence_")) %>%
  pivot_longer(
    cols      = starts_with("CDR3_nt_sequence_"),
    names_to  = NULL,
    values_to = "cdr3_nt_full"
  ) %>%
  filter(!is.na(cdr3_nt_full), cdr3_nt_full != "") %>%
  mutate(
    # Extract TRB sequence: grab everything after "TRB:" up to the next ";" or end
    cdr3_nt = str_match(cdr3_nt_full, "TRB:([^;]+)")[, 2]
  ) %>%
  filter(!is.na(cdr3_nt)) %>%           # drop entries with no TRB chain
  distinct(cdr3_nt, .keep_all = TRUE)   # if same TRB NT appears in 2 clones, keep first

# Named vector for fast lookup: bare TRB NT sequence -> Global_Clone_ID
cdr3_to_clone <- setNames(clone_lookup$Global_Clone_ID, clone_lookup$cdr3_nt)

message(sprintf("Reference loaded: %d unique TRB CDR3 NT sequences across %d global clones",
                length(cdr3_to_clone), n_distinct(clone_lookup$Global_Clone_ID)))


# ── 3. READ AND CLEAN ONE TRUST4 TRB FILE ─────────────────────────────────────

read_trust4_trb <- function(filepath, well_name) {
  
  df <- read.csv(filepath, stringsAsFactors = FALSE)
  
  # Columns to drop (positional — same as your original code)
  # Adjust if your TRUST4 version produces different column layouts
  cols_to_drop <- c(1, 6, 7, 10, 11, 13, 15, 17, 19, 21)
  cols_to_drop <- cols_to_drop[cols_to_drop <= ncol(df)]
  df <- df[, -cols_to_drop]
  
  # Remove rows with stop codons
  if ("stop_codon" %in% colnames(df)) {
    df <- df[df$stop_codon != "Stop Codon", ]
  }
  
  # Collapse duplicate CDR3 sequences by summing their read counts
  if ("cdr3_sequence" %in% colnames(df) && "read_count" %in% colnames(df)) {
    df <- df %>%
      group_by(cdr3_sequence) %>%
      mutate(read_count = sum(read_count, na.rm = TRUE)) %>%
      slice_head(n = 1) %>%
      ungroup()
  }
  
  # Apply minimum UMI filter
  df <- df %>% filter(read_count >= min_umi)
  
  # Calculate within-well frequency
  df <- df %>%
    mutate(
      fraction = read_count / sum(read_count, na.rm = TRUE),
      well     = well_name
    )
  
  return(df)
}


# ── 4. READ ALL WELLS ─────────────────────────────────────────────────────────

well_data_list <- list()

# ── D2-D5: import all wells found in the folder
for (subfolder in list.dirs(main_folder_d2_d5, recursive = FALSE, full.names = TRUE)) {
  well_name <- basename(subfolder)
  trb_file  <- file.path(subfolder, "TRB",
                         paste0(well_name, "_TRB.UMI_1.immune_viewer_report.csv"))
  if (!file.exists(trb_file)) {
    warning(sprintf("TRB file not found for well: %s", well_name))
    next
  }
  df <- read_trust4_trb(trb_file, well_name)
  well_data_list[[well_name]] <- df
  message(sprintf("  [D2-5] Loaded: %-30s  %d rows", well_name, nrow(df)))
}

# ── D1: import only the wells listed in d1_wells_to_import
for (subfolder in list.dirs(main_folder_d1, recursive = FALSE, full.names = TRUE)) {
  well_name <- basename(subfolder)
  if (!well_name %in% d1_wells_to_import) next
  trb_file  <- file.path(subfolder, "TRB",
                         paste0(well_name, "_TRB.UMI_1.immune_viewer_report.csv"))
  if (!file.exists(trb_file)) {
    warning(sprintf("TRB file not found for D1 well: %s", well_name))
    next
  }
  df <- read_trust4_trb(trb_file, well_name)
  well_data_list[[well_name]] <- df
  message(sprintf("  [D1]   Loaded: %-30s  %d rows", well_name, nrow(df)))
}

message(sprintf("\n%d wells loaded successfully", length(well_data_list)))


# ── 5. ASSIGN GLOBAL CLONE IDs ────────────────────────────────────────────────
#
# Match each row's CDR3 AA sequence against the lookup built in step 2.
# The CDR3 AA sequence column name from TRUST4 is cdr3_amino_acid_sequence.

well_data_list <- lapply(well_data_list, function(df) {
  df %>%
    mutate(
      Global_Clone_ID = cdr3_to_clone[cdr3_sequence]
    )
})


# ── 6. EXPORT ANNOTATED PER-WELL DATA ─────────────────────────────────────────
# One CSV per well showing all clones with their Global Clone ID (if matched)

annotated_dir <- file.path(output_dir, "annotated_wells")
dir.create(annotated_dir, showWarnings = FALSE)

for (well_name in names(well_data_list)) {
  write.csv(
    well_data_list[[well_name]],
    file      = file.path(annotated_dir, paste0(well_name, "_TRB_annotated.csv")),
    row.names = FALSE
  )
}


# ── 7. BUILD SUMMARY TABLES ACROSS ALL WELLS ─────────────────────────────────
#
# For each well, sum fraction and read_count for each Global Clone ID.
# Clones without a Global Clone ID match are excluded from the summary tables.
#
# Output: one row per Global Clone ID, one column per well.

# Combine all wells into one long table, keeping only matched clones
all_wells_long <- bind_rows(well_data_list) %>%
  filter(!is.na(Global_Clone_ID))

# Get all global clones present across any well, ordered by total read count
clone_order <- all_wells_long %>%
  group_by(Global_Clone_ID) %>%
  summarise(total_reads = sum(read_count), .groups = "drop") %>%
  arrange(desc(total_reads)) %>%
  pull(Global_Clone_ID)

# Well order: use the order files were read in
well_order <- names(well_data_list)

# ── Frequency table
freq_summary <- all_wells_long %>%
  group_by(Global_Clone_ID, well) %>%
  summarise(fraction = sum(fraction), .groups = "drop") %>%
  pivot_wider(
    names_from  = well,
    values_from = fraction,
    values_fill = 0
  ) %>%
  # Order rows by total read count, columns by input well order
  slice(match(clone_order, Global_Clone_ID)) %>%
  select(Global_Clone_ID, any_of(well_order))

# ── UMI count table
umi_summary <- all_wells_long %>%
  group_by(Global_Clone_ID, well) %>%
  summarise(read_count = sum(read_count), .groups = "drop") %>%
  pivot_wider(
    names_from  = well,
    values_from = read_count,
    values_fill = 0
  ) %>%
  slice(match(clone_order, Global_Clone_ID)) %>%
  select(Global_Clone_ID, any_of(well_order))

# Export
write.csv(freq_summary, file.path(output_dir, "TRB_clone_frequencies.csv"),  row.names = FALSE)
write.csv(umi_summary,  file.path(output_dir, "TRB_clone_UMI_counts.csv"),   row.names = FALSE)

message(sprintf(
  "\nSummary tables written:\n  %d global clones matched across %d wells",
  nrow(freq_summary), length(well_order)
))


# ── 8. CLR TRANSFORMATION ─────────────────────────────────────────────────────
#
# Centred log-ratio transform the UMI count table.
# CLR is computed per well (column) across all matched clones in that well,
# with a pseudocount added before transformation to handle zeros.

pseudocount <- 1

clr_wide <- umi_summary %>%
  mutate(across(-Global_Clone_ID, ~ {
    x  <- . + pseudocount
    p  <- x / sum(x)
    as.numeric(clr(acomp(p)))
  }))

write.csv(clr_wide, file.path(output_dir, "TRB_clone_UMI_counts_CLR.csv"), row.names = FALSE)


# ── 9. HEATMAP HELPER (defined here, called after ordering is established) ────
#
# make_freq_heatmap() builds a frequency heatmap from a UMI count table subset.
# make_clr_heatmap() builds a CLR heatmap from a UMI count table subset.
# well_order_vec  : character vector giving the desired left-to-right column order
# clone_order_vec : character vector giving the desired top-to-bottom row order
# title           : plot title string

make_freq_heatmap <- function(umi_df, well_order_vec, clone_order_vec,
                              title = "TRB clone frequency") {
  
  umi_df <- umi_df %>%
    filter(Global_Clone_ID %in% clone_order_vec) %>%
    select(Global_Clone_ID, all_of(well_order_vec))
  
  freq_df <- umi_df %>%
    mutate(across(-Global_Clone_ID, ~ . / sum(.)))
  
  freq_long <- freq_df %>%
    pivot_longer(
      cols      = -Global_Clone_ID,
      names_to  = "well",
      values_to = "Frequency"
    ) %>%
    mutate(
      Global_Clone_ID = factor(Global_Clone_ID, levels = rev(clone_order_vec)),
      well            = factor(well,            levels = well_order_vec)
    )
  
  ggplot(freq_long, aes(x = well, y = Global_Clone_ID, fill = Frequency)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_gradient(
      low  = "white",
      high = "firebrick",
      name = "Frequency"
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y     = element_text(size = 8),
      axis.title      = element_blank(),
      legend.position = "right"
    ) +
    ggtitle(title)
}

make_clr_heatmap <- function(umi_df, well_order_vec, clone_order_vec,
                             title = "CLR-transformed TRB clone abundance") {
  
  umi_df <- umi_df %>%
    filter(Global_Clone_ID %in% clone_order_vec) %>%
    select(Global_Clone_ID, all_of(well_order_vec))
  
  clr_df <- umi_df %>%
    mutate(across(-Global_Clone_ID, ~ {
      x <- . + pseudocount
      p <- x / sum(x)
      as.numeric(clr(acomp(p)))
    }))
  
  clr_long <- clr_df %>%
    pivot_longer(
      cols      = -Global_Clone_ID,
      names_to  = "well",
      values_to = "CLR"
    ) %>%
    mutate(
      Global_Clone_ID = factor(Global_Clone_ID, levels = rev(clone_order_vec)),
      well            = factor(well,            levels = well_order_vec)
    )
  
  ggplot(clr_long, aes(x = well, y = Global_Clone_ID, fill = CLR)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_gradient2(
      low      = "steelblue",
      mid      = "white",
      high     = "firebrick",
      midpoint = 0,
      name     = "CLR"
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y     = element_text(size = 8),
      axis.title      = element_blank(),
      legend.position = "right"
    ) +
    ggtitle(title)
}

# Note: both heatmap functions are called in section 12 after master ordering
# is established.

# Note: the master heatmap and per-donor heatmaps are saved in section 11
# after master_well_order and master_clone_order are defined.

message("\nDone. All outputs written to: ", output_dir)


# ── 10. DONOR-GROUPED OUTPUTS ─────────────────────────────────────────────────
#
# Define the donor groups: which wells belong to each group, and which
# Global Clone ID prefix to keep for that group's rows.

donor_groups <- list(
  
  D1_pp65 = list(
    wells        = c("Full-11-Control",
                     "Full-10-Prolif",
                     "Full-7-TCR-low",  "Full-9-TCR-high",
                     "Full-2-107a-pos", "Full-3-107a-neg",
                     "Full-5-IFNg-pos", "Full-6-IFNg-neg"),
    clone_prefix = "D1_pp65"
  ),
  
  D2_pp65 = list(
    wells        = c("1-D2-p-Control", "2-D2-p-Prolif",
                     "3-D2-p-TCR-low", "4-D2-p-TCR-high"),
    clone_prefix = "D2_pp65"
  ),
  
  D3_pp65 = list(
    wells        = c("5-D3-p-Control", "6-D3-p-Prolif",
                     "7-D3-p-TCR-low", "8-D3-p-TCR-high"),
    clone_prefix = "D3_pp65"
  ),
  
  D4_pp65 = list(
    wells        = c("9-D4-p-Control",   "10-D4-p-Prolif",
                     "11-D4-p-TCR-low",  "12-D4-p-TCR-high",
                     "13-D4-p-107a-neg", "14-D4-p-107a-pos"),
    clone_prefix = "D4_pp65"
  ),
  
  D4_IE1 = list(
    wells        = c("15-D4-I-Control",   "16-D4-I-Prolif",
                     "17-D4-I-TCR-low",   "18-D4-I-TCR-high",
                     "19-D4-I-107a-neg",  "20-D4-I-107a-pos"),
    clone_prefix = "D4_IE1"
  ),
  
  D5_IE1 = list(
    wells        = c("21-D5-I-Control", "22-D5-I-Prolif",
                     "23-D5-I-TCR-low", "24-D5-I-TCR-high"),
    clone_prefix = "D5_IE1"
  )
)

# Helper: apply pseudocount CLR per column
clr_transform <- function(df) {
  df %>%
    mutate(across(-Global_Clone_ID, ~ {
      x <- . + pseudocount
      p <- x / sum(x)
      as.numeric(clr(acomp(p)))
    }))
}

# Helper: compute frequency table from UMI count table
umi_to_freq <- function(df) {
  df %>%
    mutate(across(-Global_Clone_ID, ~ . / sum(.)))
}

grouped_dir <- file.path(output_dir, "donor_grouped")
dir.create(grouped_dir, showWarnings = FALSE)

for (group_name in names(donor_groups)) {
  grp        <- donor_groups[[group_name]]
  grp_wells  <- grp$wells
  grp_prefix <- grp$clone_prefix
  
  # Subset columns (wells) and rows (clones) for this group
  grp_wells_present <- intersect(grp_wells, colnames(umi_summary))
  
  grp_umi <- umi_summary %>%
    filter(str_starts(Global_Clone_ID, grp_prefix)) %>%
    select(Global_Clone_ID, all_of(grp_wells_present)) %>%
    mutate(clone_num = as.integer(str_extract(Global_Clone_ID, "\\d+$"))) %>%
    arrange(clone_num) %>%
    select(-clone_num)
  
  if (nrow(grp_umi) == 0) {
    warning(sprintf("No clones found for group %s — skipping", group_name))
    next
  }
  
  grp_clr  <- clr_transform(grp_umi)
  grp_freq <- umi_to_freq(grp_umi)
  
  write.csv(grp_umi,  file.path(grouped_dir, paste0(group_name, "_UMI_counts.csv")),  row.names = FALSE)
  write.csv(grp_clr,  file.path(grouped_dir, paste0(group_name, "_CLR.csv")),          row.names = FALSE)
  write.csv(grp_freq, file.path(grouped_dir, paste0(group_name, "_frequencies.csv")), row.names = FALSE)
  
  message(sprintf("  Wrote %d clones x %d wells for group: %s",
                  nrow(grp_umi), length(grp_wells_present), group_name))
}


# ── 11. MASTER FILE: ALL WELLS + ALL CLONES, DONOR-ORDERED ───────────────────
#
# Columns ordered: D2 wells, D3 wells, D4-pp65 wells, D4-IE1 wells, D5 wells
# Rows ordered:    D2 clones, D3 clones, D4-pp65 clones, D4-IE1 clones, D5 clones
#                  within each group sorted by descending total count

master_well_order <- unlist(lapply(donor_groups, `[[`, "wells"))
master_well_order <- intersect(master_well_order, colnames(umi_summary))

master_clone_order <- bind_rows(
  lapply(names(donor_groups), function(g) {
    grp <- donor_groups[[g]]
    umi_summary %>%
      filter(str_starts(Global_Clone_ID, grp$clone_prefix)) %>%
      mutate(clone_num = as.integer(str_extract(Global_Clone_ID, "\\d+$"))) %>%
      arrange(clone_num) %>%
      select(Global_Clone_ID)
  })
) %>% pull(Global_Clone_ID)

# Any clones not captured by a prefix
remaining <- setdiff(umi_summary$Global_Clone_ID, master_clone_order)
master_clone_order <- c(master_clone_order, remaining)

# Apply ordering to all three table types
apply_master_order <- function(df) {
  df %>%
    slice(match(master_clone_order, Global_Clone_ID)) %>%
    filter(!is.na(Global_Clone_ID)) %>%
    select(Global_Clone_ID, all_of(master_well_order))
}

master_umi  <- apply_master_order(umi_summary)
master_clr  <- apply_master_order(clr_wide)
master_freq <- apply_master_order(freq_summary)

# Write combined master CSV with labelled sections for UMI, CLR, and frequency
write_master_combined <- function(umi, clr_df, freq, filepath) {
  
  make_header <- function(label) {
    h <- as.data.frame(matrix("", nrow = 1, ncol = ncol(umi)))
    colnames(h) <- colnames(umi)
    h[1, 1] <- label
    h
  }
  
  blank <- make_header("")
  
  # Convert all numeric columns to character so bind_rows doesn't complain
  to_char <- function(df) mutate(df, across(everything(), as.character))
  
  combined <- bind_rows(
    make_header("--- UMI counts ---"),
    to_char(umi),
    blank,
    make_header("--- CLR-transformed counts ---"),
    to_char(clr_df),
    blank,
    make_header("--- Frequencies ---"),
    to_char(freq)
  )
  
  write.csv(combined, filepath, row.names = FALSE)
}

write_master_combined(
  master_umi, master_clr, master_freq,
  file.path(output_dir, "TRB_master_all_donors.csv")
)

message(sprintf(
  "\nMaster file written: %d clones x %d wells",
  nrow(master_umi), ncol(master_umi) - 1
))


# ── 12. HEATMAP EXPORTS ───────────────────────────────────────────────────────

freq_heatmap_dir <- file.path(output_dir, "heatmaps", "frequency")
clr_heatmap_dir  <- file.path(output_dir, "heatmaps", "CLR")
dir.create(freq_heatmap_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(clr_heatmap_dir,  showWarnings = FALSE, recursive = TRUE)

save_both_heatmaps <- function(umi_df, well_order_vec, clone_order_vec,
                               base_filename, title_suffix) {
  w <- max(8,  length(well_order_vec)  * 0.5 + 2)
  h <- max(6,  length(clone_order_vec) * 0.3 + 1.5)
  
  ggsave(
    filename = file.path(freq_heatmap_dir, paste0(base_filename, "_frequency.png")),
    plot     = make_freq_heatmap(umi_df, well_order_vec, clone_order_vec,
                                 title = paste0("TRB clone frequency — ", title_suffix)),
    width = w, height = h, dpi = 300
  )
  
  ggsave(
    filename = file.path(clr_heatmap_dir, paste0(base_filename, "_CLR.png")),
    plot     = make_clr_heatmap(umi_df, well_order_vec, clone_order_vec,
                                title = paste0("CLR-transformed TRB clone abundance — ", title_suffix)),
    width = w, height = h, dpi = 300
  )
}

# ── Master heatmap (all donors, master ordering)
save_both_heatmaps(
  umi_df          = umi_summary,
  well_order_vec  = master_well_order,
  clone_order_vec = master_clone_order,
  base_filename   = "TRB_heatmap_all_donors",
  title_suffix    = "all donors"
)
message("  Saved master heatmaps")

# ── Per-donor/specificity heatmaps
for (group_name in names(donor_groups)) {
  grp        <- donor_groups[[group_name]]
  grp_wells  <- intersect(grp$wells, colnames(umi_summary))
  grp_clones <- master_clone_order[str_starts(master_clone_order, grp$clone_prefix)]
  
  if (length(grp_clones) == 0 || length(grp_wells) == 0) next
  
  save_both_heatmaps(
    umi_df          = umi_summary,
    well_order_vec  = grp_wells,
    clone_order_vec = grp_clones,
    base_filename   = paste0("TRB_heatmap_", group_name),
    title_suffix    = group_name
  )
  message(sprintf("  Saved heatmaps for %s", group_name))
}

message("\nAll done.")


# ── 13. PER-ASSAY CLR BAR GRAPHS ─────────────────────────────────────────────
#
# One PNG per assay per donor/specificity group.
# Each assay plot includes the Control well as the first bar, followed by the
# assay-specific wells. Each clone gets its own facet panel.
#
# Assay definitions (pattern matched against well names):
#   Proliferation : Control + Prolif
#   TCR level     : Control + TCR-low + TCR-high
#   CD107a        : Control + 107a-neg + 107a-pos  (only for D4)

bar_dir <- file.path(output_dir, "barplots")
dir.create(bar_dir, showWarnings = FALSE)

# ── Assay definitions: name -> regex pattern for non-control wells
assay_defs <- list(
  Proliferation = "Prolif",
  TCR_level     = "TCR",
  CD107a        = "107a",
  IFNg          = "IFNg"
)

# ── Well fill colours per assay (Control always grey, assay wells get a colour)
assay_well_colors <- list(
  Proliferation = c("Control" = "#AAAAAA", "Prolif"    = "#4E9AC7"),
  TCR_level     = c("Control" = "#AAAAAA", "TCR-low"   = "#F5C176", "TCR-high" = "#E07B3A"),
  CD107a        = c("Control" = "#AAAAAA", "107a-neg"  = "#A8D5A2", "107a-pos" = "#57A65A"),
  IFNg          = c("Control" = "#AAAAAA", "IFNg-neg"  = "#D4A8D5", "IFNg-pos" = "#9B4FA6")
)

# ── Short well labels: handles both D2-5 format ("1-D2-p-Control")
# and D1 format ("Full-11-Control")
short_label <- function(well_name) {
  well_name %>%
    str_remove("^\\d+-D\\d+-[pI]-") %>%    # D2-5 format: "1-D2-p-"
    str_remove("^(Full|Half)-\\d+-")        # D1 format:   "Full-11-"
}

# ── Build and save one bar plot for a given set of wells and clones
make_assay_bar_plot <- function(umi_df, assay_wells, grp_clones,
                                color_map, title) {
  
  # Subset rows (clones) and columns (wells for this assay incl. Control)
  grp_umi <- umi_df %>%
    filter(Global_Clone_ID %in% grp_clones) %>%
    select(Global_Clone_ID, all_of(assay_wells))
  
  # CLR transform per well independently (each well normalised to itself only)
  grp_clr <- grp_umi %>%
    mutate(across(-Global_Clone_ID, ~ {
      x <- . + pseudocount
      as.numeric(clr(acomp(x / sum(x))))
    }))
  
  # Well display metadata
  well_meta <- tibble(
    well      = assay_wells,
    label     = short_label(assay_wells),
    x_pos     = seq_along(assay_wells)
  )
  
  # Pivot to long and join metadata
  plot_df <- grp_clr %>%
    mutate(Global_Clone_ID = factor(Global_Clone_ID, levels = grp_clones)) %>%
    pivot_longer(
      cols      = -Global_Clone_ID,
      names_to  = "well",
      values_to = "CLR"
    ) %>%
    left_join(well_meta, by = "well") %>%
    mutate(label = factor(label, levels = well_meta$label))
  
  y_lim <- max(abs(plot_df$CLR), na.rm = TRUE) * 1.15
  
  # Extend color_map to cover any labels not explicitly listed (fallback grey)
  all_labels   <- well_meta$label
  named_colors <- color_map
  missing      <- setdiff(all_labels, names(named_colors))
  if (length(missing) > 0) {
    named_colors[missing] <- "#CCCCCC"
  }
  
  ggplot(plot_df, aes(x = label, y = CLR, fill = label)) +
    geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
    facet_wrap(~ Global_Clone_ID, nrow = 1) +
    scale_fill_manual(values = named_colors, guide = "none") +
    scale_y_continuous(
      limits = c(-y_lim, y_lim),
      expand = expansion(mult = 0.05)
    ) +
    labs(
      title = title,
      y     = "CLR value",
      x     = NULL
    ) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x      = element_text(size = 8, angle = 45, hjust = 1),
      axis.text.y      = element_text(size = 7),
      strip.text       = element_text(size = 8, face = "bold"),
      strip.background = element_blank(),
      panel.spacing    = unit(0.8, "lines"),
      plot.title       = element_text(size = 11, face = "bold")
    )
}

# ── Loop: one plot per donor group × assay
for (group_name in names(donor_groups)) {
  grp        <- donor_groups[[group_name]]
  grp_wells  <- intersect(grp$wells, colnames(umi_summary))
  grp_clones <- master_clone_order[str_starts(master_clone_order, grp$clone_prefix)]
  
  if (length(grp_clones) == 0 || length(grp_wells) == 0) next
  
  # Identify the Control well for this group (there should be exactly one)
  control_well <- grp_wells[str_detect(grp_wells, "Control")]
  
  for (assay_name in names(assay_defs)) {
    pattern      <- assay_defs[[assay_name]]
    assay_wells  <- grp_wells[str_detect(grp_wells, pattern)]
    
    # Skip if this assay has no wells in this donor group (e.g. no CD107a for D2)
    if (length(assay_wells) == 0) next
    
    # Prepend the Control well (Control + assay-specific wells in order)
    plot_wells <- c(control_well, assay_wells)
    
    color_map <- assay_well_colors[[assay_name]]
    title     <- paste0("CLR — ", group_name, " — ", assay_name)
    
    p <- make_assay_bar_plot(
      umi_df      = umi_summary,
      assay_wells = plot_wells,
      grp_clones  = grp_clones,
      color_map   = color_map,
      title       = title
    )
    
    n_rows   <- 1
    n_wells  <- length(plot_wells)
    
    ggsave(
      filename = file.path(bar_dir,
                           paste0("TRB_CLR_barplot_", group_name, "_", assay_name, ".png")),
      plot     = p,
      width    = max(5, length(grp_clones) * n_wells * 0.5 + 1.5),
      height   = 3.5,
      dpi      = 300
    )
    message(sprintf("  Saved: %s — %s  (%d clones, %d wells)",
                    group_name, assay_name, length(grp_clones), n_wells))
  }
}

message("\nBar plots done.")


# ── 14. MASTER ASSAY SCORES ───────────────────────────────────────────────────
#
# For each clone × donor group, compute pairwise assay scores using CLR values
# recomputed fresh from only the two wells being contrasted. Using per-assay
# CLR rather than global CLR ensures each score is a pure pairwise contrast
# uninfluenced by other wells.
#
# Scores defined (CLR recomputed across the relevant well pair each time):
#   Proliferation score : CLR(Prolif)   - CLR(Control)
#   TCR score           : CLR(TCR-high) - CLR(TCR-low)
#   CD107a score        : CLR(107a-pos) - CLR(107a-neg)   [D1, D4 only]
#   IFNg score          : CLR(IFNg-pos) - CLR(IFNg-neg)   [D1 only]
#
# Also computed: frequency-weighted versions of each score, where the raw score
# is multiplied by the clone's mean frequency across both wells of the contrast.
# This downweights scores from very rare clones where the ratio is noisy.

# Helper: given a UMI data frame, compute CLR for two specified wells and
# return their difference (pos - neg) plus a frequency-weighted version.
pairwise_clr_diff <- function(umi_df, well_pos, well_neg) {
  
  if (!all(c(well_pos, well_neg) %in% colnames(umi_df))) return(NULL)
  
  pair <- umi_df %>%
    select(Global_Clone_ID, all_of(c(well_pos, well_neg)))
  
  # CLR transform per well independently (each well normalised to itself only)
  clr_pair <- pair %>%
    mutate(across(-Global_Clone_ID, ~ {
      x <- . + pseudocount
      as.numeric(clr(acomp(x / sum(x))))
    }))
  
  raw_score <- clr_pair[[well_pos]] - clr_pair[[well_neg]]
  
  freq_pair <- pair %>%
    mutate(across(-Global_Clone_ID, ~ . / sum(.)))
  mean_freq <- (freq_pair[[well_pos]] + freq_pair[[well_neg]]) / 2
  
  tibble(
    Global_Clone_ID = clr_pair$Global_Clone_ID,
    raw_score       = raw_score,
    weighted_score  = raw_score * mean_freq
  )
}

# Helper: find the actual well name in a vector that matches a short pattern
find_well <- function(wells, pattern) {
  matched <- wells[str_detect(wells, pattern)]
  if (length(matched) == 0) return(NA_character_)
  matched[1]
}

# Define contrasts per donor group
group_contrasts <- list(
  D1_pp65 = list(
    list(pos = "Prolif",   neg = "Control", name = "Proliferation"),
    list(pos = "TCR-high", neg = "TCR-low", name = "TCR"),
    list(pos = "107a-pos", neg = "107a-neg", name = "CD107a"),
    list(pos = "IFNg-pos", neg = "IFNg-neg", name = "IFNg")
  ),
  D2_pp65 = list(
    list(pos = "Prolif",   neg = "Control", name = "Proliferation"),
    list(pos = "TCR-high", neg = "TCR-low", name = "TCR")
  ),
  D3_pp65 = list(
    list(pos = "Prolif",   neg = "Control", name = "Proliferation"),
    list(pos = "TCR-high", neg = "TCR-low", name = "TCR")
  ),
  D4_pp65 = list(
    list(pos = "Prolif",   neg = "Control", name = "Proliferation"),
    list(pos = "TCR-high", neg = "TCR-low", name = "TCR"),
    list(pos = "107a-pos", neg = "107a-neg", name = "CD107a")
  ),
  D4_IE1 = list(
    list(pos = "Prolif",   neg = "Control", name = "Proliferation"),
    list(pos = "TCR-high", neg = "TCR-low", name = "TCR"),
    list(pos = "107a-pos", neg = "107a-neg", name = "CD107a")
  ),
  D5_IE1 = list(
    list(pos = "Prolif",   neg = "Control", name = "Proliferation"),
    list(pos = "TCR-high", neg = "TCR-low", name = "TCR")
  )
)

scores_list <- list()

for (group_name in names(group_contrasts)) {
  grp        <- donor_groups[[group_name]]
  grp_wells  <- intersect(grp$wells, colnames(umi_summary))
  grp_clones <- master_clone_order[str_starts(master_clone_order, grp$clone_prefix)]
  
  if (length(grp_clones) == 0) next
  
  grp_umi <- umi_summary %>%
    filter(Global_Clone_ID %in% grp_clones) %>%
    select(Global_Clone_ID, all_of(grp_wells))
  
  score_df <- tibble(
    Global_Clone_ID = grp_clones,
    Donor_Group     = group_name
  )
  
  for (contrast in group_contrasts[[group_name]]) {
    well_pos <- find_well(grp_wells, contrast$pos)
    well_neg <- find_well(grp_wells, contrast$neg)
    
    if (is.na(well_pos) || is.na(well_neg)) {
      message(sprintf("  Skipping %s for %s (wells not found)", contrast$name, group_name))
      next
    }
    
    result <- pairwise_clr_diff(grp_umi, well_pos, well_neg)
    
    score_df <- score_df %>%
      left_join(
        result %>% select(
          Global_Clone_ID,
          !!paste0(contrast$name, "_score")          := raw_score,
          !!paste0(contrast$name, "_score_weighted") := weighted_score
        ),
        by = "Global_Clone_ID"
      )
  }
  
  scores_list[[group_name]] <- score_df
  message(sprintf("  Computed assay scores for %s", group_name))
}

master_assay_scores <- bind_rows(scores_list) %>%
  relocate(Donor_Group, .before = everything())

# Drop weighted score columns — global clones are already size-filtered
master_assay_scores <- master_assay_scores %>%
  select(-ends_with("_score_weighted"))

# Add within-donor range-based percentiles for each score column.
# Percentile = (score - min) / (max - min) * 100
# This reflects where each clone falls on the actual score spectrum for that
# donor group, rather than just its rank position among clones.
# NAs (assays not run for that donor) are left as NA.
score_cols <- names(master_assay_scores) %>%
  str_subset("_score$")

master_assay_scores <- master_assay_scores %>%
  group_by(Donor_Group) %>%
  mutate(across(
    all_of(score_cols),
    ~ {
      mn <- min(., na.rm = TRUE)
      mx <- max(., na.rm = TRUE)
      if (mx == mn) NA_real_ else (. - mn) / (mx - mn) * 100
    },
    .names = "{.col}_percentile"
  )) %>%
  ungroup() %>%
  # Interleave: raw score then its percentile, for each assay
  relocate(Global_Clone_ID, .before = everything()) %>%
  relocate(Donor_Group, .after = Global_Clone_ID) %>%
  relocate(
    unlist(lapply(score_cols, function(s) c(s, paste0(s, "_percentile")))),
    .after = Donor_Group
  )

write.csv(
  master_assay_scores,
  file.path(output_dir, "Master_Assay_Scores.csv"),
  row.names = FALSE,
  na        = ""
)

message(sprintf(
  "\nMaster_Assay_Scores.csv written: %d clones across %d donor groups",
  nrow(master_assay_scores),
  n_distinct(master_assay_scores$Donor_Group)
))


# ── 15. CORRELATION PLOTS ─────────────────────────────────────────────────────
#
# For each donor group, produce:
#   (a) Scatter plots of TCR_score_percentile (x) vs each other assay score
#       percentile (y), one plot per available assay pairing.
#   (b) Scatter plot of TCR_score_percentile (x) vs Control well frequency (y).
#
# Each plot shows individual clones as labelled points, a linear regression
# line with 95% confidence ribbon, and an annotation with R² and p-value
# (p-value tests whether slope differs from zero).

corr_dir <- file.path(output_dir, "correlation_plots")
dir.create(corr_dir, showWarnings = FALSE)

# Helper: build one scatter + regression plot
make_corr_plot <- function(df, x_col, y_col, x_label, y_label, title) {
  
  plot_df <- df %>%
    select(Global_Clone_ID, x = all_of(x_col), y = all_of(y_col)) %>%
    filter(!is.na(x), !is.na(y))
  
  if (nrow(plot_df) < 3) {
    message(sprintf("  Skipping %s — fewer than 3 complete observations", title))
    return(NULL)
  }
  
  fit   <- lm(y ~ x, data = plot_df)
  r2    <- summary(fit)$r.squared
  pval  <- tryCatch(
    summary(fit)$coefficients["x", "Pr(>|t|)"],
    error = function(e) NA_real_
  )
  pval_text <- if (is.na(pval))  "p = NA" else
    if (pval < 0.001) "p < 0.001" else
      sprintf("p = %.3f", pval)
  annot     <- sprintf("R² = %.3f\n%s", r2, pval_text)
  
  ggplot(plot_df, aes(x = x, y = y)) +
    geom_smooth(method = "lm", se = FALSE,
                colour = "steelblue",
                linewidth = 0.8) +
    geom_point(size = 2.5, colour = "grey30") +
    ggrepel::geom_text_repel(
      aes(label = str_replace(Global_Clone_ID, ".*-clone-(\\d+)$", "Clone \\1")),
      size          = 2.5,
      max.overlaps  = 20,
      segment.size  = 0.3,
      segment.color = "grey60"
    ) +
    annotate("text",
             x = -Inf, y = Inf,
             hjust = -0.1, vjust = 1.4,
             label = annot,
             size  = 3.2,
             colour = "grey20") +
    labs(title = title, x = x_label, y = y_label) +
    theme_classic(base_size = 11) +
    theme(
      plot.title  = element_text(size = 10, face = "bold"),
      axis.title  = element_text(size = 9)
    )
}

# ── Build control-well frequency lookup: Global_Clone_ID -> frequency in Control
# Find control well name for each donor group from the well list
control_freq_lookup <- bind_rows(
  lapply(names(donor_groups), function(g) {
    grp         <- donor_groups[[g]]
    ctrl_well   <- intersect(grp$wells, colnames(umi_summary)) %>%
      .[str_detect(., "Control")] %>% .[1]
    if (is.na(ctrl_well)) return(NULL)
    
    umi_summary %>%
      select(Global_Clone_ID, count = all_of(ctrl_well)) %>%
      mutate(
        freq        = count / sum(count),
        Donor_Group = g,
        ctrl_well   = ctrl_well
      ) %>%
      select(Donor_Group, Global_Clone_ID, ctrl_freq = freq)
  })
)

# ── Loop over donor groups and generate plots
for (group_name in names(donor_groups)) {
  
  grp_scores <- master_assay_scores %>%
    filter(Donor_Group == group_name)
  
  if (nrow(grp_scores) == 0) next
  
  # Percentile columns available for this group (exclude TCR itself and any
  # columns that are all NA for this group)
  pct_cols <- names(grp_scores) %>%
    str_subset("_score_percentile$") %>%
    setdiff("TCR_score_percentile") %>%
    keep(~ !all(is.na(grp_scores[[.x]])))
  
  # (a) TCR percentile vs other assay percentiles
  for (y_col in pct_cols) {
    assay_label <- y_col %>%
      str_remove("_score_percentile$") %>%
      str_replace("_", " ")
    
    p <- make_corr_plot(
      df      = grp_scores,
      x_col   = "TCR_score_percentile",
      y_col   = y_col,
      x_label = "TCR score percentile",
      y_label = paste0(assay_label, " score percentile"),
      title   = paste0(group_name, "  |  TCR vs ", assay_label)
    )
    
    if (!is.null(p)) {
      ggsave(
        filename = file.path(corr_dir,
                             paste0(group_name, "_TCR_vs_", assay_label, ".png")),
        plot   = p, width = 5, height = 4.5, dpi = 300
      )
      message(sprintf("  Saved: %s — TCR vs %s", group_name, assay_label))
    }
  }
  
  # (b) TCR percentile vs Control well frequency
  ctrl_freq <- control_freq_lookup %>%
    filter(Donor_Group == group_name)
  
  plot_df_freq <- grp_scores %>%
    left_join(ctrl_freq, by = c("Donor_Group", "Global_Clone_ID"))
  
  if (!all(is.na(plot_df_freq$ctrl_freq))) {
    p_freq <- make_corr_plot(
      df      = plot_df_freq,
      x_col   = "TCR_score_percentile",
      y_col   = "ctrl_freq",
      x_label = "TCR score percentile",
      y_label = "Control well frequency",
      title   = paste0(group_name, "  |  TCR vs Control frequency")
    )
    
    if (!is.null(p_freq)) {
      ggsave(
        filename = file.path(corr_dir,
                             paste0(group_name, "_TCR_vs_Control_frequency.png")),
        plot   = p_freq, width = 5, height = 4.5, dpi = 300
      )
      message(sprintf("  Saved: %s — TCR vs Control frequency", group_name))
    }
  }
}

message("\nCorrelation plots done.")


# ── 16. METACLUSTER COMPOSITION CORRELATION PLOTS ────────────────────────────
#
# Reads the metacluster composition CSV exported from the Seurat notebook and
# plots TCR_score_percentile (x) vs percentage of cells in each metacluster (y),
# one plot per metacluster per donor group.
# Linear regression with R² and p-value annotation, clone labels as "Clone N".

# Path to the CSV exported from the Seurat notebook
metacluster_csv <- "Global_Clone_metacluster_composition.csv"

metacluster_comp <- read_csv(metacluster_csv, show_col_types = FALSE)

metacluster_levels <- c("Tcm", "Tem", "T-CMC1", "NK-like", "ISG+")
pct_cols           <- paste0("pct_", metacluster_levels)

meta_corr_dir <- file.path(output_dir, "correlation_plots", "metacluster_composition")
dir.create(meta_corr_dir, showWarnings = FALSE, recursive = TRUE)

for (group_name in names(donor_groups)) {
  
  grp_scores <- master_assay_scores %>%
    filter(Donor_Group == group_name) %>%
    select(Global_Clone_ID, TCR_score_percentile)
  
  if (nrow(grp_scores) == 0 || all(is.na(grp_scores$TCR_score_percentile))) next
  
  # Join TCR percentile with metacluster composition
  plot_base <- grp_scores %>%
    left_join(metacluster_comp, by = "Global_Clone_ID") %>%
    filter(!is.na(TCR_score_percentile))
  
  # One plot per metacluster
  for (pct_col in pct_cols) {
    
    cluster_name <- str_remove(pct_col, "^pct_")
    
    plot_df <- plot_base %>%
      select(Global_Clone_ID,
             x = TCR_score_percentile,
             y = all_of(pct_col)) %>%
      filter(!is.na(y))
    
    if (nrow(plot_df) < 3) next
    
    fit       <- lm(y ~ x, data = plot_df)
    r2        <- summary(fit)$r.squared
    pval      <- tryCatch(
      summary(fit)$coefficients["x", "Pr(>|t|)"],
      error = function(e) NA_real_
    )
    pval_text <- if (is.na(pval))  "p = NA" else
      if (pval < 0.001) "p < 0.001" else
        sprintf("p = %.3f", pval)
    annot     <- sprintf("R² = %.3f\n%s", r2, pval_text)
    
    p <- ggplot(plot_df, aes(x = x, y = y)) +
      geom_smooth(method = "lm", se = FALSE,
                  colour = "steelblue", linewidth = 0.8) +
      geom_point(size = 2.5, colour = "grey30") +
      ggrepel::geom_text_repel(
        aes(label = str_replace(Global_Clone_ID, ".*-clone-(\\d+)$", "Clone \\1")),
        size          = 2.5,
        max.overlaps  = 20,
        segment.size  = 0.3,
        segment.color = "grey60"
      ) +
      annotate("text",
               x = -Inf, y = Inf,
               hjust = -0.1, vjust = 1.4,
               label = annot,
               size  = 3.2,
               colour = "grey20") +
      labs(
        title = paste0(group_name, "  |  TCR vs % ", cluster_name),
        x     = "TCR score percentile",
        y     = paste0("% cells in ", cluster_name)
      ) +
      theme_classic(base_size = 11) +
      theme(
        plot.title = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 9)
      )
    
    ggsave(
      filename = file.path(meta_corr_dir,
                           paste0(group_name, "_TCR_vs_pct_", cluster_name, ".png")),
      plot   = p, width = 5, height = 4.5, dpi = 300
    )
    message(sprintf("  Saved: %s — TCR vs pct %s", group_name, cluster_name))
  }
}




###############################################################################
# correlation_plots_for_figure.R
#
# Regenerates the correlation plots from bulk_tcrseq_analysis.R (sections 15-16)
# in a figure-ready form.
#
# Differences from the original:
#   * Starts from the exported UMI count table + metacluster composition CSV
#     instead of re-reading TRUST4 output
#   * Writes to Bulk_TCR_seq/Final_analysis_results/correlation_plots_for_figure/
#   * Clones are pooled by specificity (pp65, IE1) rather than one set of plots
#     per donor group; points are coloured by donor
#   * Assay axes use the raw CLR-difference scores (set use_percentiles = TRUE
#     to switch back to the within-donor percentiles)
#   * Percentage axes fixed 0-100; raw score axes are shared across all panels
#     of a given specificity
#   * Adds CD107a score vs metacluster composition plots (donors with CD107a
#     data only - detected automatically)
#   * Title, axis title, axis tick label and R2 / p annotation font sizes are
#     1.2x the original; clone point labels are unchanged
###############################################################################

# ── 0. PACKAGES ───────────────────────────────────────────────────────────────
library(tidyverse)
library(compositions)
library(ggrepel)


# ── 1. CONFIGURATION ──────────────────────────────────────────────────────────

results_dir     <- "Bulk_TCR_seq/Final_analysis_results"
umi_file        <- file.path(results_dir, "TRB_clone_UMI_counts.csv")
metacluster_csv <- "Bulk_TCR_seq/Global_Clone_metacluster_composition.csv"

output_dir <- file.path(results_dir, "correlation_plots_for_figure")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

pseudocount <- 1

# ── Assay axis units: raw CLR-difference scores, or within-donor percentiles
use_percentiles <- FALSE

# ── Fixed axis limits
pct_limits        <- c(0, 100)    # metacluster percentages and clone dominance
percentile_limits <- c(0, 100)    # assay score percentiles (if used)

# Raw scores have no natural bounds. TRUE = one shared range per score across
# all panels of a specificity, so panels can be compared side by side.
shared_score_limits <- TRUE

# ── Font sizes: 1.2x the originals
font_scale       <- 1.2
title_size       <- 10  * font_scale   # was 10
axis_title_size  <- 9   * font_scale   # was 9
axis_text_size   <- 8.8 * font_scale   # was theme_classic(11) default (~8.8)
annot_size       <- 3.2 * font_scale   # R2 / p annotation, was 3.2
clone_label_size <- 2.5                # clone point labels - original size, not scaled

# ── Canvas (inches). Wider than the original 5 to make room for the right-hand
# legend, so the plotting panel itself stays close to its original size.
plot_width <- 6

# y_axis_stretch lengthens the panel (i.e. the y axis itself) rather than the
# whole canvas: non_panel_height is the approximate vertical space taken by the
# title, the x-axis text and title, and the margins, none of which stretch.
plot_height_base <- 4.5
y_axis_stretch   <- 1.2
non_panel_height <- 1.1
plot_height <- non_panel_height + (plot_height_base - non_panel_height) * y_axis_stretch

plot_dpi <- 300

# ── Point colours, fixed per donor so a given donor is the same colour in
# every panel regardless of which donors have data for that contrast
donor_colours <- c(
  D1 = "#4E9AC7",
  D2 = "#E07B3A",
  D3 = "#57A65A",
  D4 = "#9B4FA6",
  D5 = "#D9B23A"
)

# ── Options
se_ribbon <- FALSE   # TRUE = draw the 95% CI ribbon around the lm line


# ── 2. DONOR GROUPS AND CONTRASTS ─────────────────────────────────────────────

donor_groups <- list(
  D1_pp65 = list(
    wells = c("Full-11-Control", "Full-10-Prolif",
              "Full-7-TCR-low",  "Full-9-TCR-high",
              "Full-2-107a-pos", "Full-3-107a-neg",
              "Full-5-IFNg-pos", "Full-6-IFNg-neg"),
    clone_prefix = "D1_pp65"
  ),
  D2_pp65 = list(
    wells = c("1-D2-p-Control", "2-D2-p-Prolif",
              "3-D2-p-TCR-low", "4-D2-p-TCR-high"),
    clone_prefix = "D2_pp65"
  ),
  D3_pp65 = list(
    wells = c("5-D3-p-Control", "6-D3-p-Prolif",
              "7-D3-p-TCR-low", "8-D3-p-TCR-high"),
    clone_prefix = "D3_pp65"
  ),
  D4_pp65 = list(
    wells = c("9-D4-p-Control",   "10-D4-p-Prolif",
              "11-D4-p-TCR-low",  "12-D4-p-TCR-high",
              "13-D4-p-107a-neg", "14-D4-p-107a-pos"),
    clone_prefix = "D4_pp65"
  ),
  D4_IE1 = list(
    wells = c("15-D4-I-Control",  "16-D4-I-Prolif",
              "17-D4-I-TCR-low",  "18-D4-I-TCR-high",
              "19-D4-I-107a-neg", "20-D4-I-107a-pos"),
    clone_prefix = "D4_IE1"
  ),
  D5_IE1 = list(
    wells = c("21-D5-I-Control", "22-D5-I-Prolif",
              "23-D5-I-TCR-low", "24-D5-I-TCR-high"),
    clone_prefix = "D5_IE1"
  )
)

group_contrasts <- list(
  D1_pp65 = list(
    list(pos = "Prolif",   neg = "Control",   name = "Proliferation"),
    list(pos = "TCR-high", neg = "TCR-low",   name = "TCR"),
    list(pos = "107a-pos", neg = "107a-neg",  name = "CD107a"),
    list(pos = "IFNg-pos", neg = "IFNg-neg",  name = "IFNg")
  ),
  D2_pp65 = list(
    list(pos = "Prolif",   neg = "Control",   name = "Proliferation"),
    list(pos = "TCR-high", neg = "TCR-low",   name = "TCR")
  ),
  D3_pp65 = list(
    list(pos = "Prolif",   neg = "Control",   name = "Proliferation"),
    list(pos = "TCR-high", neg = "TCR-low",   name = "TCR")
  ),
  D4_pp65 = list(
    list(pos = "Prolif",   neg = "Control",   name = "Proliferation"),
    list(pos = "TCR-high", neg = "TCR-low",   name = "TCR"),
    list(pos = "107a-pos", neg = "107a-neg",  name = "CD107a")
  ),
  D4_IE1 = list(
    list(pos = "Prolif",   neg = "Control",   name = "Proliferation"),
    list(pos = "TCR-high", neg = "TCR-low",   name = "TCR"),
    list(pos = "107a-pos", neg = "107a-neg",  name = "CD107a")
  ),
  D5_IE1 = list(
    list(pos = "Prolif",   neg = "Control",   name = "Proliferation"),
    list(pos = "TCR-high", neg = "TCR-low",   name = "TCR")
  )
)

# ── Pooled plotting groups: all clones of a given specificity together
combined_groups <- list(
  pp65 = list(label = "pp65", donors = c("D1_pp65", "D2_pp65", "D3_pp65", "D4_pp65")),
  IE1  = list(label = "IE1",  donors = c("D4_IE1", "D5_IE1"))
)

metacluster_levels <- c("Tcm", "Tem", "T-CMC1", "NK-like", "ISG+")


# ── 3. IMPORT EXPORTED TABLES ─────────────────────────────────────────────────
# read_csv (not read.csv) so well names like "1-D2-p-Control" are not mangled.

umi_summary <- read_csv(umi_file, show_col_types = FALSE)

message(sprintf("UMI table loaded: %d clones x %d wells",
                nrow(umi_summary), ncol(umi_summary) - 1))

if (file.exists(metacluster_csv)) {
  metacluster_comp <- read_csv(metacluster_csv, show_col_types = FALSE)
  pct_cols <- intersect(paste0("pct_", metacluster_levels), names(metacluster_comp))
  metacluster_comp <- metacluster_comp %>% select(Global_Clone_ID, all_of(pct_cols))
  message(sprintf("Metacluster composition loaded: %d clones, %d metaclusters",
                  nrow(metacluster_comp), length(pct_cols)))
} else {
  metacluster_comp <- tibble(Global_Clone_ID = character())
  pct_cols <- character()
  warning("Metacluster composition CSV not found - metacluster plots will be skipped")
}


# ── 4. ASSAY SCORES AND WITHIN-DONOR PERCENTILES ─────────────────────────────
# Identical maths to section 14 of the original script.

clones_for_group <- function(prefix) {
  umi_summary %>%
    filter(str_starts(Global_Clone_ID, prefix)) %>%
    mutate(clone_num = as.integer(str_extract(Global_Clone_ID, "\\d+$"))) %>%
    arrange(clone_num) %>%
    pull(Global_Clone_ID)
}

find_well <- function(wells, pattern) {
  matched <- wells[str_detect(wells, pattern)]
  if (length(matched) == 0) return(NA_character_)
  matched[1]
}

pairwise_clr_diff <- function(umi_df, well_pos, well_neg) {
  if (!all(c(well_pos, well_neg) %in% colnames(umi_df))) return(NULL)

  pair <- umi_df %>% select(Global_Clone_ID, all_of(c(well_pos, well_neg)))

  # CLR per well independently, so each score is a pure pairwise contrast
  clr_pair <- pair %>%
    mutate(across(-Global_Clone_ID, ~ {
      x <- . + pseudocount
      as.numeric(clr(acomp(x / sum(x))))
    }))

  tibble(
    Global_Clone_ID = clr_pair$Global_Clone_ID,
    raw_score       = clr_pair[[well_pos]] - clr_pair[[well_neg]]
  )
}

scores_list <- list()

for (group_name in names(group_contrasts)) {
  grp        <- donor_groups[[group_name]]
  grp_wells  <- intersect(grp$wells, colnames(umi_summary))
  grp_clones <- clones_for_group(grp$clone_prefix)

  if (length(grp_clones) == 0 || length(grp_wells) == 0) {
    warning(sprintf("No clones or wells found for %s - skipping", group_name))
    next
  }

  grp_umi <- umi_summary %>%
    filter(Global_Clone_ID %in% grp_clones) %>%
    select(Global_Clone_ID, all_of(grp_wells))

  score_df <- tibble(Global_Clone_ID = grp_clones, Donor_Group = group_name)

  for (contrast in group_contrasts[[group_name]]) {
    well_pos <- find_well(grp_wells, contrast$pos)
    well_neg <- find_well(grp_wells, contrast$neg)

    if (is.na(well_pos) || is.na(well_neg)) {
      message(sprintf("  Skipping %s for %s (wells not found)",
                      contrast$name, group_name))
      next
    }

    result <- pairwise_clr_diff(grp_umi, well_pos, well_neg)

    score_df <- score_df %>%
      left_join(
        result %>% select(Global_Clone_ID,
                          !!paste0(contrast$name, "_score") := raw_score),
        by = "Global_Clone_ID"
      )
  }

  scores_list[[group_name]] <- score_df
  message(sprintf("  Computed assay scores for %s", group_name))
}

master_assay_scores <- bind_rows(scores_list)

# Within-donor range-based percentile: (score - min) / (max - min) * 100
# Kept so use_percentiles can be flipped without touching anything else.
score_cols <- names(master_assay_scores) %>% str_subset("_score$")

master_assay_scores <- master_assay_scores %>%
  group_by(Donor_Group) %>%
  mutate(across(
    all_of(score_cols),
    ~ {
      mn <- suppressWarnings(min(., na.rm = TRUE))
      mx <- suppressWarnings(max(., na.rm = TRUE))
      if (!is.finite(mn) || !is.finite(mx) || mx == mn) NA_real_
      else (. - mn) / (mx - mn) * 100
    },
    .names = "{.col}_percentile"
  )) %>%
  ungroup()


# ── 5. CONTROL-WELL FREQUENCY LOOKUP ─────────────────────────────────────────

control_freq_lookup <- bind_rows(
  lapply(names(donor_groups), function(g) {
    grp       <- donor_groups[[g]]
    wells     <- intersect(grp$wells, colnames(umi_summary))
    ctrl_well <- wells[str_detect(wells, "Control")][1]
    if (is.na(ctrl_well)) return(NULL)

    umi_summary %>%
      select(Global_Clone_ID, count = all_of(ctrl_well)) %>%
      mutate(ctrl_freq = count / sum(count), Donor_Group = g) %>%
      select(Donor_Group, Global_Clone_ID, ctrl_freq)
  })
)


# ── 6. ASSEMBLE ONE TABLE PER CLONE ──────────────────────────────────────────

clone_table <- master_assay_scores %>%
  left_join(control_freq_lookup, by = c("Donor_Group", "Global_Clone_ID")) %>%
  left_join(metacluster_comp,    by = "Global_Clone_ID") %>%
  mutate(
    ctrl_pct    = ctrl_freq * 100,
    donor_short = str_remove(Donor_Group, "_(pp65|IE1)$"),
    point_label = paste0("Clone ", str_extract(Global_Clone_ID, "\\d+$"))
  )


# ── 7. PLOT BUILDER ──────────────────────────────────────────────────────────

safe_name <- function(x) {
  x %>% str_replace_all("[^A-Za-z0-9]+", "_") %>% str_remove("^_|_$")
}

# Padded range, used for axes with no fixed limits (i.e. raw scores)
padded_range <- function(values, pad_frac = 0.08) {
  v <- values[is.finite(values)]
  if (length(v) == 0) return(c(-1, 1))
  rng <- range(v)
  if (diff(rng) == 0) rng <- rng + c(-1, 1)
  pad <- diff(rng) * pad_frac
  c(rng[1] - pad, rng[2] + pad)
}

make_corr_plot <- function(plot_df, x_lab, y_lab, title,
                           x_lim, y_lim, x_breaks, y_breaks) {

  fit  <- lm(y ~ x, data = plot_df)
  r2   <- summary(fit)$r.squared
  pval <- tryCatch(summary(fit)$coefficients["x", "Pr(>|t|)"],
                   error = function(e) NA_real_)

  pval_text <- if (is.na(pval)) "p = NA" else
    if (pval < 0.001) "p < 0.001" else sprintf("p = %.3f", pval)
  annot <- sprintf("R² = %.3f\n%s", r2, pval_text)

  ggplot(plot_df, aes(x = x, y = y)) +
    geom_smooth(method = "lm", formula = y ~ x, se = se_ribbon,
                colour = "steelblue", fill = "grey85", linewidth = 0.8) +
    geom_point(aes(colour = donor_short), size = 2.5) +
    geom_text_repel(
      aes(label = point_label),
      size          = clone_label_size,
      max.overlaps  = 20,
      segment.size  = 0.3,
      segment.color = "grey60"
    ) +
    annotate("text",
             x     = x_lim[1] + 0.03 * diff(x_lim),
             y     = y_lim[2] - 0.03 * diff(y_lim),
             hjust = 0, vjust = 1,
             label = annot,
             size  = annot_size,
             colour = "grey20") +
    scale_colour_manual(values = donor_colours, name = NULL, drop = TRUE) +
    guides(colour = guide_legend(override.aes = list(size = 3.5))) +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(breaks = y_breaks) +
    # coord_cartesian rather than scale limits so no point is dropped from the fit
    coord_cartesian(xlim = x_lim, ylim = y_lim) +
    labs(title = title, x = x_lab, y = y_lab) +
    theme_classic(base_size = 11 * font_scale) +
    theme(
      plot.title      = element_text(size = title_size, face = "bold"),
      axis.title      = element_text(size = axis_title_size),
      axis.text       = element_text(size = axis_text_size),
      legend.text     = element_text(size = axis_text_size),
      legend.position = "right"
    )
}


# ── 8. PLOT SPECIFICATIONS ───────────────────────────────────────────────────
#
#   x = TCR score    vs  other assay scores
#   x = TCR score    vs  control-well frequency
#   x = TCR score    vs  % cells in each metacluster
#   x = CD107a score vs  % cells in each metacluster   [new]

score_suffix     <- if (use_percentiles) "_score_percentile" else "_score"
score_lab_suffix <- if (use_percentiles) " score percentile" else " score"
score_limits     <- if (use_percentiles) percentile_limits else NULL

build_specs <- function(available_cols) {

  specs <- list()

  add <- function(specs, x_col, x_short, x_lab, x_lim,
                  y_col, y_short, y_lab, y_lim, y_tag) {
    if (!all(c(x_col, y_col) %in% available_cols)) return(specs)
    append(specs, list(list(
      x_col = x_col, x_short = x_short, x_lab = x_lab, x_lim = x_lim,
      y_col = y_col, y_short = y_short, y_lab = y_lab, y_lim = y_lim,
      y_tag = y_tag
    )))
  }

  # (a) TCR vs other assay scores
  for (assay in c("Proliferation", "CD107a", "IFNg")) {
    specs <- add(specs,
                 paste0("TCR", score_suffix), "TCR",
                 paste0("TCR", score_lab_suffix), score_limits,
                 paste0(assay, score_suffix), assay,
                 paste0(assay, score_lab_suffix), score_limits,
                 safe_name(paste0(assay, score_suffix)))
  }

  # (b) TCR vs clone dominance within the control well (%)
  specs <- add(specs,
               paste0("TCR", score_suffix), "TCR",
               paste0("TCR", score_lab_suffix), score_limits,
               "ctrl_pct", "% Dominance",
               "% of clone within donor/specificity",
               pct_limits, "Control_frequency")

  # (c) and (d) TCR and CD107a vs metacluster composition
  for (x_assay in c("TCR", "CD107a")) {
    for (pct_col in pct_cols) {
      cluster_name <- str_remove(pct_col, "^pct_")
      specs <- add(specs,
                   paste0(x_assay, score_suffix), x_assay,
                   paste0(x_assay, score_lab_suffix), score_limits,
                   pct_col, paste0("% ", cluster_name),
                   paste0("% cells in ", cluster_name), pct_limits,
                   paste0("pct_", safe_name(cluster_name)))
    }
  }

  specs
}

plot_specs <- build_specs(names(clone_table))


# ── 9. GENERATE PLOTS ────────────────────────────────────────────────────────

stats_rows <- list()

for (grp_name in names(combined_groups)) {

  grp       <- combined_groups[[grp_name]]
  grp_df    <- clone_table %>% filter(Donor_Group %in% grp$donors)
  all_short <- unique(str_remove(grp$donors, "_(pp65|IE1)$"))

  if (nrow(grp_df) == 0) {
    warning(sprintf("No clones for combined group %s - skipping", grp_name))
    next
  }

  # One shared range per score column across every panel of this specificity
  shared_lims <- list()
  if (shared_score_limits) {
    for (cc in intersect(paste0(c("TCR", "Proliferation", "CD107a", "IFNg"),
                                score_suffix), names(grp_df))) {
      shared_lims[[cc]] <- padded_range(grp_df[[cc]])
    }
  }

  resolve_lim <- function(spec_lim, col, values) {
    if (!is.null(spec_lim)) return(spec_lim)
    if (!is.null(shared_lims[[col]])) return(shared_lims[[col]])
    padded_range(values)
  }

  for (spec in plot_specs) {

    plot_df <- grp_df %>%
      select(Global_Clone_ID, Donor_Group, donor_short, point_label,
             x = all_of(spec$x_col), y = all_of(spec$y_col)) %>%
      filter(!is.na(x), !is.na(y))

    if (nrow(plot_df) < 3) next

    x_lim <- resolve_lim(spec$x_lim, spec$x_col, plot_df$x)
    y_lim <- resolve_lim(spec$y_lim, spec$y_col, plot_df$y)

    # Fixed axes get evenly spaced breaks; free (raw score) axes use defaults
    x_breaks <- if (!is.null(spec$x_lim))
      seq(spec$x_lim[1], spec$x_lim[2], length.out = 5) else waiver()
    y_breaks <- if (!is.null(spec$y_lim))
      seq(spec$y_lim[1], spec$y_lim[2], length.out = 5) else waiver()

    # Which donors actually contributed data for this contrast
    donors_present <- plot_df %>% pull(donor_short) %>% unique() %>% sort()

    donor_tag <- if (length(donors_present) < length(all_short)) {
      paste0(" (", paste(donors_present, collapse = ", "), ")")
    } else ""

    title <- paste0(grp$label, donor_tag, "  |  ",
                    spec$x_short, " vs ", spec$y_short)

    p <- make_corr_plot(
      plot_df  = plot_df,
      x_lab    = spec$x_lab,
      y_lab    = spec$y_lab,
      title    = title,
      x_lim    = x_lim,
      y_lim    = y_lim,
      x_breaks = x_breaks,
      y_breaks = y_breaks
    )

    filename <- paste0(grp$label, "_", safe_name(spec$x_short),
                       "_vs_", spec$y_tag, ".png")

    ggsave(
      filename = file.path(output_dir, filename),
      plot     = p,
      width    = plot_width,
      height   = plot_height,
      dpi      = plot_dpi
    )

    # Record fit statistics
    fit  <- lm(y ~ x, data = plot_df)
    pval <- tryCatch(summary(fit)$coefficients["x", "Pr(>|t|)"],
                     error = function(e) NA_real_)

    stats_rows[[length(stats_rows) + 1]] <- tibble(
      Specificity = grp$label,
      Donors      = paste(donors_present, collapse = ", "),
      X           = spec$x_short,
      Y           = spec$y_short,
      n_clones    = nrow(plot_df),
      R2          = summary(fit)$r.squared,
      p_value     = pval,
      File        = filename
    )

    message(sprintf("  Saved: %s (%d clones)", filename, nrow(plot_df)))
  }
}

if (length(stats_rows) > 0) {
  correlation_stats <- bind_rows(stats_rows)
  write.csv(correlation_stats,
            file.path(output_dir, "correlation_stats.csv"),
            row.names = FALSE)
  message(sprintf("\nWrote correlation_stats.csv (%d plots)",
                  nrow(correlation_stats)))
}

message("\nDone. Plots written to: ", output_dir)



###############################################################################
# major_clone_barplots_for_figure.R
#
# Two figure panels per donor, for the pp65-specific major clones of D1-D4:
#
#   1. "Frequency of major clones' TCRs from each sort (donor N, pp65)"
#      Grouped bars: Control | TCR-low | TCR-high, one bar per major clone
#      within each group, coloured by the global clone palette.
#
#   2. "Relative TCR avidity score"
#      The same bar style, one group, showing each clone's raw TCR score
#      (CLR(TCR-high) - CLR(TCR-low)).
#
# Reads the exported UMI count table; no TRUST4 re-parsing.
###############################################################################

# ── 0. PACKAGES ───────────────────────────────────────────────────────────────
library(tidyverse)
library(compositions)


# ── 1. CONFIGURATION ──────────────────────────────────────────────────────────

results_dir <- "Bulk_TCR_seq/Final_analysis_results"
umi_file    <- file.path(results_dir, "TRB_clone_UMI_counts.csv")

output_dir <- file.path(results_dir, "barplots_for_figure")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

pseudocount <- 1

# ── Donors to plot. Wells are named by the label that appears on the x axis,
# in the left-to-right order they should be drawn.
donor_panels <- list(
  list(
    prefix = "D1_pp65", label = "donor 1, pp65",
    wells  = c("Control"  = "Full-11-Control",
               "TCR-low"  = "Full-7-TCR-low",
               "TCR-high" = "Full-9-TCR-high")
  ),
  list(
    prefix = "D2_pp65", label = "donor 2, pp65",
    wells  = c("Control"  = "1-D2-p-Control",
               "TCR-low"  = "3-D2-p-TCR-low",
               "TCR-high" = "4-D2-p-TCR-high")
  ),
  list(
    prefix = "D3_pp65", label = "donor 3, pp65",
    wells  = c("Control"  = "5-D3-p-Control",
               "TCR-low"  = "7-D3-p-TCR-low",
               "TCR-high" = "8-D3-p-TCR-high")
  ),
  list(
    prefix = "D4_pp65", label = "donor 4, pp65",
    wells  = c("Control"  = "9-D4-p-Control",
               "TCR-low"  = "11-D4-p-TCR-low",
               "TCR-high" = "12-D4-p-TCR-high")
  )
)

# ── Which clones count as "major".
# Any prefix not listed here uses every clone carrying that prefix, ordered by
# clone number. Override individual donors like:
#   major_clones_override <- list("D3_pp65" = paste0("D3_pp65-clone-", 1:6))
major_clones_override <- list()

# ── Denominator for the frequencies in panel 1
#   "major_clones"      : each clone / sum of the major clones in that well
#   "all_donor_clones"  : each clone / sum of all this donor's clones in that
#                         well (matches the donor-grouped frequency exports)
# Identical if the major clones are the only clones for that donor.
frequency_denominator <- "major_clones"

# ── Share y-axis ranges across donors so the panels can be compared directly.
# FALSE gives each donor its own range.
shared_y_across_donors <- FALSE

# ─────────────────────────────────────────────────────────────────────────────
# GLOBAL CLONE PALETTE (same one used for the pie charts)
#
# Keyed by legend label ("Clone 1"). Full Global_Clone_IDs and bare numbers
# also resolve, as does an unnamed vector applied in clone order. The
# "Other (<1%)" entry is ignored here since only named clones are plotted.
# ─────────────────────────────────────────────────────────────────────────────

clone_colors_pp65 <- c(
  "Clone 1"     = "#1b9e77",
  "Clone 2"     = "#d95f02",
  "Clone 3"     = "#7570b3",
  "Clone 4"     = "#e7298a",
  "Clone 5"     = "#66a61e",
  "Clone 6"     = "#e6ab02",
  "Clone 7"     = "#a6761d",
  "Clone 8"     = "#666666",
  "Clone 9"     = "#1f78b4",
  "Clone 10"    = "#b2df8a",
  "Clone 11"    = "#33a02c",
  "Clone 12"    = "#fb9a99",
  "Clone 13"    = "#e31a1c",
  "Clone 14"    = "#fdbf6f",
  "Clone 15"    = "#ff7f00",
  "Clone 16"    = "#cab2d6",
  "Clone 17"    = "#6a3d9a",
  "Clone 18"    = "#ffff99",
  "Clone 19"    = "#b15928",
  "Clone 20"    = "#8dd3c7",
  "Other (<1%)" = "grey70"
)

clone_colors_ie1 <- c(
  "Clone 1"     = "#1b9e77",
  "Clone 2"     = "#d95f02",
  "Clone 3"     = "#7570b3",
  "Clone 4"     = "#e7298a",
  "Clone 5"     = "#66a61e",
  "Clone 6"     = "#e6ab02",
  "Clone 7"     = "#a6761d",
  "Other (<1%)" = "grey70"
)

palette_for <- function(prefix) {
  if (str_detect(prefix, "IE1")) clone_colors_ie1 else clone_colors_pp65
}

# ── Axis / title text.
# A "%s" in either title template is filled with the donor label above;
# leave it out for a fixed title.
panel1_title_template <- "Frequency of major clones' TCRs from each sort (%s)"
panel1_ylab           <- "% of clonal TCRs"

panel2_title_template <- "Relative TCR avidity score (%s)"
panel2_ylab           <- "TCR score"
score_group_label     <- "TCR-high vs TCR-low"   # single x-axis group in panel 2

# ── X-axis tick labels for panel 1, keyed by the well names in donor_panels.
# Anything not listed here is drawn using its well name as-is.
sort_display_labels <- c(
  "Control"  = "Control sort",
  "TCR-low"  = "TCR-low sort",
  "TCR-high" = "TCR-high sort"
)

# ── Bar geometry
group_width <- 0.80   # width of a whole group; smaller = bigger gap BETWEEN groups
bar_padding <- 0.05   # sliver between bars WITHIN a group

# ── Dashed separator drawn in the gap between groups (multi-group panels only)
group_divider    <- TRUE
divider_colour   <- "grey60"
divider_linetype <- "dashed"
divider_width    <- 0.4

# ── Panel 1 x axis: every bar gets its own tick and rotated "Clone N" label,
# with the sort name below them. The drop of the sort names is computed from the
# clone label size and angle, so it follows any change made here.
x_expand           <- 0.20   # padding at each end of the x axis
clone_label_size   <- 7      # pt, the rotated "Clone N" labels
clone_label_angle  <- 45
group_label_gap_pt <- 10     # clear space between the two rows of labels
char_width_frac    <- 0.62   # average glyph width as a fraction of font size
text_height_frac   <- 0.72   # glyph height as a fraction of font size
axis_offset_pt     <- 5      # tick length + axis text margin below the panel
panel1_legend      <- FALSE  # the bars are labelled, so the legend is redundant

# ── Minimum drawn bar height, so a clone that is absent from a sort — or too
# rare to render — still marks its slot. A fraction of the y range; the value
# is a floor on the drawn height only, not on the exported numbers.
zero_stub      <- TRUE
zero_stub_frac <- 0.010

# ── Font sizes: 1.2x the originals, matching the correlation plots
font_scale      <- 1.2
title_size      <- 10  * font_scale
axis_title_size <- 9   * font_scale
axis_text_size  <- 8.8 * font_scale

# ── Canvas (inches)
plot_height_base <- 4.5
y_axis_stretch   <- 1.2
non_panel_height <- 1.1
height_scale     <- 0.50   # short, wide panels
width_scale      <- 1.125  # 1.5 x 0.75
plot_height <- (non_panel_height +
                (plot_height_base - non_panel_height) * y_axis_stretch) * height_scale

bar_width_in    <- 0.18   # horizontal space per bar, before width_scale
axis_width_in   <- 0.80   # y-axis text and title
legend_width_in <- 1.30   # clone legend
plot_dpi        <- 300


# ── 2. IMPORT ─────────────────────────────────────────────────────────────────
# read_csv (not read.csv) so well names like "1-D2-p-Control" are not mangled.

umi_summary <- read_csv(umi_file, show_col_types = FALSE)

message(sprintf("UMI table loaded: %d clones x %d wells",
                nrow(umi_summary), ncol(umi_summary) - 1))


# ── 3. HELPERS ────────────────────────────────────────────────────────────────

clone_label_of <- function(clone_ids) {
  paste0("Clone ", str_extract(clone_ids, "\\d+$"))
}

display_sort <- function(x) {
  out <- unname(sort_display_labels[x])
  ifelse(is.na(out), x, out)
}

clones_with_prefix <- function(prefix) {
  umi_summary %>%
    filter(str_starts(Global_Clone_ID, prefix)) %>%
    mutate(clone_num = as.integer(str_extract(Global_Clone_ID, "\\d+$"))) %>%
    arrange(clone_num) %>%
    pull(Global_Clone_ID)
}

resolve_clone_colours <- function(pal, clone_ids, labels) {

  if (is.null(pal)) stop("No palette supplied for these clones")

  # Unnamed palette: applied in clone order
  if (is.null(names(pal))) {
    if (length(pal) < length(clone_ids)) {
      stop(sprintf("Palette has %d colours but %d clones are being plotted",
                   length(pal), length(clone_ids)))
    }
    return(setNames(unname(pal)[seq_along(clone_ids)], labels))
  }

  nums <- str_extract(clone_ids, "\\d+$")
  out  <- rep(NA_character_, length(clone_ids))

  for (i in seq_along(clone_ids)) {
    keys <- c(clone_ids[i], labels[i], nums[i])
    hit  <- keys[keys %in% names(pal)]
    if (length(hit) > 0) out[i] <- unname(pal[[hit[1]]])
  }

  if (any(is.na(out))) {
    stop("No colour found in the palette for: ",
         paste(labels[is.na(out)], collapse = ", "))
  }

  setNames(out, labels)
}

# Width that fits the bars, and never narrower than the title needs
panel_width <- function(n_groups, n_clones, title, legend = TRUE) {
  legend_in <- if (legend) legend_width_in else 0.2
  bars  <- axis_width_in + legend_in +
           n_groups * n_clones * bar_width_in + n_groups * 0.15
  title <- nchar(title) * title_size * 0.55 / 72 + 0.4
  max(4.5, bars, title) * width_scale
}

# How far the sort names must sit below the panel to clear the rotated clone
# labels, plus the canvas space that needs reserving for both rows.
group_label_geometry <- function(labels) {
  rotated <- max(nchar(labels)) * clone_label_size * char_width_frac *
             sin(clone_label_angle * pi / 180)
  # Distance from the panel bottom to the top of the sort name: the clone labels
  # already start axis_offset_pt below the panel and run rotated pt further.
  drop   <- axis_offset_pt + rotated + group_label_gap_pt
  margin <- drop + 1.4 * axis_text_size
  list(
    # vjust is measured in glyph heights, not font sizes, hence text_height_frac
    vjust     = 1 + drop / (axis_text_size * text_height_frac),
    margin_pt = margin,
    extra_in  = (rotated + margin) / 72
  )
}

bar_theme <- function() {
  theme_classic(base_size = 11 * font_scale) +
    theme(
      plot.title      = element_text(size = title_size, face = "bold", hjust = 0.5),
      plot.title.position = "panel",
      axis.title      = element_text(size = axis_title_size),
      axis.text       = element_text(size = axis_text_size),
      legend.text     = element_text(size = axis_text_size),
      legend.position = "right"
    )
}

bar_layer <- function() {
  geom_col(
    position = position_dodge2(padding = bar_padding, preserve = "single"),
    width    = group_width,
    colour   = NA
  )
}

make_title <- function(template, label) {
  if (str_detect(template, "%s")) sprintf(template, label) else template
}

# Centre of every bar slot, on an x axis where group g sits at g. Used both to
# place the bars and as the axis breaks, so each bar gets its own tick.
bar_tick_positions <- function(n_groups, n_clones) {
  slot <- group_width / n_clones
  as.vector(vapply(
    seq_len(n_groups),
    function(g) g - group_width / 2 + (seq_len(n_clones) - 0.5) * slot,
    numeric(n_clones)
  ))
}


# ── 4. PASS ONE: BUILD THE DATA FOR EVERY DONOR ──────────────────────────────

panel_data <- list()

for (panel in donor_panels) {

  prefix <- panel$prefix
  wells  <- panel$wells

  missing_wells <- setdiff(unname(wells), colnames(umi_summary))
  if (length(missing_wells) > 0) {
    warning(sprintf("%s: wells not found (%s) - skipping",
                    prefix, paste(missing_wells, collapse = ", ")))
    next
  }

  donor_clones <- clones_with_prefix(prefix)
  if (length(donor_clones) == 0) {
    warning(sprintf("%s: no clones found - skipping", prefix))
    next
  }

  major <- if (!is.null(major_clones_override[[prefix]])) {
    major_clones_override[[prefix]]
  } else {
    donor_clones
  }

  missing_clones <- setdiff(major, umi_summary$Global_Clone_ID)
  if (length(missing_clones) > 0) {
    stop(sprintf("%s: clones not in the UMI table: %s",
                 prefix, paste(missing_clones, collapse = ", ")))
  }

  labels <- clone_label_of(major)

  denom_clones <- switch(
    frequency_denominator,
    "major_clones"     = major,
    "all_donor_clones" = donor_clones,
    stop("frequency_denominator must be 'major_clones' or 'all_donor_clones'")
  )

  # ── Panel 1: per-sort frequencies
  freq_df <- umi_summary %>%
    filter(Global_Clone_ID %in% denom_clones) %>%
    select(Global_Clone_ID, all_of(unname(wells))) %>%
    mutate(across(-Global_Clone_ID, ~ . / sum(.) * 100)) %>%
    filter(Global_Clone_ID %in% major) %>%
    pivot_longer(-Global_Clone_ID, names_to = "well", values_to = "Percent") %>%
    mutate(
      Sort  = factor(display_sort(names(wells)[match(well, wells)]),
                     levels = display_sort(names(wells))),
      Clone = factor(clone_label_of(Global_Clone_ID), levels = labels)
    )

  # ── Panel 2: raw TCR score.
  # CLR per well across this donor's clones, exactly as in the
  # Master_Assay_Scores calculation, then TCR-high minus TCR-low.
  score_pair <- umi_summary %>%
    filter(Global_Clone_ID %in% donor_clones) %>%
    select(Global_Clone_ID, all_of(unname(wells[c("TCR-high", "TCR-low")]))) %>%
    mutate(across(-Global_Clone_ID, ~ {
      x <- . + pseudocount
      as.numeric(clr(acomp(x / sum(x))))
    }))

  score_df <- tibble(
    Global_Clone_ID = score_pair$Global_Clone_ID,
    TCR_score       = score_pair[[unname(wells["TCR-high"])]] -
                      score_pair[[unname(wells["TCR-low"])]]
  ) %>%
    filter(Global_Clone_ID %in% major) %>%
    mutate(
      Group = factor(score_group_label),
      Clone = factor(clone_label_of(Global_Clone_ID), levels = labels)
    )

  panel_data[[prefix]] <- list(
    prefix      = prefix,
    label       = panel$label,
    n_groups    = length(wells),
    n_clones    = length(major),
    fill_values = resolve_clone_colours(palette_for(prefix), major, labels),
    freq_df     = freq_df,
    score_df    = score_df
  )

  message(sprintf("  %s: %d clones (%s)",
                  prefix, length(major), paste(labels, collapse = ", ")))
}

if (length(panel_data) == 0) stop("No donor panels could be built")


# ── 5. SHARED Y RANGES ───────────────────────────────────────────────────────

if (shared_y_across_donors) {
  max_freq  <- max(map_dbl(panel_data, ~ max(.x$freq_df$Percent, na.rm = TRUE)))
  max_score <- max(map_dbl(panel_data, ~ max(abs(.x$score_df$TCR_score), na.rm = TRUE)))
  freq_ylim  <- c(0, max_freq * 1.05)
  score_ylim <- c(-max_score * 1.10, max_score * 1.10)   # symmetric about zero
} else {
  freq_ylim  <- NULL
  score_ylim <- NULL
}


# ── 6. PASS TWO: DRAW AND SAVE ───────────────────────────────────────────────

for (pd in panel_data) {

  # ── Panel 1: frequency by sort
  title1 <- make_title(panel1_title_template, pd$label)

  # Explicit numeric x positions so every bar can carry its own tick and label:
  # group g sits at g, and its bars fill group_width centred on it.
  n_clones <- pd$n_clones
  slot     <- group_width / n_clones

  y_upper <- if (!is.null(freq_ylim)) freq_ylim[2]
             else max(pd$freq_df$Percent, na.rm = TRUE)

  # Floor every bar at a thin stub, so a clone that is absent from a sort — or
  # present at a level too small to render as a visible bar — still marks its
  # slot. Only the drawn height is floored; the exported values are untouched.
  stub_h <- if (zero_stub) y_upper * zero_stub_frac else 0

  freq_plot_df <- pd$freq_df %>%
    mutate(
      x_pos        = as.integer(Sort) - group_width / 2 +
                     (as.integer(Clone) - 0.5) * slot,
      Percent_plot = pmax(Percent, stub_h)
    )

  p1 <- ggplot(freq_plot_df, aes(x = x_pos, y = Percent_plot, fill = Clone))

  if (group_divider && pd$n_groups > 1) {
    p1 <- p1 + geom_vline(
      xintercept = seq_len(pd$n_groups - 1) + 0.5,
      colour     = divider_colour,
      linetype   = divider_linetype,
      linewidth  = divider_width
    )
  }

  glab <- group_label_geometry(levels(pd$freq_df$Clone))

  p1 <- p1 +
    geom_col(width = slot * (1 - bar_padding), colour = NA) +
    scale_fill_manual(values = pd$fill_values, name = NULL, drop = FALSE) +
    scale_x_continuous(
      breaks = bar_tick_positions(pd$n_groups, n_clones),
      labels = rep(levels(pd$freq_df$Clone), times = pd$n_groups),
      expand = expansion(add = x_expand)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    # Sort names, below the rotated clone labels. vjust is in units of the
    # label's own height, so it tracks the font size rather than the data range.
    annotate("text",
             x     = seq_len(pd$n_groups),
             y     = -Inf,
             label = levels(pd$freq_df$Sort),
             vjust = glab$vjust,
             size  = axis_text_size / .pt) +
    labs(title = title1, x = NULL, y = panel1_ylab) +
    bar_theme() +
    theme(
      axis.text.x     = element_text(size  = clone_label_size,
                                     angle = clone_label_angle,
                                     hjust = 1, vjust = 1),
      legend.position = if (panel1_legend) "right" else "none",
      plot.margin     = margin(t = 5.5, r = 5.5, b = glab$margin_pt, l = 5.5)
    ) +
    # clip = "off" lets the sort names draw below the panel
    coord_cartesian(ylim = c(0, y_upper), clip = "off")

  ggsave(
    filename = file.path(output_dir,
                         paste0(pd$prefix, "_major_clone_frequency_by_sort.png")),
    plot   = p1,
    width  = panel_width(pd$n_groups, pd$n_clones, title1, legend = panel1_legend),
    height = plot_height + glab$extra_in,
    dpi    = plot_dpi
  )

  # ── Panel 2: raw TCR score
  title2 <- make_title(panel2_title_template, pd$label)

  p2 <- ggplot(pd$score_df, aes(x = Group, y = TCR_score, fill = Clone)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
    bar_layer() +
    scale_fill_manual(values = pd$fill_values, name = NULL, drop = FALSE) +
    labs(title = title2, x = NULL, y = panel2_ylab) +
    bar_theme() +
    # centred over the whole plot, legend included, rather than over the panel
    theme(plot.title.position = "plot")

  if (!is.null(score_ylim)) p2 <- p2 + coord_cartesian(ylim = score_ylim)

  ggsave(
    filename = file.path(output_dir,
                         paste0(pd$prefix, "_major_clone_TCR_score.png")),
    plot   = p2,
    width  = panel_width(1, pd$n_clones, title2),
    height = plot_height,
    dpi    = plot_dpi
  )

  # ── Plotted values, for the figure legend
  pd$freq_df %>%
    select(Global_Clone_ID, Clone, Sort, Percent) %>%
    arrange(Sort, Clone) %>%
    write.csv(file.path(output_dir,
                        paste0(pd$prefix, "_major_clone_frequency_by_sort.csv")),
              row.names = FALSE)

  pd$score_df %>%
    select(Global_Clone_ID, Clone, TCR_score) %>%
    arrange(Clone) %>%
    write.csv(file.path(output_dir,
                        paste0(pd$prefix, "_major_clone_TCR_score.csv")),
              row.names = FALSE)

  message(sprintf("  Saved both panels for %s", pd$prefix))
}

message("\nDone. Plots written to: ", output_dir)




###############################################################################
# correlation_plots_for_figure_v2.R   (RE-EXPORT)
#
# Same plots as correlation_plots_for_figure.R, regenerated with:
#   * no per-point "Clone N" labels
#   * larger points (2.5 -> 3.5)
#   * all text 1.5x the current figure versions: titles, axis titles, tick
#     labels, legend text and the R2 / p annotation
#   * axis titles "Relative TCR avidity score" / "Relative proliferation score"
#   * optional two-line titles (specificity/donors, then the comparison)
#   * legend pulled in closer to the plotting panel
# Scores, axis limits, colours and file names are unchanged.
# Output: re-exports_for_figures/correlation_plots_for_figure/
###############################################################################

# ── 0. PACKAGES ───────────────────────────────────────────────────────────────
library(tidyverse)
library(compositions)


# ── 1. CONFIGURATION ──────────────────────────────────────────────────────────

results_dir     <- "Bulk_TCR_seq/Final_analysis_results"
umi_file        <- file.path(results_dir, "TRB_clone_UMI_counts.csv")
metacluster_csv <- "Bulk_TCR_seq/Global_Clone_metacluster_composition.csv"

reexport_dir <- "re-exports_for_figures"
output_dir   <- file.path(reexport_dir, "correlation_plots_for_figure")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

pseudocount <- 1

# ── Assay axis units: raw CLR-difference scores, or within-donor percentiles
use_percentiles <- FALSE

# ── Fixed axis limits
pct_limits        <- c(0, 100)    # metacluster percentages and clone dominance
percentile_limits <- c(0, 100)    # assay score percentiles (if used)

# Raw scores have no natural bounds. TRUE = one shared range per score across
# all panels of a specificity, so panels can be compared side by side.
shared_score_limits <- TRUE

# ── Axis titles for the assay scores (" percentile" is appended when
# use_percentiles = TRUE)
score_axis_labels <- c(
  TCR           = "Relative TCR avidity score",
  Proliferation = "Relative proliferation score",
  CD107a        = "CD107a score",
  IFNg          = "IFNg score"
)

score_axis_label <- function(assay) {
  lab <- score_axis_labels[[assay]]
  if (use_percentiles) paste0(lab, " percentile") else lab
}

# ── Text sizes. text_scale multiplies the CURRENT figure versions, which were
# already 1.2x the sizes in the original bulk_tcrseq_analysis.R.
font_scale_prev <- 1.2
text_scale      <- 1.5
fs              <- font_scale_prev * text_scale   # overall factor vs the original

base_size        <- 11 * fs     # pt, passed to theme_classic()
half_line        <- base_size / 2
title_size       <- 10  * fs    # pt
axis_title_size  <- 9   * fs    # pt
axis_text_size   <- 8.8 * fs    # pt (theme_classic(11) default tick text ~8.8)
legend_text_size <- 8.8 * fs    # pt
annot_size       <- 3.2 * fs    # R2 / p annotation (ggplot text-size units)

# ── Points
point_size        <- 3.5        # was 2.5
point_alpha       <- 1          # lower (e.g. 0.85) if overlapping points hide each other
legend_point_size <- point_size

# ── Titles: TRUE = "pp65 (D1, D4)" on line 1, "TCR vs % Tem" on line 2.
# At the larger font a single-line title can run off the canvas.
two_line_titles <- FALSE

# ── Legend placement.
# ggplot's default gap between the panel and the legend (legend.box.spacing
# plus the legend's own left margin) scales with base_size, so at the enlarged
# text it was ~30 pt. legend_gap_pt replaces both with a single gap.
legend_gap_pt <- 6              # raise to push the legend back out to the right

# FALSE = canvas size unchanged; the freed-up space goes to the panel, which
#         gets slightly wider. Nothing else on the plot moves.
# TRUE  = panel width unchanged; the canvas is trimmed by the freed-up space
#         instead. Check the longest single-line titles if you use this.
trim_canvas_for_legend <- FALSE

# ── Canvas (inches).
# keep_panel_size = TRUE keeps the plotting panel the same physical size as the
# current figure versions and grows the canvas around it to fit the larger
# text, so text and points get bigger RELATIVE to the panel.
# FALSE keeps the current canvas size instead (the panel shrinks to make room).
keep_panel_size <- TRUE

base_width       <- 6                                  # current canvas width
base_height      <- 1.1 + (4.5 - 1.1) * 1.2            # current canvas height (5.18)
non_panel_width  <- 1.6   # approx: y-axis text + title, legend, margins (current text)
non_panel_height <- 1.1   # approx: title, x-axis text + title, margins (current text)
title_line_in    <- title_size * 1.2 / 72              # one extra title line

if (keep_panel_size) {
  plot_width  <- (base_width  - non_panel_width)  + non_panel_width  * text_scale
  plot_height <- (base_height - non_panel_height) + non_panel_height * text_scale +
                 (if (two_line_titles) title_line_in else 0)
} else {
  plot_width  <- base_width
  plot_height <- base_height
}

# Space freed up by the tighter legend gap (default = base_size + half_line)
default_legend_gap_pt <- base_size + half_line
legend_saving_in      <- max(0, default_legend_gap_pt - legend_gap_pt) / 72
if (trim_canvas_for_legend) plot_width <- plot_width - legend_saving_in

plot_dpi <- 300

message(sprintf("Canvas: %.2f x %.2f in  (legend gap %.0f pt, %.2f in freed)",
                plot_width, plot_height, legend_gap_pt, legend_saving_in))

# ── Point colours, fixed per donor so a given donor is the same colour in
# every panel regardless of which donors have data for that contrast
donor_colours <- c(
  D1 = "#4E9AC7",
  D2 = "#E07B3A",
  D3 = "#57A65A",
  D4 = "#9B4FA6",
  D5 = "#D9B23A"
)

# ── Options
se_ribbon <- FALSE   # TRUE = draw the 95% CI ribbon around the lm line


# ── 2. DONOR GROUPS AND CONTRASTS ─────────────────────────────────────────────

donor_groups <- list(
  D1_pp65 = list(
    wells = c("Full-11-Control", "Full-10-Prolif",
              "Full-7-TCR-low",  "Full-9-TCR-high",
              "Full-2-107a-pos", "Full-3-107a-neg",
              "Full-5-IFNg-pos", "Full-6-IFNg-neg"),
    clone_prefix = "D1_pp65"
  ),
  D2_pp65 = list(
    wells = c("1-D2-p-Control", "2-D2-p-Prolif",
              "3-D2-p-TCR-low", "4-D2-p-TCR-high"),
    clone_prefix = "D2_pp65"
  ),
  D3_pp65 = list(
    wells = c("5-D3-p-Control", "6-D3-p-Prolif",
              "7-D3-p-TCR-low", "8-D3-p-TCR-high"),
    clone_prefix = "D3_pp65"
  ),
  D4_pp65 = list(
    wells = c("9-D4-p-Control",   "10-D4-p-Prolif",
              "11-D4-p-TCR-low",  "12-D4-p-TCR-high",
              "13-D4-p-107a-neg", "14-D4-p-107a-pos"),
    clone_prefix = "D4_pp65"
  ),
  D4_IE1 = list(
    wells = c("15-D4-I-Control",  "16-D4-I-Prolif",
              "17-D4-I-TCR-low",  "18-D4-I-TCR-high",
              "19-D4-I-107a-neg", "20-D4-I-107a-pos"),
    clone_prefix = "D4_IE1"
  ),
  D5_IE1 = list(
    wells = c("21-D5-I-Control", "22-D5-I-Prolif",
              "23-D5-I-TCR-low", "24-D5-I-TCR-high"),
    clone_prefix = "D5_IE1"
  )
)

group_contrasts <- list(
  D1_pp65 = list(
    list(pos = "Prolif",   neg = "Control",   name = "Proliferation"),
    list(pos = "TCR-high", neg = "TCR-low",   name = "TCR"),
    list(pos = "107a-pos", neg = "107a-neg",  name = "CD107a"),
    list(pos = "IFNg-pos", neg = "IFNg-neg",  name = "IFNg")
  ),
  D2_pp65 = list(
    list(pos = "Prolif",   neg = "Control",   name = "Proliferation"),
    list(pos = "TCR-high", neg = "TCR-low",   name = "TCR")
  ),
  D3_pp65 = list(
    list(pos = "Prolif",   neg = "Control",   name = "Proliferation"),
    list(pos = "TCR-high", neg = "TCR-low",   name = "TCR")
  ),
  D4_pp65 = list(
    list(pos = "Prolif",   neg = "Control",   name = "Proliferation"),
    list(pos = "TCR-high", neg = "TCR-low",   name = "TCR"),
    list(pos = "107a-pos", neg = "107a-neg",  name = "CD107a")
  ),
  D4_IE1 = list(
    list(pos = "Prolif",   neg = "Control",   name = "Proliferation"),
    list(pos = "TCR-high", neg = "TCR-low",   name = "TCR"),
    list(pos = "107a-pos", neg = "107a-neg",  name = "CD107a")
  ),
  D5_IE1 = list(
    list(pos = "Prolif",   neg = "Control",   name = "Proliferation"),
    list(pos = "TCR-high", neg = "TCR-low",   name = "TCR")
  )
)

# ── Pooled plotting groups: all clones of a given specificity together
combined_groups <- list(
  pp65 = list(label = "pp65", donors = c("D1_pp65", "D2_pp65", "D3_pp65", "D4_pp65")),
  IE1  = list(label = "IE1",  donors = c("D4_IE1", "D5_IE1"))
)

metacluster_levels <- c("Tcm", "Tem", "T-CMC1", "NK-like", "ISG+")


# ── 3. IMPORT EXPORTED TABLES ─────────────────────────────────────────────────
# read_csv (not read.csv) so well names like "1-D2-p-Control" are not mangled.

umi_summary <- read_csv(umi_file, show_col_types = FALSE)

message(sprintf("UMI table loaded: %d clones x %d wells",
                nrow(umi_summary), ncol(umi_summary) - 1))

if (file.exists(metacluster_csv)) {
  metacluster_comp <- read_csv(metacluster_csv, show_col_types = FALSE)
  pct_cols <- intersect(paste0("pct_", metacluster_levels), names(metacluster_comp))
  metacluster_comp <- metacluster_comp %>% select(Global_Clone_ID, all_of(pct_cols))
  message(sprintf("Metacluster composition loaded: %d clones, %d metaclusters",
                  nrow(metacluster_comp), length(pct_cols)))
} else {
  metacluster_comp <- tibble(Global_Clone_ID = character())
  pct_cols <- character()
  warning("Metacluster composition CSV not found - metacluster plots will be skipped")
}


# ── 4. ASSAY SCORES AND WITHIN-DONOR PERCENTILES ─────────────────────────────
# Identical maths to section 14 of the original script.

clones_for_group <- function(prefix) {
  umi_summary %>%
    filter(str_starts(Global_Clone_ID, prefix)) %>%
    mutate(clone_num = as.integer(str_extract(Global_Clone_ID, "\\d+$"))) %>%
    arrange(clone_num) %>%
    pull(Global_Clone_ID)
}

find_well <- function(wells, pattern) {
  matched <- wells[str_detect(wells, pattern)]
  if (length(matched) == 0) return(NA_character_)
  matched[1]
}

pairwise_clr_diff <- function(umi_df, well_pos, well_neg) {
  if (!all(c(well_pos, well_neg) %in% colnames(umi_df))) return(NULL)

  pair <- umi_df %>% select(Global_Clone_ID, all_of(c(well_pos, well_neg)))

  # CLR per well independently, so each score is a pure pairwise contrast
  clr_pair <- pair %>%
    mutate(across(-Global_Clone_ID, ~ {
      x <- . + pseudocount
      as.numeric(clr(acomp(x / sum(x))))
    }))

  tibble(
    Global_Clone_ID = clr_pair$Global_Clone_ID,
    raw_score       = clr_pair[[well_pos]] - clr_pair[[well_neg]]
  )
}

scores_list <- list()

for (group_name in names(group_contrasts)) {
  grp        <- donor_groups[[group_name]]
  grp_wells  <- intersect(grp$wells, colnames(umi_summary))
  grp_clones <- clones_for_group(grp$clone_prefix)

  if (length(grp_clones) == 0 || length(grp_wells) == 0) {
    warning(sprintf("No clones or wells found for %s - skipping", group_name))
    next
  }

  grp_umi <- umi_summary %>%
    filter(Global_Clone_ID %in% grp_clones) %>%
    select(Global_Clone_ID, all_of(grp_wells))

  score_df <- tibble(Global_Clone_ID = grp_clones, Donor_Group = group_name)

  for (contrast in group_contrasts[[group_name]]) {
    well_pos <- find_well(grp_wells, contrast$pos)
    well_neg <- find_well(grp_wells, contrast$neg)

    if (is.na(well_pos) || is.na(well_neg)) {
      message(sprintf("  Skipping %s for %s (wells not found)",
                      contrast$name, group_name))
      next
    }

    result <- pairwise_clr_diff(grp_umi, well_pos, well_neg)

    score_df <- score_df %>%
      left_join(
        result %>% select(Global_Clone_ID,
                          !!paste0(contrast$name, "_score") := raw_score),
        by = "Global_Clone_ID"
      )
  }

  scores_list[[group_name]] <- score_df
  message(sprintf("  Computed assay scores for %s", group_name))
}

master_assay_scores <- bind_rows(scores_list)

# Within-donor range-based percentile: (score - min) / (max - min) * 100
# Kept so use_percentiles can be flipped without touching anything else.
score_cols <- names(master_assay_scores) %>% str_subset("_score$")

master_assay_scores <- master_assay_scores %>%
  group_by(Donor_Group) %>%
  mutate(across(
    all_of(score_cols),
    ~ {
      mn <- suppressWarnings(min(., na.rm = TRUE))
      mx <- suppressWarnings(max(., na.rm = TRUE))
      if (!is.finite(mn) || !is.finite(mx) || mx == mn) NA_real_
      else (. - mn) / (mx - mn) * 100
    },
    .names = "{.col}_percentile"
  )) %>%
  ungroup()


# ── 5. CONTROL-WELL FREQUENCY LOOKUP ─────────────────────────────────────────

control_freq_lookup <- bind_rows(
  lapply(names(donor_groups), function(g) {
    grp       <- donor_groups[[g]]
    wells     <- intersect(grp$wells, colnames(umi_summary))
    ctrl_well <- wells[str_detect(wells, "Control")][1]
    if (is.na(ctrl_well)) return(NULL)

    umi_summary %>%
      select(Global_Clone_ID, count = all_of(ctrl_well)) %>%
      mutate(ctrl_freq = count / sum(count), Donor_Group = g) %>%
      select(Donor_Group, Global_Clone_ID, ctrl_freq)
  })
)


# ── 6. ASSEMBLE ONE TABLE PER CLONE ──────────────────────────────────────────

clone_table <- master_assay_scores %>%
  left_join(control_freq_lookup, by = c("Donor_Group", "Global_Clone_ID")) %>%
  left_join(metacluster_comp,    by = "Global_Clone_ID") %>%
  mutate(
    ctrl_pct    = ctrl_freq * 100,
    donor_short = str_remove(Donor_Group, "_(pp65|IE1)$")
  )


# ── 7. PLOT BUILDER ──────────────────────────────────────────────────────────

safe_name <- function(x) {
  x %>% str_replace_all("[^A-Za-z0-9]+", "_") %>% str_remove("^_|_$")
}

# Padded range, used for axes with no fixed limits (i.e. raw scores)
padded_range <- function(values, pad_frac = 0.08) {
  v <- values[is.finite(values)]
  if (length(v) == 0) return(c(-1, 1))
  rng <- range(v)
  if (diff(rng) == 0) rng <- rng + c(-1, 1)
  pad <- diff(rng) * pad_frac
  c(rng[1] - pad, rng[2] + pad)
}

make_corr_plot <- function(plot_df, x_lab, y_lab, title,
                           x_lim, y_lim, x_breaks, y_breaks) {

  fit  <- lm(y ~ x, data = plot_df)
  r2   <- summary(fit)$r.squared
  pval <- tryCatch(summary(fit)$coefficients["x", "Pr(>|t|)"],
                   error = function(e) NA_real_)

  pval_text <- if (is.na(pval)) "p = NA" else
    if (pval < 0.001) "p < 0.001" else sprintf("p = %.3f", pval)
  annot <- sprintf("R² = %.3f\n%s", r2, pval_text)

  ggplot(plot_df, aes(x = x, y = y)) +
    geom_smooth(method = "lm", formula = y ~ x, se = se_ribbon,
                colour = "steelblue", fill = "grey85", linewidth = 0.8) +
    geom_point(aes(colour = donor_short), size = point_size, alpha = point_alpha) +
    annotate("text",
             x     = x_lim[1] + 0.03 * diff(x_lim),
             y     = y_lim[2] - 0.03 * diff(y_lim),
             hjust = 0, vjust = 1,
             label = annot,
             size  = annot_size,
             colour = "grey20") +
    scale_colour_manual(values = donor_colours, name = NULL, drop = TRUE) +
    guides(colour = guide_legend(override.aes = list(size = legend_point_size,
                                                     alpha = 1))) +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(breaks = y_breaks) +
    # coord_cartesian rather than scale limits so no point is dropped from the fit
    coord_cartesian(xlim = x_lim, ylim = y_lim) +
    labs(title = title, x = x_lab, y = y_lab) +
    theme_classic(base_size = base_size) +
    theme(
      plot.title         = element_text(size = title_size, face = "bold",
                                        lineheight = 1.05),
      axis.title         = element_text(size = axis_title_size),
      axis.text          = element_text(size = axis_text_size),
      legend.text        = element_text(size = legend_text_size),
      legend.position    = "right",
      # Legend pulled in toward the panel: one explicit gap instead of the
      # base_size-scaled default spacing plus the legend's own left margin,
      # and a key just wide enough for the dot so it doesn't float in padding
      legend.box.spacing = unit(legend_gap_pt, "pt"),
      legend.margin      = margin(t = half_line, r = half_line,
                                  b = half_line, l = 0),
      legend.key.width   = unit(legend_point_size * .pt + 3, "pt")
    )
}


# ── 8. PLOT SPECIFICATIONS ───────────────────────────────────────────────────
#
#   x = TCR score    vs  other assay scores
#   x = TCR score    vs  control-well frequency
#   x = TCR score    vs  % cells in each metacluster
#   x = CD107a score vs  % cells in each metacluster

score_suffix <- if (use_percentiles) "_score_percentile" else "_score"
score_limits <- if (use_percentiles) percentile_limits else NULL

build_specs <- function(available_cols) {

  specs <- list()

  add <- function(specs, x_col, x_short, x_lab, x_lim,
                  y_col, y_short, y_lab, y_lim, y_tag) {
    if (!all(c(x_col, y_col) %in% available_cols)) return(specs)
    append(specs, list(list(
      x_col = x_col, x_short = x_short, x_lab = x_lab, x_lim = x_lim,
      y_col = y_col, y_short = y_short, y_lab = y_lab, y_lim = y_lim,
      y_tag = y_tag
    )))
  }

  # (a) TCR vs other assay scores
  for (assay in c("Proliferation", "CD107a", "IFNg")) {
    specs <- add(specs,
                 paste0("TCR", score_suffix), "TCR",
                 score_axis_label("TCR"), score_limits,
                 paste0(assay, score_suffix), assay,
                 score_axis_label(assay), score_limits,
                 safe_name(paste0(assay, score_suffix)))
  }

  # (b) TCR vs clone dominance within the control well (%)
  specs <- add(specs,
               paste0("TCR", score_suffix), "TCR",
               score_axis_label("TCR"), score_limits,
               "ctrl_pct", "% Dominance",
               "% of clone within donor/specificity",
               pct_limits, "Control_frequency")

  # (c) and (d) TCR and CD107a vs metacluster composition
  for (x_assay in c("TCR", "CD107a")) {
    for (pct_col in pct_cols) {
      cluster_name <- str_remove(pct_col, "^pct_")
      specs <- add(specs,
                   paste0(x_assay, score_suffix), x_assay,
                   score_axis_label(x_assay), score_limits,
                   pct_col, paste0("% ", cluster_name),
                   paste0("% cells in ", cluster_name), pct_limits,
                   paste0("pct_", safe_name(cluster_name)))
    }
  }

  specs
}

plot_specs <- build_specs(names(clone_table))


# ── 9. GENERATE PLOTS ────────────────────────────────────────────────────────

stats_rows <- list()

for (grp_name in names(combined_groups)) {

  grp       <- combined_groups[[grp_name]]
  grp_df    <- clone_table %>% filter(Donor_Group %in% grp$donors)
  all_short <- unique(str_remove(grp$donors, "_(pp65|IE1)$"))

  if (nrow(grp_df) == 0) {
    warning(sprintf("No clones for combined group %s - skipping", grp_name))
    next
  }

  # One shared range per score column across every panel of this specificity
  shared_lims <- list()
  if (shared_score_limits) {
    for (cc in intersect(paste0(c("TCR", "Proliferation", "CD107a", "IFNg"),
                                score_suffix), names(grp_df))) {
      shared_lims[[cc]] <- padded_range(grp_df[[cc]])
    }
  }

  resolve_lim <- function(spec_lim, col, values) {
    if (!is.null(spec_lim)) return(spec_lim)
    if (!is.null(shared_lims[[col]])) return(shared_lims[[col]])
    padded_range(values)
  }

  for (spec in plot_specs) {

    plot_df <- grp_df %>%
      select(Global_Clone_ID, Donor_Group, donor_short,
             x = all_of(spec$x_col), y = all_of(spec$y_col)) %>%
      filter(!is.na(x), !is.na(y))

    if (nrow(plot_df) < 3) next

    x_lim <- resolve_lim(spec$x_lim, spec$x_col, plot_df$x)
    y_lim <- resolve_lim(spec$y_lim, spec$y_col, plot_df$y)

    # Fixed axes get evenly spaced breaks; free (raw score) axes use defaults
    x_breaks <- if (!is.null(spec$x_lim))
      seq(spec$x_lim[1], spec$x_lim[2], length.out = 5) else waiver()
    y_breaks <- if (!is.null(spec$y_lim))
      seq(spec$y_lim[1], spec$y_lim[2], length.out = 5) else waiver()

    # Which donors actually contributed data for this contrast
    donors_present <- plot_df %>% pull(donor_short) %>% unique() %>% sort()

    donor_tag <- if (length(donors_present) < length(all_short)) {
      paste0(" (", paste(donors_present, collapse = ", "), ")")
    } else ""

    title_line1 <- paste0(grp$label, donor_tag)
    title_line2 <- paste0(spec$x_short, " vs ", spec$y_short)
    title <- if (two_line_titles) {
      paste0(title_line1, "\n", title_line2)
    } else {
      paste0(title_line1, "  |  ", title_line2)
    }

    p <- make_corr_plot(
      plot_df  = plot_df,
      x_lab    = spec$x_lab,
      y_lab    = spec$y_lab,
      title    = title,
      x_lim    = x_lim,
      y_lim    = y_lim,
      x_breaks = x_breaks,
      y_breaks = y_breaks
    )

    filename <- paste0(grp$label, "_", safe_name(spec$x_short),
                       "_vs_", spec$y_tag, ".png")

    ggsave(
      filename = file.path(output_dir, filename),
      plot     = p,
      width    = plot_width,
      height   = plot_height,
      dpi      = plot_dpi,
      bg       = "white"
    )

    # Record fit statistics
    fit  <- lm(y ~ x, data = plot_df)
    pval <- tryCatch(summary(fit)$coefficients["x", "Pr(>|t|)"],
                     error = function(e) NA_real_)

    stats_rows[[length(stats_rows) + 1]] <- tibble(
      Specificity = grp$label,
      Donors      = paste(donors_present, collapse = ", "),
      X           = spec$x_short,
      Y           = spec$y_short,
      n_clones    = nrow(plot_df),
      R2          = summary(fit)$r.squared,
      p_value     = pval,
      File        = filename
    )

    message(sprintf("  Saved: %s (%d clones)", filename, nrow(plot_df)))
  }
}

if (length(stats_rows) > 0) {
  correlation_stats <- bind_rows(stats_rows)
  write.csv(correlation_stats,
            file.path(output_dir, "correlation_stats.csv"),
            row.names = FALSE)
  message(sprintf("\nWrote correlation_stats.csv (%d plots)",
                  nrow(correlation_stats)))
}

message("\nDone. Plots written to: ", output_dir)

