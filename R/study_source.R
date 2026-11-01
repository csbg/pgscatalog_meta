# Display names and cluster labels for the studies behind a sample set.
# normalize_study() is the name shown in a table. The cluster label is what
# metafor::vcalc groups on. A missing cohort is unnamed:<sample-set id>,
# one token per sample set, so two unnamed sample sets stay independent.

normalize_study <- function(x) {
  x0 <- trimws(as.character(x))
  x0[is.na(x) | !nzchar(x0)] <- "(blank)"
  xl <- tolower(x0)
  out <- x0
  out[xl %in% c("ukb", "uk biobank") | grepl("^uk\\s*biobank", xl)] <- "UK Biobank"
  out[xl %in% c("g&h", "g and h") | grepl("genes\\s*&\\s*health|genes and health", xl)] <- "Genes & Health"
  out[xl %in% c("allofus", "all of us") | grepl("all\\s*of\\s*us", xl)] <- "All of Us"
  out
}

cohort_is_unnamed <- function(x) {
  x0 <- trimws(as.character(x))
  is.na(x) | !nzchar(x0) | tolower(x0) %in% c("na", "(blank)")
}

sample_set_source_token <- function(cohort_txt, sampleset_id) {
  sampleset_id <- as.character(sampleset_id)
  unnamed <- cohort_is_unnamed(cohort_txt)
  if (any(unnamed & !nzchar(sampleset_id))) {
    stop("An unnamed cohort has no sample-set id.", call. = FALSE)
  }
  token <- normalize_study(cohort_txt)
  token[unnamed] <- paste0("unnamed:", sampleset_id[unnamed])
  token
}

source_cluster_label <- function(tokens) {
  tokens <- unique(as.character(tokens))
  tokens <- tokens[!is.na(tokens) & nzchar(tokens)]
  if (!length(tokens)) stop("A score side has no source token.", call. = FALSE)
  paste(sort(tokens), collapse = " + ")
}

component_ids <- function(n, edges) {
  if (n == 0L) return(integer())
  parent <- seq_len(n)
  find_root <- function(i) {
    if (parent[[i]] != i) parent[[i]] <<- find_root(parent[[i]])
    parent[[i]]
  }
  if (nrow(edges) > 0L) {
    for (k in seq_len(nrow(edges))) {
      a <- find_root(edges[k, 1L])
      b <- find_root(edges[k, 2L])
      if (a != b) parent[[b]] <- a
    }
  }
  roots <- vapply(seq_len(n), find_root, integer(1L))
  match(roots, unique(roots))
}

pair_edge_matrix <- function(sets) {
  n <- length(sets)
  empty <- matrix(integer(), ncol = 2L)
  if (n < 2L) return(empty)
  ii <- integer()
  jj <- integer()
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      if (length(intersect(sets[[i]], sets[[j]])) > 0L) {
        ii <- c(ii, i)
        jj <- c(jj, j)
      }
    }
  }
  if (!length(ii)) empty else cbind(ii, jj)
}

group_labels <- function(ids, comp) {
  ids <- as.character(ids)
  labs <- tapply(ids, comp, function(x) paste(sort(unique(x)), collapse = "+"))
  as.vector(unname(labs[as.character(comp)]))
}

# pairs: one row per score in a comparison, with the delta-method standard errors.
# keys: sample-set ids that entered that score and covariate class.
# sample_sets: sampleset_id and cohort_txt from the catalog sample-set table.
build_source_sharing <- function(pairs, keys, sample_sets) {
  need <- c(
    "trait_label", "trained_bucket", "target_ancestry", "pgs_id", "covariate_class",
    "delta", "delta_se", "eur_se", "tgt_se"
  )
  missing <- setdiff(need, names(pairs))
  if (length(missing)) stop("pairs is missing: ", paste(missing, collapse = ", "), call. = FALSE)

  eur <- dplyr::transmute(
    pairs,
    trait_label, trained_bucket, target_ancestry, pgs_id, covariate_class,
    side = "European", ancestry = "European"
  )
  tgt <- dplyr::transmute(
    pairs,
    trait_label, trained_bucket, target_ancestry, pgs_id, covariate_class,
    side = "Target", ancestry = .data$target_ancestry
  )
  sides <- dplyr::bind_rows(eur, tgt)
  long <- dplyr::inner_join(
    sides,
    dplyr::select(keys, "pgs_id", "trait_label", "ancestry_display", "covariate_class", "sampleset_id"),
    by = c("pgs_id", "trait_label", "ancestry" = "ancestry_display", "covariate_class"),
    relationship = "many-to-many"
  )
  if (nrow(long) == 0L) stop("No sample sets joined to the paired scores.", call. = FALSE)

  present <- sides |>
    dplyr::distinct(.data$trait_label, .data$trained_bucket, .data$target_ancestry, .data$pgs_id, .data$side)
  got <- long |>
    dplyr::distinct(.data$trait_label, .data$trained_bucket, .data$target_ancestry, .data$pgs_id, .data$side)
  if (nrow(present) != nrow(got)) {
    stop("A score side has no sample set in the kept evaluations.", call. = FALSE)
  }

  sample_sets <- dplyr::distinct(sample_sets, .data$sampleset_id, .keep_all = TRUE)
  unknown <- setdiff(unique(long$sampleset_id), sample_sets$sampleset_id)
  if (length(unknown)) {
    stop("Sample set missing from the catalog table: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  long <- dplyr::left_join(
    long,
    dplyr::select(sample_sets, "sampleset_id", "cohort_txt"),
    by = "sampleset_id"
  )
  long$source_token <- sample_set_source_token(long$cohort_txt, long$sampleset_id)
  long$display_name <- normalize_study(long$cohort_txt)
  long$cohort_show <- ifelse(
    is.na(long$cohort_txt) | !nzchar(trimws(as.character(long$cohort_txt))),
    "NA",
    trimws(as.character(long$cohort_txt))
  )

  long <- dplyr::distinct(
    long,
    .data$trait_label, .data$trained_bucket, .data$target_ancestry,
    .data$pgs_id, .data$side, .data$sampleset_id,
    .keep_all = TRUE
  )

  scores <- long |>
    dplyr::arrange(.data$sampleset_id, .by_group = FALSE) |>
    dplyr::group_by(.data$trait_label, .data$trained_bucket, .data$target_ancestry, .data$pgs_id, .data$side) |>
    dplyr::summarise(
      sampleset_ids = paste(.data$sampleset_id, collapse = ";"),
      cohort_txt = paste(.data$cohort_show, collapse = ";"),
      display_name = paste(.data$display_name, collapse = ";"),
      cluster_label = source_cluster_label(.data$source_token),
      tokens = list(sort(unique(.data$source_token))),
      pss = list(sort(unique(.data$sampleset_id))),
      .groups = "drop"
    )

  eur_s <- scores |>
    dplyr::filter(.data$side == "European") |>
    dplyr::select(
      "trait_label", "trained_bucket", "target_ancestry", "pgs_id",
      eur_cluster = "cluster_label",
      eur_tokens = "tokens",
      eur_pss = "pss",
      eur_sampleset_ids = "sampleset_ids",
      eur_cohort = "cohort_txt",
      eur_display = "display_name"
    )
  tgt_s <- scores |>
    dplyr::filter(.data$side == "Target") |>
    dplyr::select(
      "trait_label", "trained_bucket", "target_ancestry", "pgs_id",
      tgt_cluster = "cluster_label",
      tgt_tokens = "tokens",
      tgt_pss = "pss",
      tgt_sampleset_ids = "sampleset_ids",
      tgt_cohort = "cohort_txt",
      tgt_display = "display_name"
    )
  wide <- dplyr::inner_join(
    eur_s, tgt_s,
    by = c("trait_label", "trained_bucket", "target_ancestry", "pgs_id")
  )
  wide <- dplyr::inner_join(
    wide,
    dplyr::select(
      pairs,
      "trait_label", "trained_bucket", "target_ancestry", "pgs_id",
      "delta", "delta_se", "eur_se", "tgt_se"
    ),
    by = c("trait_label", "trained_bucket", "target_ancestry", "pgs_id")
  )
  if (nrow(wide) != nrow(pairs)) {
    stop("Source labels do not cover every paired score.", call. = FALSE)
  }

  pieces <- dplyr::group_split(wide, .data$trait_label, .data$trained_bucket, .data$target_ancestry)
  summaries <- vector("list", length(pieces))
  model_parts <- vector("list", length(pieces))
  for (i in seq_along(pieces)) {
    one <- sharing_one_comparison(pieces[[i]], long)
    summaries[[i]] <- one$summary
    model_parts[[i]] <- one$scores
  }

  list(
    comparisons = dplyr::bind_rows(summaries),
    scores = dplyr::select(
      scores,
      "trait_label", "trained_bucket", "target_ancestry", "pgs_id", "side",
      "sampleset_ids", "cohort_txt", "display_name", "cluster_label"
    ),
    model_scores = dplyr::bind_rows(model_parts)
  )
}

sharing_one_comparison <- function(df, long) {
  df <- df[order(df$pgs_id), , drop = FALSE]
  n <- nrow(df)
  ids <- df$pgs_id
  eur_edges <- pair_edge_matrix(df$eur_tokens)
  tgt_edges <- pair_edge_matrix(df$tgt_tokens)
  either <- rbind(eur_edges, tgt_edges)
  if (nrow(either) > 0L) either <- unique(either)
  sharing <- component_ids(n, if (nrow(either)) either else matrix(integer(), ncol = 2L))
  eur_comp <- component_ids(n, eur_edges)
  tgt_comp <- component_ids(n, tgt_edges)

  n_partial <- 0L
  eur_shared <- FALSE
  tgt_shared <- FALSE
  pss_eur <- FALSE
  pss_tgt <- FALSE
  source_adds_link <- FALSE
  if (n >= 2L) {
    for (i in seq_len(n - 1L)) {
      for (j in (i + 1L):n) {
        eur_tok <- length(intersect(df$eur_tokens[[i]], df$eur_tokens[[j]])) > 0L
        tgt_tok <- length(intersect(df$tgt_tokens[[i]], df$tgt_tokens[[j]])) > 0L
        eur_id <- length(intersect(df$eur_pss[[i]], df$eur_pss[[j]])) > 0L
        tgt_id <- length(intersect(df$tgt_pss[[i]], df$tgt_pss[[j]])) > 0L
        eur_shared <- eur_shared || eur_tok
        tgt_shared <- tgt_shared || tgt_tok
        pss_eur <- pss_eur || eur_id
        pss_tgt <- pss_tgt || tgt_id
        partial <- (eur_tok && df$eur_cluster[[i]] != df$eur_cluster[[j]]) ||
          (tgt_tok && df$tgt_cluster[[i]] != df$tgt_cluster[[j]])
        if (partial) n_partial <- n_partial + 1L
        if ((eur_tok && !eur_id) || (tgt_tok && !tgt_id)) source_adds_link <- TRUE
      }
    }
  }

  share_side <- if (eur_shared && tgt_shared) {
    "both"
  } else if (eur_shared) {
    "European"
  } else if (tgt_shared) {
    "Target"
  } else {
    "neither"
  }

  cohorts <- long$cohort_show[
    long$trait_label == df$trait_label[[1]] &
      long$trained_bucket == df$trained_bucket[[1]] &
      long$target_ancestry == df$target_ancestry[[1]]
  ]
  consortia <- sort(unique(cohorts[grepl("|", cohorts, fixed = TRUE)]))

  summary <- tibble::tibble(
    trait_label = df$trait_label[[1]],
    trained_bucket = df$trained_bucket[[1]],
    target_ancestry = df$target_ancestry[[1]],
    n_scores = n,
    eur_shared_sample_set_id = pss_eur,
    tgt_shared_sample_set_id = pss_tgt,
    eur_shared_source = eur_shared,
    tgt_shared_source = tgt_shared,
    share_side = share_side,
    source_adds_link = source_adds_link,
    n_sharing_groups = length(unique(sharing)),
    n_partial_pairs = n_partial,
    consortium_strings = paste(consortia, collapse = "; "),
    goes_through_covariance_model = eur_shared || tgt_shared
  )
  scores <- df
  scores$sharing_group <- group_labels(ids, sharing)
  scores$eur_group <- group_labels(ids, eur_comp)
  scores$tgt_group <- group_labels(ids, tgt_comp)
  scores$eur_tokens <- NULL
  scores$tgt_tokens <- NULL
  scores$eur_pss <- NULL
  scores$tgt_pss <- NULL
  list(summary = summary, scores = scores)
}

unnamed_token_collisions <- function(scores) {
  # A token unnamed:<id> may appear on many scores that use that sample set.
  # It must never be the token of a different sample set.
  hits <- 0L
  for (i in seq_len(nrow(scores))) {
    ids <- strsplit(scores$sampleset_ids[[i]], ";", fixed = TRUE)[[1]]
    labs <- strsplit(scores$cluster_label[[i]], " + ", fixed = TRUE)[[1]]
    unnamed <- labs[startsWith(labs, "unnamed:")]
    for (tok in unnamed) {
      pss <- sub("^unnamed:", "", tok)
      if (!pss %in% ids) hits <- hits + 1L
    }
  }
  hits
}
