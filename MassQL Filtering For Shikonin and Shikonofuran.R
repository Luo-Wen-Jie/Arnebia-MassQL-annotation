packages <- c("readxl", "openxlsx", "dplyr", "stringr", "purrr", "tidyr")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

txt_file <- ""

table_s4_file <- ""

out_file <- ""

rt_tol <- 0.03
ppm_tol <- 10
frag_tol_da <- 0.01
min_required_ions <- 3
top_n <- 20

massql_rules <- list(
  Shikonin_Type_I = c(269.081, 251.070, 241.084),
  Shikonin_Type_II = c(285.076, 267.065, 257.082),
  Shikonin_Type_III = c(283.097, 268.073, 255.109),
  Shikonin_Type_IV = c(299.092, 281.082, 271.098),
  Shikonofuran_Type_I = c(255.109, 237.091, 227.107),
  Shikonofuran_Type_II = c(273.120, 255.102, 237.088)
)

ppm_error <- function(observed, theoretical) {
  abs(observed - theoretical) / theoretical * 1e6
}

parse_sample_compound <- function(x) {
  x <- as.character(x)
  
  rt <- str_match(x, "^\\s*([0-9]+\\.?[0-9]*)_")[, 2]
  mz <- str_match(x, "_([0-9]+\\.?[0-9]*)\\s*(m/z|mz|n|p)?")[, 2]
  
  tibble(
    Sample_Compound = x,
    Sample_RT = as.numeric(rt),
    Sample_mz = as.numeric(mz)
  )
}

normalize_id <- function(x) {
  x %>%
    as.character() %>%
    str_to_lower() %>%
    str_replace_all("\\s+", "") %>%
    str_replace_all("\\(|\\)", "") %>%
    str_replace_all("m/z", "mz")
}

extract_feature_from_text <- function(x) {
  m <- str_match(x, "([0-9]+\\.?[0-9]*)_([0-9]+\\.?[0-9]*)(m/z|mz|n|p)?")
  
  tibble(
    Feature_RT = as.numeric(m[, 2]),
    Feature_mz = as.numeric(m[, 3])
  )
}

parse_one_block <- function(block_lines, spec_id) {
  
  name_line <- block_lines[str_detect(block_lines, "^Name\\s*:")]
  comment_line <- block_lines[str_detect(block_lines, "^Comment\\s*:")]
  precursor_line <- block_lines[str_detect(block_lines, "^PrecursorMZ\\s*:")]
  
  name <- ifelse(
    length(name_line) > 0,
    str_replace(name_line[1], "^Name\\s*:\\s*", ""),
    NA_character_
  )
  
  comment <- ifelse(
    length(comment_line) > 0,
    str_replace(comment_line[1], "^Comment\\s*:\\s*", ""),
    NA_character_
  )
  
  precursor_mz <- ifelse(
    length(precursor_line) > 0,
    as.numeric(str_extract(precursor_line[1], "[0-9]+\\.?[0-9]*")),
    NA_real_
  )
  
  id_text <- paste(name, comment, sep = " ; ")
  feature_info <- extract_feature_from_text(id_text)
  
  peak_lines <- block_lines[str_detect(block_lines, "^\\s*[0-9]+\\.?[0-9]*\\s+[0-9Ee.+-]+")]
  
  if (length(peak_lines) > 0) {
    peaks <- tibble(raw = peak_lines) %>%
      mutate(
        mz = as.numeric(str_match(raw, "^\\s*([0-9]+\\.?[0-9]*)")[, 2]),
        intensity = as.numeric(str_match(raw, "^\\s*[0-9]+\\.?[0-9]*\\s+([0-9Ee.+-]+)")[, 2])
      ) %>%
      filter(!is.na(mz), !is.na(intensity)) %>%
      select(mz, intensity)
  } else {
    peaks <- tibble(mz = numeric(), intensity = numeric())
  }
  
  list(
    spectrum_info = tibble(
      Spectrum_ID = spec_id,
      Name = name,
      Comment = comment,
      PrecursorMZ = precursor_mz,
      Feature_RT = feature_info$Feature_RT[1],
      Feature_mz = feature_info$Feature_mz[1],
      ID_text = id_text,
      ID_norm = normalize_id(id_text)
    ),
    peaks = peaks %>%
      mutate(Spectrum_ID = spec_id, .before = 1)
  )
}

read_msp_like_txt <- function(txt_file) {
  
  lines <- readLines(txt_file, warn = FALSE, encoding = "UTF-8")
  lines <- str_replace_all(lines, "\r", "")
  
  name_idx <- which(str_detect(lines, "^Name\\s*:"))
  
  if (length(name_idx) == 0) {
    stop("No spectra found. The txt file does not contain lines starting with 'Name:'.")
  }
  
  end_idx <- c(name_idx[-1] - 1, length(lines))
  
  parsed <- map2(name_idx, end_idx, function(start, end) {
    block <- lines[start:end]
    parse_one_block(block, spec_id = start)
  })
  
  spectrum_info <- map_dfr(parsed, "spectrum_info")
  peak_table <- map_dfr(parsed, "peaks")
  
  list(
    spectrum_info = spectrum_info,
    peak_table = peak_table
  )
}

find_matched_spectrum <- function(sample_id, sample_rt, sample_mz, spectrum_info) {
  
  sample_norm <- normalize_id(sample_id)
  
  exact_hit <- spectrum_info %>%
    filter(str_detect(ID_norm, fixed(sample_norm))) %>%
    mutate(
      Match_method = "Exact ID match",
      RT_error = abs(Feature_RT - sample_rt),
      mz_error_ppm = ppm_error(Feature_mz, sample_mz),
      Match_score = 0
    )
  
  if (nrow(exact_hit) > 0) {
    return(exact_hit %>% slice(1))
  }
  
  numeric_hit <- spectrum_info %>%
    mutate(
      RT_error = abs(Feature_RT - sample_rt),
      Feature_mz_error_ppm = ppm_error(Feature_mz, sample_mz),
      Precursor_mz_error_ppm = ppm_error(PrecursorMZ, sample_mz),
      mz_error_ppm = pmin(Feature_mz_error_ppm, Precursor_mz_error_ppm, na.rm = TRUE)
    ) %>%
    filter(
      !is.na(RT_error),
      !is.na(mz_error_ppm),
      is.finite(mz_error_ppm),
      RT_error <= rt_tol,
      mz_error_ppm <= ppm_tol
    ) %>%
    mutate(
      Match_method = "RT + m/z match",
      Match_score = RT_error / rt_tol + mz_error_ppm / ppm_tol
    ) %>%
    arrange(Match_score)
  
  if (nrow(numeric_hit) > 0) {
    return(numeric_hit %>% slice(1))
  }
  
  tibble(
    Spectrum_ID = NA_integer_,
    Name = NA_character_,
    Comment = NA_character_,
    PrecursorMZ = NA_real_,
    Feature_RT = NA_real_,
    Feature_mz = NA_real_,
    ID_text = NA_character_,
    ID_norm = NA_character_,
    Match_method = "No matched spectrum",
    RT_error = NA_real_,
    mz_error_ppm = NA_real_,
    Match_score = NA_real_
  )
}

match_rule_one_spectrum <- function(peaks_top20, rule_name, target_ions, frag_tol_da) {
  
  if (nrow(peaks_top20) == 0) {
    return(tibble(
      Rule = rule_name,
      Target_ions = paste(target_ions, collapse = "; "),
      Matched_ion_count = 0,
      Matched_targets = NA_character_,
      Matched_observed_mz = NA_character_,
      Matched_delta_Da = NA_character_,
      Pass = FALSE
    ))
  }
  
  matched_list <- map_dfr(target_ions, function(target) {
    
    hit <- peaks_top20 %>%
      mutate(
        Target_mz = target,
        Delta_Da = abs(mz - target)
      ) %>%
      filter(Delta_Da <= frag_tol_da) %>%
      arrange(Delta_Da, desc(intensity)) %>%
      slice(1)
    
    if (nrow(hit) == 0) {
      tibble(
        Target_mz = target,
        Observed_mz = NA_real_,
        Delta_Da = NA_real_,
        intensity = NA_real_
      )
    } else {
      tibble(
        Target_mz = target,
        Observed_mz = hit$mz[1],
        Delta_Da = hit$Delta_Da[1],
        intensity = hit$intensity[1]
      )
    }
  })
  
  matched <- matched_list %>%
    filter(!is.na(Observed_mz))
  
  tibble(
    Rule = rule_name,
    Target_ions = paste(target_ions, collapse = "; "),
    Matched_ion_count = nrow(matched),
    Matched_targets = ifelse(
      nrow(matched) > 0,
      paste(round(matched$Target_mz, 4), collapse = "; "),
      NA_character_
    ),
    Matched_observed_mz = ifelse(
      nrow(matched) > 0,
      paste(round(matched$Observed_mz, 4), collapse = "; "),
      NA_character_
    ),
    Matched_delta_Da = ifelse(
      nrow(matched) > 0,
      paste(round(matched$Delta_Da, 5), collapse = "; "),
      NA_character_
    ),
    Pass = nrow(matched) >= min_required_ions
  )
}

table_s4 <- read_excel(table_s4_file, sheet = "SHK and SHF")

sample_col <- names(table_s4)[str_trim(str_to_lower(names(table_s4))) == "sample compound"]

if (length(sample_col) == 0) {
  stop("Cannot find column named 'Sample Compound' in sheet 'SHK and SHF'.")
}

table_s4 <- table_s4 %>%
  mutate(Row_ID = row_number(), .before = 1)

sample_info <- parse_sample_compound(table_s4[[sample_col]])

table_s4_parsed <- bind_cols(table_s4, sample_info %>% select(Sample_RT, Sample_mz))

txt_data <- read_msp_like_txt(txt_file)

spectrum_info <- txt_data$spectrum_info
peak_table <- txt_data$peak_table

matched_spectra <- pmap_dfr(
  list(
    sample_id = table_s4_parsed[[sample_col]],
    sample_rt = table_s4_parsed$Sample_RT,
    sample_mz = table_s4_parsed$Sample_mz
  ),
  function(sample_id, sample_rt, sample_mz) {
    find_matched_spectrum(sample_id, sample_rt, sample_mz, spectrum_info)
  }
) %>%
  mutate(Row_ID = table_s4_parsed$Row_ID, .before = 1)

table_s4_matched <- table_s4_parsed %>%
  left_join(
    matched_spectra %>%
      select(
        Row_ID,
        Spectrum_ID,
        Match_method,
        Name,
        Comment,
        PrecursorMZ,
        Feature_RT,
        Feature_mz,
        RT_error,
        mz_error_ppm
      ),
    by = "Row_ID"
  )

top20_ms2 <- table_s4_matched %>%
  select(Row_ID, all_of(sample_col), Spectrum_ID) %>%
  filter(!is.na(Spectrum_ID)) %>%
  left_join(peak_table, by = "Spectrum_ID") %>%
  group_by(Row_ID) %>%
  arrange(desc(intensity), .by_group = TRUE) %>%
  slice_head(n = top_n) %>%
  mutate(Fragment_rank_by_intensity = row_number()) %>%
  ungroup() %>%
  arrange(Row_ID, Fragment_rank_by_intensity)

rule_detail <- map_dfr(table_s4_matched$Row_ID, function(rid) {
  
  peaks_i <- top20_ms2 %>%
    filter(Row_ID == rid) %>%
    select(mz, intensity, Fragment_rank_by_intensity)
  
  map_dfr(names(massql_rules), function(rule_name) {
    match_rule_one_spectrum(
      peaks_top20 = peaks_i,
      rule_name = rule_name,
      target_ions = massql_rules[[rule_name]],
      frag_tol_da = frag_tol_da
    )
  }) %>%
    mutate(Row_ID = rid, .before = 1)
})

filter_summary <- rule_detail %>%
  group_by(Row_ID) %>%
  summarise(
    MassQL_Pass = any(Pass),
    Passed_rules = ifelse(
      any(Pass),
      paste(Rule[Pass], collapse = "; "),
      NA_character_
    ),
    Best_matched_ion_count = max(Matched_ion_count, na.rm = TRUE),
    Matched_diagnostic_ions = paste(
      unique(na.omit(Matched_observed_mz[Matched_ion_count > 0])),
      collapse = " | "
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Filter_result = ifelse(MassQL_Pass, "Pass", "Fail")
  )

final_result <- table_s4_matched %>%
  left_join(filter_summary, by = "Row_ID") %>%
  mutate(
    MassQL_Pass = ifelse(is.na(MassQL_Pass), FALSE, MassQL_Pass),
    Filter_result = case_when(
      Match_method == "No matched spectrum" ~ "No MS/MS found",
      MassQL_Pass ~ "Pass",
      TRUE ~ "Fail"
    )
  )

summary_table <- tibble(
  Item = c(
    "Input SHK/SHF rows",
    "Rows with matched MS/MS spectrum",
    "Rows without matched MS/MS spectrum",
    "Rows after MassQL-like filtering",
    "Rows failed after filtering",
    "Top N fragments retained",
    "Fragment tolerance Da",
    "Minimum required diagnostic ions"
  ),
  Value = c(
    nrow(final_result),
    sum(final_result$Match_method != "No matched spectrum", na.rm = TRUE),
    sum(final_result$Match_method == "No matched spectrum", na.rm = TRUE),
    sum(final_result$Filter_result == "Pass", na.rm = TRUE),
    sum(final_result$Filter_result == "Fail", na.rm = TRUE),
    top_n,
    frag_tol_da,
    min_required_ions
  )
)

missing_spectra <- final_result %>%
  filter(Match_method == "No matched spectrum")

openxlsx::write.xlsx(
  list(
    "Filter result" = final_result,
    "Top20 MS2" = top20_ms2,
    "Rule detail" = rule_detail,
    "Missing spectra" = missing_spectra,
    "Summary" = summary_table
  ),
  file = out_file,
  overwrite = TRUE
)

cat("Done!\n")
cat("Output file:\n", out_file, "\n")