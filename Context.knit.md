---
title: "ICP - NOC Nanopore Tree of Life Metabarcoding - analysis report"
author: "Nathan"
output-file: index.html
format:
  html:
    toc: true
    toc-depth: 3
    theme: cosmo
    code-fold: true
    code-block-bg: true
editor: source
execute:
  echo: TRUE
  warning: false
  message: false
---

# 🧭 Introduction

This project evaluates the use of Oxford Nanopore Technologies (ONT) sequencing combined with a Tree of Life metabarcoding approach for broad-scale marine biodiversity assessment. Environmental DNA (eDNA) samples were collected using the **RoCSI autonomous sampling system** and subsequently analysed using the **portable MinION sequencing platform**, providing a workflow that combines autonomous marine sampling with rapid, portable molecular analysis.

<div style="display:flex; justify-content:center; gap:3rem; align-items:flex-start; margin:1.5rem auto;">

<div style="width:35%; text-align:center;"> <img src="figures/RoCSI.jpeg" style="width:100%; height:auto;"> <p style="font-size:0.9em; margin-top:0.4rem;">The NOC’s Robotic Cartridge Sampling Instrument (RoCSI).</p> </div>

<div style="width:35%; text-align:center;"> <img src="figures/MinION.jpg" style="width:100%; height:auto;"> <p style="font-size:0.9em; margin-top:0.4rem;">The Oxford Nanopore MinION sequencing device.</p> </div>

</div>

The study uses 33 seawater samples collected during the **BIO-Carbon II** cruise aboard the *RRS James Cook*, along a transect extending from northern Scotland to Cardiff and covering much of the west coast of Great Britain. The transect spans contrasting marine environments and substantial gradients in environmental conditions, including temperature, salinity and chlorophyll concentration, providing a valuable test bed for assessing the ability of ONT metabarcoding to characterise biodiversity across different marine habitats.

![](figures/Map.png){fig-cap="Map of UK showing the chlorophyll concentrations and the location of the samples" width="100%"
fig-align="center"
}

By integrating portable nanopore sequencing with autonomous eDNA sampling, the project also contributes to the development of faster and more scalable approaches to marine biodiversity monitoring. Ultimately, such workflows could support the transition from conventional retrospective eDNA analysis towards near-real-time biodiversity observations at sea.

## Metadata plots

Before looking into the metabarcoding data, let's have a quick look at the metadata.


::: {.cell}

```{.r .cell-code}
library(dplyr)
library(readxl)
library(tidyverse)
library(stringr)
library(tibble)
library(Biostrings)
library(ggOceanMaps)
library(ggspatial)

# Metadata → rownames = barcodeXX
meta <- read_xlsx("data/Metadata_JC269_UKWC.xlsx")

meta <- meta %>%
  mutate(
    Lon = as.numeric(Lon),
    Lat = as.numeric(Lat),
    Temperature = as.numeric(Temperature),
    Fluorescence = as.numeric(Fluorescence)
  ) %>%
  pivot_longer(
    cols = c(Temperature, Fluorescence),
    names_to = "Variable",
    values_to = "Value"
  )


ggplot(meta, aes(x = Lat, y = Value)) +
  geom_point(size = 2.5, colour = "steelblue") +
  geom_line(colour = "steelblue", alpha = 0.6) +
  facet_wrap(~Variable, ncol = 1, scales = "free_y") +
  labs(
    x = "Longitude (°)",
    y = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
```

::: {.cell-output-display}
![](Context_files/figure-html/unnamed-chunk-1-1.png){width=672}
:::
:::



::: {.cell}

```{.r .cell-code}
library(ggOceanMaps)
library(ggspatial)

# Define map limits
lon1 <- min(meta$Lon, na.rm = TRUE) - 2
lon2 <- max(meta$Lon, na.rm = TRUE) + 2
lat1 <- min(meta$Lat, na.rm = TRUE) - 2
lat2 <- max(meta$Lat, na.rm = TRUE) + 2

# Sample locations
stations <- meta %>%
  distinct(ID, Lon, Lat, Station)

# Plot
basemap(
  limits = c(lon1, lon2, lat1, lat2),
  bathy.style = NULL
) +
  geom_spatial_point(
    data = stations,
    aes(x = Lon, y = Lat),
    color = "red",
    size = 2
  ) +
  coord_sf(expand = FALSE) +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_bw(base_size = 13)
```

::: {.cell-output-display}
![](Context_files/figure-html/map-1.png){width=672}
:::
:::


# 🧬 Sequencing and bioinformatics overview


Six metabarcoding markers were used to characterise complementary components of marine eDNA biodiversity. For each marker, **33 environmental samples and two negative controls** (extraction and PCR controls) were analysed, representing **210 sequencing libraries**.

- **16S rRNA V4:** Bacteria and archaea
- **18S rRNA V4:** Eukaryotes, including protists and phytoplankton
- **COI:** Metazoans
- **12S rRNA (MiFish):** Fish
- **ITS:** Primarily fungi
- **rbcL:** Photosynthetic organisms, particularly diatoms

Amplicons were sequenced at the **Centre for Genomic Research (CGR), University of Liverpool**, using Illumina MiSeq i100 paired-end sequencing (2 × 300 bp).

## Metabarcoding markers and primers


::: {.cell}

```{.r .cell-code}
library(knitr)
library(kableExtra)

primers <- tibble::tribble(
  ~Marker, ~`Fragment (bp)`, ~Primer, ~`Sequence (5′–3′)`, ~Reference,
  "16S V4", "~252", "515F_mod_CGR", "GTGYCAGCMGCCGCGGTAA", "Parada et al.",
  "16S V4", "~252", "806R_mod_CGR", "GGACTACNVGGGTWTCTAAT", "Apprill et al.",
  "18S V4", "~400", "TAReuk454FWD1_CGR", "CCAGCASCYGCGGTAATTCC", "Stoeck et al., 2010",
  "18S V4", "~400", "TAReukREV3_CGR", "ACTTTCGTTCTTGATYRA", "Stoeck et al., 2010",
  "COI", "~313", "mlCOIintF_CGR", "GGWACWGGWTGAACWGTWTAYCCYCC", "Leray, 2013",
  "COI", "~313", "jgHCO2198_CGR", "TAIACYTCIGGRTGICCRAARAAYCA", "Folmer et al., 1994",
  "MiFish (12S)", "~170", "MiFish-U-F_CGR", "GTCGGTAAAACTCGTGCCAGC", "Miya et al., 2015",
  "MiFish (12S)", "~170", "MiFish-U-R_mod_CGR", "CATAGTGGGGTATCTAATCCCAGTTTGT", "Stoeckle et al., 2022",
  "ITS2", "~250", "ITS3", "GCATCGATGAAGAACGCAGC", "White et al., 1990",
  "ITS2", "~250", "ITS4", "TCCTCCGCTTATTGATATGC", "White et al., 1990",
  "rbcL", "~312", "708F-DEG", "AGGTGAAGYWAAAGGTTCWTAYTTAAA", "Vasselon et al., 2017",
  "rbcL", "~312", "R3-DEG", "CCTTCTAATTTACCWACWACWG", "Vasselon et al., 2017"
)

primers |>
  kbl(
    format = "html",
    escape = TRUE,
    align = c("l", "c", "l", "l", "l")
  ) |>
  kable_styling(
    bootstrap_options = c("striped", "condensed"),
    full_width = FALSE,
    position = "center"
  ) |>
  column_spec(1, width = "115px", bold = TRUE) |>
  column_spec(2, width = "95px") |>
  column_spec(3, width = "195px") |>
  column_spec(4, width = "285px", extra_css = "font-family: monospace; font-size: 0.9em;") |>
  column_spec(5, width = "150px") |>
  collapse_rows(
    columns = 1:2,
    valign = "middle",
    row_group_label_position = "stack"
  )
```

::: {.cell-output-display}
`````{=html}
<table class="table table-striped table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> Marker </th>
   <th style="text-align:center;"> Fragment (bp) </th>
   <th style="text-align:left;"> Primer </th>
   <th style="text-align:left;"> Sequence (5′–3′) </th>
   <th style="text-align:left;"> Reference </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;width: 115px; font-weight: bold;vertical-align: middle !important;" rowspan="2"> 16S V4 </td>
   <td style="text-align:center;width: 95px; vertical-align: middle !important;" rowspan="2"> ~252 </td>
   <td style="text-align:left;width: 195px; "> 515F_mod_CGR </td>
   <td style="text-align:left;width: 285px; font-family: monospace; font-size: 0.9em;"> GTGYCAGCMGCCGCGGTAA </td>
   <td style="text-align:left;width: 150px; "> Parada et al. </td>
  </tr>
  <tr>
   
   
   <td style="text-align:left;width: 195px; "> 806R_mod_CGR </td>
   <td style="text-align:left;width: 285px; font-family: monospace; font-size: 0.9em;"> GGACTACNVGGGTWTCTAAT </td>
   <td style="text-align:left;width: 150px; "> Apprill et al. </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 115px; font-weight: bold;vertical-align: middle !important;" rowspan="2"> 18S V4 </td>
   <td style="text-align:center;width: 95px; vertical-align: middle !important;" rowspan="2"> ~400 </td>
   <td style="text-align:left;width: 195px; "> TAReuk454FWD1_CGR </td>
   <td style="text-align:left;width: 285px; font-family: monospace; font-size: 0.9em;"> CCAGCASCYGCGGTAATTCC </td>
   <td style="text-align:left;width: 150px; "> Stoeck et al., 2010 </td>
  </tr>
  <tr>
   
   
   <td style="text-align:left;width: 195px; "> TAReukREV3_CGR </td>
   <td style="text-align:left;width: 285px; font-family: monospace; font-size: 0.9em;"> ACTTTCGTTCTTGATYRA </td>
   <td style="text-align:left;width: 150px; "> Stoeck et al., 2010 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 115px; font-weight: bold;vertical-align: middle !important;" rowspan="2"> COI </td>
   <td style="text-align:center;width: 95px; vertical-align: middle !important;" rowspan="2"> ~313 </td>
   <td style="text-align:left;width: 195px; "> mlCOIintF_CGR </td>
   <td style="text-align:left;width: 285px; font-family: monospace; font-size: 0.9em;"> GGWACWGGWTGAACWGTWTAYCCYCC </td>
   <td style="text-align:left;width: 150px; "> Leray, 2013 </td>
  </tr>
  <tr>
   
   
   <td style="text-align:left;width: 195px; "> jgHCO2198_CGR </td>
   <td style="text-align:left;width: 285px; font-family: monospace; font-size: 0.9em;"> TAIACYTCIGGRTGICCRAARAAYCA </td>
   <td style="text-align:left;width: 150px; "> Folmer et al., 1994 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 115px; font-weight: bold;vertical-align: middle !important;" rowspan="2"> MiFish (12S) </td>
   <td style="text-align:center;width: 95px; vertical-align: middle !important;" rowspan="2"> ~170 </td>
   <td style="text-align:left;width: 195px; "> MiFish-U-F_CGR </td>
   <td style="text-align:left;width: 285px; font-family: monospace; font-size: 0.9em;"> GTCGGTAAAACTCGTGCCAGC </td>
   <td style="text-align:left;width: 150px; "> Miya et al., 2015 </td>
  </tr>
  <tr>
   
   
   <td style="text-align:left;width: 195px; "> MiFish-U-R_mod_CGR </td>
   <td style="text-align:left;width: 285px; font-family: monospace; font-size: 0.9em;"> CATAGTGGGGTATCTAATCCCAGTTTGT </td>
   <td style="text-align:left;width: 150px; "> Stoeckle et al., 2022 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 115px; font-weight: bold;vertical-align: middle !important;" rowspan="2"> ITS2 </td>
   <td style="text-align:center;width: 95px; vertical-align: middle !important;" rowspan="2"> ~250 </td>
   <td style="text-align:left;width: 195px; "> ITS3 </td>
   <td style="text-align:left;width: 285px; font-family: monospace; font-size: 0.9em;"> GCATCGATGAAGAACGCAGC </td>
   <td style="text-align:left;width: 150px; "> White et al., 1990 </td>
  </tr>
  <tr>
   
   
   <td style="text-align:left;width: 195px; "> ITS4 </td>
   <td style="text-align:left;width: 285px; font-family: monospace; font-size: 0.9em;"> TCCTCCGCTTATTGATATGC </td>
   <td style="text-align:left;width: 150px; "> White et al., 1990 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 115px; font-weight: bold;vertical-align: middle !important;" rowspan="2"> rbcL </td>
   <td style="text-align:center;width: 95px; vertical-align: middle !important;" rowspan="2"> ~312 </td>
   <td style="text-align:left;width: 195px; "> 708F-DEG </td>
   <td style="text-align:left;width: 285px; font-family: monospace; font-size: 0.9em;"> AGGTGAAGYWAAAGGTTCWTAYTTAAA </td>
   <td style="text-align:left;width: 150px; "> Vasselon et al., 2017 </td>
  </tr>
  <tr>
   
   
   <td style="text-align:left;width: 195px; "> R3-DEG </td>
   <td style="text-align:left;width: 285px; font-family: monospace; font-size: 0.9em;"> CCTTCTAATTTACCWACWACWG </td>
   <td style="text-align:left;width: 150px; "> Vasselon et al., 2017 </td>
  </tr>
</tbody>
</table>

`````
:::
:::


**Library preparation:** CGR overhang sequences were added to the 5′ ends of the forward (`ACACTCTTTCCCTACACGACGCTCTTCCGATCT`) and reverse (`GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT`) primers. These overhangs are excluded from the sequences in the table.

## Bioinformatics processing parameters

Primer trimming was performed with **Cutadapt v5.2**, followed by ASV inference and taxonomic assignment using **DADA2**.

For all markers, Cutadapt used an **error rate of 0.1**, minimum primer **overlap of 15 bp**, **no indels** and `--discard-untrimmed`. Both forward and reverse primer matches were required.


::: {.cell}

```{.r .cell-code}
params <- tribble(
  ~Marker, ~`Truncation R1/R2 (bp)`, ~`Taxonomic reference`,
  "16S",    "250 / 250", "SILVA 138.2",
  "18S",    "270 / 255", "PR2 5.1.0",
  "COI",    "260 / 260", "MIDORI2 GB273",
  "MiFish", "None",      "MitoFish (Jul 2026)",
  "ITS",    "270 / 270", "UNITE (Feb 2025)",
  "rbcL",   "250 / 250", "Diat.barcode 16.3"
)

params |>
  kbl(
    format = "html",
    align = c("l", "c", "l"),
    escape = TRUE
  ) |>
  kable_styling(
    bootstrap_options = c("striped", "condensed"),
    full_width = FALSE,
    position = "center"
  ) |>
  column_spec(1, bold = TRUE, width = "90px") |>
  column_spec(2, width = "190px") |>
  column_spec(3, width = "190px")
```

::: {.cell-output-display}
`````{=html}
<table class="table table-striped table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> Marker </th>
   <th style="text-align:center;"> Truncation R1/R2 (bp) </th>
   <th style="text-align:left;"> Taxonomic reference </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;width: 90px; font-weight: bold;"> 16S </td>
   <td style="text-align:center;width: 190px; "> 250 / 250 </td>
   <td style="text-align:left;width: 190px; "> SILVA 138.2 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 90px; font-weight: bold;"> 18S </td>
   <td style="text-align:center;width: 190px; "> 270 / 255 </td>
   <td style="text-align:left;width: 190px; "> PR2 5.1.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 90px; font-weight: bold;"> COI </td>
   <td style="text-align:center;width: 190px; "> 260 / 260 </td>
   <td style="text-align:left;width: 190px; "> MIDORI2 GB273 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 90px; font-weight: bold;"> MiFish </td>
   <td style="text-align:center;width: 190px; "> None </td>
   <td style="text-align:left;width: 190px; "> MitoFish (Jul 2026) </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 90px; font-weight: bold;"> ITS </td>
   <td style="text-align:center;width: 190px; "> 270 / 270 </td>
   <td style="text-align:left;width: 190px; "> UNITE (Feb 2025) </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 90px; font-weight: bold;"> rbcL </td>
   <td style="text-align:center;width: 190px; "> 250 / 250 </td>
   <td style="text-align:left;width: 190px; "> Diat.barcode 16.3 </td>
  </tr>
</tbody>
</table>

`````
:::
:::



**DADA2 settings shared across markers:**

- **Filtering:** `maxEE = c(2,3)`, `maxN = 0`, `minLen = 200`, `rm.phix = TRUE`.
- **Denoising:** Independent forward/reverse error models, with `pool = "pseudo"`.
- **Paired-end merging:** `minOverlap = 90`, `trimOverhang = TRUE`.
- **Chimera removal:** `method = "consensus"`, `allowOneOff = FALSE`.
- **Taxonomy:** `assignTaxonomy()` using marker-specific reference databases; 16S species-level assignment additionally used `addSpecies()`.

Truncation lengths are expressed as R1/R2 in base pairs. "None" indicates no fixed-length truncation. Other unspecified settings used DADA2 defaults.



::: {.cell}

```{.r .cell-code}
data_dir <- "data/Illumina"
markers <- c("16S", "18S", "COI", "MiFish", "ITS", "rbcL")

for (m in markers) {
  track <- readRDS(file.path(data_dir, paste0(m, "_read-count-tracking.rds")))
  seqtab <- readRDS(file.path(data_dir, paste0(m, "_seqtab_nochim.rds")))
  taxa   <- readRDS(file.path(data_dir, paste0(m, "_taxa.rds")))

  assign(paste0("track_", m), track)
  assign(paste0("seqtab_", m), seqtab)
  assign(paste0("taxa_", m), taxa)
}
```
:::



::: {.cell}

```{.r .cell-code}
qc_samples <- map_dfr(markers, function(m) {

  track <- get(paste0("track_", m))
  seqtab <- get(paste0("seqtab_", m))

  # Check sample matching
  stopifnot(setequal(rownames(track), rownames(seqtab)))

  # Match tracking rows to ASV table rows
  track <- track[rownames(seqtab), , drop = FALSE]

  tibble(
    marker = m,
    sample = rownames(seqtab),
    input_reads = track[, "dada2_input"],
    filtered_reads = track[, "filtered"],
    merged_reads = track[, "merged"],
    final_reads = track[, "nonchim"],
    retention = track[, "final_perc_reads_retained"],
    n_ASVs = rowSums(seqtab > 0),
    sample_type = case_when(
      str_detect(sample, "Extract-ctrl$") ~ "Extraction control",
      str_detect(sample, "NTC$") ~ "PCR control",
      TRUE ~ "Biological"
    )
  )
})
```
:::




::: {.cell}

```{.r .cell-code}
# 1. Combine DADA2 tracking data -----------------------------------------

qc_tracking <- purrr::map_dfr(markers, function(m) {

  get(paste0("track_", m)) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample") %>%
    mutate(
      marker = m,
      sample_type = case_when(
        str_detect(sample, "Extract-ctrl$") ~ "Extraction control",
        str_detect(sample, "NTC$") ~ "PCR control",
        TRUE ~ "Biological"
      )
    )
})

# 2. Calculate median reads at each processing step -----------------------

qc_steps <- qc_tracking %>%
  filter(sample_type == "Biological") %>%
  group_by(marker) %>%
  summarise(
    across(
      c(dada2_input, filtered, dada_f, merged, nonchim),
      median
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -marker,
    names_to = "step",
    values_to = "median_reads"
  ) %>%
  mutate(
    step = factor(
      step,
      levels = c("dada2_input", "filtered", "dada_f",
                 "merged", "nonchim")
    )
  )

# 3. Figure: DADA2 read retention ----------------------------------------

p_retention <- ggplot(
  qc_steps,
  aes(x = step, y = median_reads, group = 1)
) +
  geom_line(colour = "#0072B2", linewidth = 0.8) +
  geom_point(colour = "#0072B2", size = 2.5) +
  facet_wrap(~ marker, scales = "free_y") +
  scale_x_discrete(labels = c(
    dada2_input = "Input",
    filtered = "Filtered",
    dada_f = "Denoised F",
    merged = "Merged",
    nonchim = "Nonchim"
  )) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(
    x = NULL,
    y = "Median reads per sample",
    title = "DADA2 read retention across markers",
    subtitle = "Biological samples only"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 40, hjust = 1)
  )

print(p_retention)

# 4. Figure: Final sequencing depth per sample ----------------------------

p_depth <- qc_samples %>%
  mutate(
    marker = factor(marker, levels = markers),
    sample_type = factor(
      sample_type,
      levels = c("Biological", "Extraction control", "PCR control")
    ),
    sample_order = as.integer(str_extract(sample, "^[0-9]+"))
  ) %>%
  ggplot(aes(
    x = sample_order,
    y = final_reads,
    fill = sample_type
  )) +
  geom_col(width = 0.85) +
  facet_wrap(~ marker, scales = "free_y") +
  scale_fill_manual(values = c(
    "Biological" = "#0072B2",
    "Extraction control" = "#D55E00",
    "PCR control" = "#CC79A7"
  )) +
  scale_x_continuous(breaks = c(1, 10, 20, 30, 34, 35)) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(
    x = "Sample number",
    y = "Non-chimeric reads",
    fill = NULL,
    title = "Final sequencing depth per sample"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

print(p_depth)
```

::: {.cell-output-display}
![](Context_files/figure-html/fig-illumina-qc-1.png){#fig-illumina-qc-1 width=1152}
:::

::: {.cell-output-display}
![](Context_files/figure-html/fig-illumina-qc-2.png){#fig-illumina-qc-2 width=1152}
:::
:::




::: {.cell}

```{.r .cell-code}
# Calculate percentage retention per marker
qc_steps_perc <- qc_steps %>%
  group_by(marker) %>%
  mutate(
    retention = 100 * median_reads /
      median_reads[step == "dada2_input"]
  ) %>%
  ungroup()

# Combined plot

p_retention_combined <- ggplot(
  qc_steps_perc,
  aes(x = step, y = retention, colour = marker, group = marker)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_x_discrete(labels = c(
    dada2_input = "Input",
    filtered = "Filtered",
    dada_f = "Denoised F",
    merged = "Merged",
    nonchim = "Nonchim"
  )) +
  scale_y_continuous(
    breaks = seq(80, 100, 5),
    labels = scales::label_percent(scale = 1)
  ) +
  coord_cartesian(ylim = c(80, 100)) +
  scale_colour_brewer(palette = "Dark2") +
  labs(
    x = NULL,
    y = "Reads retained (%)",
    colour = "Marker",
    title = "DADA2 read retention across markers"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

p_retention_combined
```

::: {.cell-output-display}
![Median DADA2 read retention across markers, expressed relative to the median input read count.](Context_files/figure-html/fig-dada2-combined-1.png){#fig-dada2-combined width=864}
:::
:::



## Quality control by marker {.tabset}

::: {.panel-tabset}

### 16S


::: {.cell}
::: {.cell-output-display}
![](Context_files/figure-html/unnamed-chunk-4-1.png){width=864}
:::
:::


### 18S


::: {.cell}
::: {.cell-output-display}
![](Context_files/figure-html/unnamed-chunk-5-1.png){width=864}
:::
:::


### COI


::: {.cell}
::: {.cell-output-display}
![](Context_files/figure-html/unnamed-chunk-6-1.png){width=864}
:::
:::


### MiFish


::: {.cell}
::: {.cell-output-display}
![](Context_files/figure-html/unnamed-chunk-7-1.png){width=864}
:::
:::


### ITS


::: {.cell}
::: {.cell-output-display}
![](Context_files/figure-html/unnamed-chunk-8-1.png){width=864}
:::
:::


### rbcL


::: {.cell}
::: {.cell-output-display}
![](Context_files/figure-html/unnamed-chunk-9-1.png){width=864}
:::
:::


:::



::: {.cell}

```{.r .cell-code}
library(tidyverse)
library(vegan)

markers <- c("16S", "18S", "COI", "MiFish", "ITS", "rbcL")

# Calculate rarefaction curves -------------------------------------------

rarefaction <- purrr::map_dfr(markers, function(m) {

  seqtab <- get(paste0("seqtab_", m))

  # Exclude both negative controls
  seqtab <- seqtab[
    !grepl("Extract-ctrl$|NTC$", rownames(seqtab)),
    , drop = FALSE
  ]

  # Remove zero-depth samples
  seqtab <- seqtab[rowSums(seqtab) > 0, , drop = FALSE]

  # Depths covering the full range of each sample
  max_depth <- max(rowSums(seqtab))
  depths <- unique(round(seq(1, max_depth, length.out = 60)))

  purrr::map_dfr(seq_len(nrow(seqtab)), function(i) {

    sample_depth <- rowSums(seqtab[i, , drop = FALSE])

    d <- unique(c(
      depths[depths <= sample_depth],
      sample_depth
    ))

    tibble(
      marker = m,
      sample = rownames(seqtab)[i],
      depth = d,
      richness = vapply(
        d,
        function(n) as.numeric(
          vegan::rarefy(seqtab[i, , drop = FALSE], sample = n)
        ),
        numeric(1)
      )
    )
  })
})

# Summary curves:
# Recalculate on a common depth grid for each marker so that
# all biological samples contribute at every plotted depth.

rarefaction_summary <- purrr::map_dfr(markers, function(m) {

  x <- get(paste0("seqtab_", m))

  # Exclude negative controls
  x <- x[!grepl("Extract-ctrl$|NTC$", rownames(x)), , drop = FALSE]

  # Exclude COI sample 22
  if (m == "COI") {
    x <- x[!grepl("^22_COI_", rownames(x)), , drop = FALSE]
  }

  depths <- unique(round(seq(1, min(rowSums(x)), length.out = 60)))

  purrr::map_dfr(depths, function(d) {

    richness <- as.numeric(vegan::rarefy(x, sample = d))

    tibble(
      marker = m,
      depth = d,
      mean_richness = mean(richness),
      sd_richness = sd(richness)
    )
  })
})

# Combined figure --------------------------------------------------------

p_rarefaction_combined <- ggplot(
  rarefaction_summary,
  aes(x = depth, y = mean_richness, colour = marker, fill = marker)
) +
  geom_ribbon(
    aes(
      ymin = pmax(0, mean_richness - sd_richness),
      ymax = mean_richness + sd_richness
    ),
    alpha = 0.15,
    colour = NA
  ) +
  geom_line(linewidth = 1.1) +
  scale_colour_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_continuous(labels = scales::label_comma()) +
  labs(
    x = "Sequencing depth (reads)",
    y = "Expected ASV richness",
    colour = "Marker",
    fill = "Marker",
    title = "Rarefaction curves across markers",
    subtitle = "Mean ± SD across biological samples"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

p_rarefaction_combined
```

::: {.cell-output-display}
![](Context_files/figure-html/fig-rarefaction-combined-1.png){#fig-rarefaction-combined width=864}
:::
:::



::: {.cell}

```{.r .cell-code}
plot_rarefaction <- function(m) {

  rarefaction %>%
    filter(marker == m) %>%
    ggplot(aes(
      x = depth,
      y = richness,
      group = sample
    )) +
    geom_line(
      colour = "#0072B2",
      alpha = 0.45,
      linewidth = 0.6
    ) +
    scale_x_continuous(labels = scales::label_comma()) +
    labs(
      x = "Sequencing depth (reads)",
      y = "Expected ASV richness",
      title = paste0(m, " rarefaction curves"),
      subtitle = "Individual biological samples"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank()
    )
}
```
:::



::: {.panel-tabset}

## 16S


::: {.cell}
::: {.cell-output-display}
![](Context_files/figure-html/unnamed-chunk-10-1.png){width=864}
:::
:::


## 18S


::: {.cell}
::: {.cell-output-display}
![](Context_files/figure-html/unnamed-chunk-11-1.png){width=864}
:::
:::


## COI


::: {.cell}
::: {.cell-output-display}
![](Context_files/figure-html/unnamed-chunk-12-1.png){width=864}
:::
:::


## MiFish


::: {.cell}
::: {.cell-output-display}
![](Context_files/figure-html/unnamed-chunk-13-1.png){width=864}
:::
:::



## ITS


::: {.cell}
::: {.cell-output-display}
![](Context_files/figure-html/unnamed-chunk-14-1.png){width=864}
:::
:::


## rbcL


::: {.cell}
::: {.cell-output-display}
![](Context_files/figure-html/unnamed-chunk-15-1.png){width=864}
:::
:::


:::

# 📦 Building the phyloseq objects


::: {.cell}

```{.r .cell-code}
library(phyloseq)

markers <- c("16S", "18S", "COI", "MiFish", "ITS", "rbcL")

for (m in markers) {

  seqtab <- get(paste0("seqtab_", m))
  taxa <- get(paste0("taxa_", m))

  # Check that taxonomy and ASV sequences match
  stopifnot(setequal(colnames(seqtab), rownames(taxa)))

  # Align taxonomy to the ASV table
  taxa <- taxa[colnames(seqtab), , drop = FALSE]

  # Assign simple ASV identifiers
  asv_ids <- paste0("ASV", seq_len(ncol(seqtab)))

  sequences <- colnames(seqtab)

  colnames(seqtab) <- asv_ids
  rownames(taxa) <- asv_ids

  # Sample metadata
  metadata <- data.frame(
    SampleID = rownames(seqtab),
    SampleType = ifelse(
      grepl("Extract-ctrl$", rownames(seqtab)), "Extraction control",
      ifelse(
        grepl("NTC$", rownames(seqtab)), "PCR control", "Biological"
      )
    ),
    row.names = rownames(seqtab)
  )

  # Create phyloseq object
  ps <- phyloseq(
    otu_table(seqtab, taxa_are_rows = FALSE),
    tax_table(as.matrix(taxa)),
    sample_data(metadata),
    refseq(Biostrings::DNAStringSet(
      setNames(sequences, asv_ids)
    ))
  )

  assign(paste0("ps_", m), ps)
}

# 18S - PR2
colnames(tax_table(ps_18S)) <- c(
  "Domain", "Supergroup", "Division", "Subdivision",
  "Class", "Order", "Family", "Genus", "Species"
)

# rbcL - Diat.barcode 16.3
colnames(tax_table(ps_rbcL)) <- c(
  "Empire", "Kingdom", "Subkingdom", "Phylum",
  "Subphylum", "Class", "Order", "Family",
  "Genus", "Species"
)
```
:::


## Proportion of unassigned ASVs


::: {.cell}

```{.r .cell-code}
unclassified_summary <- purrr::map_dfr(markers, function(m) {

  ps <- get(paste0("ps_", m))

  tax_df <- as.data.frame(as(tax_table(ps), "matrix"))

  tax_df %>%
    summarise(across(
      everything(),
      ~ sum(
        is.na(.) |
        trimws(as.character(.)) == "" |
        tolower(trimws(as.character(.))) == "unclassified"
      )
    )) %>%
    pivot_longer(
      everything(),
      names_to = "Rank",
      values_to = "n_unclassified_ASVs"
    ) %>%
    mutate(
      Marker = m,
      total_ASVs = ntaxa(ps),
      percent_unclassified = 100 * n_unclassified_ASVs / total_ASVs
    )
})

unclassified_summary
```

::: {.cell-output .cell-output-stdout}

```
# A tibble: 46 × 5
   Rank       n_unclassified_ASVs Marker total_ASVs percent_unclassified
   <chr>                    <int> <chr>       <int>                <dbl>
 1 Kingdom                      6 16S          6697               0.0896
 2 Phylum                     252 16S          6697               3.76  
 3 Class                      546 16S          6697               8.15  
 4 Order                     1176 16S          6697              17.6   
 5 Family                    2357 16S          6697              35.2   
 6 Genus                     4549 16S          6697              67.9   
 7 Species                   6630 16S          6697              99.0   
 8 Domain                      80 18S         10697               0.748 
 9 Supergroup                 339 18S         10697               3.17  
10 Division                   404 18S         10697               3.78  
# ℹ 36 more rows
```


:::
:::



::: {.cell}

```{.r .cell-code}
# Preserve taxonomic rank order as stored in the taxonomy tables
rank_order <- unique(unclassified_summary$Rank)

unclassified_summary <- unclassified_summary %>%
  mutate(
    Marker = factor(Marker, levels = markers),
    Rank = factor(Rank, levels = rank_order)
  )

p_unclassified <- ggplot(
  unclassified_summary,
  aes(x = Rank, y = percent_unclassified)
) +
  geom_col(fill = "#0072B2", width = 0.75) +
  facet_wrap(~ Marker, scales = "free_x") +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    x = "Taxonomic rank",
    y = "Unclassified ASVs (%)",
    title = "Proportion of unclassified ASVs across markers"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text = element_text(face = "bold")
  )

p_unclassified
```

::: {.cell-output-display}
![](Context_files/figure-html/plot_unassigned_ASVs-1.png){width=672}
:::
:::



::: {.cell}

```{.r .cell-code}
# Preserve taxonomic rank order as stored in the taxonomy tables
rank_order <- unique(unclassified_summary$Rank)

unclassified_summary <- unclassified_summary %>%
  mutate(
    Marker = factor(Marker, levels = markers),
    Rank = factor(Rank, levels = rank_order)
  )

p_unclassified <- ggplot(
  unclassified_summary,
  aes(x = Rank, y = percent_unclassified)
) +
  geom_col(fill = "#0072B2", width = 0.75) +
  facet_wrap(~ Marker, scales = "free_x") +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    x = "Taxonomic rank",
    y = "Unclassified ASVs (%)",
    title = "Proportion of unclassified ASVs across markers"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text = element_text(face = "bold")
  )

p_unclassified
```

::: {.cell-output-display}
![](Context_files/figure-html/unnamed-chunk-16-1.png){width=1152}
:::
:::



::: {.cell}

```{.r .cell-code}
library(tidyverse)
library(phyloseq)

markers <- c("16S", "18S", "COI", "MiFish", "ITS", "rbcL")

# Define a general taxonomic order spanning all databases
rank_levels <- c(
  "Empire",
  "Domain",
  "Kingdom",
  "Supergroup",
  "Division",
  "Subdivision",
  "Subkingdom",
  "Phylum",
  "Subphylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Species"
)

# Count assigned ASVs at each rank for each marker
assigned_summary <- purrr::map_dfr(markers, function(m) {
  
  ps <- get(paste0("ps_", m))
  tax_df <- as.data.frame(as(tax_table(ps), "matrix"))
  
  tibble(
    Rank = colnames(tax_df),
    n_assigned_ASVs = sapply(tax_df, function(x) {
      sum(
        !is.na(x) &
        trimws(as.character(x)) != "" &
        !tolower(trimws(as.character(x))) %in% c("na", "unclassified")
      )
    }),
    Marker = m,
    total_ASVs = ntaxa(ps)
  )
})

# Keep only ranks that are actually present
assigned_summary <- assigned_summary %>%
  mutate(
    Rank = factor(Rank, levels = rank_levels),
    Marker = factor(Marker, levels = markers)
  ) %>%
  arrange(Marker, Rank)

# Combined plot
p_assigned_combined <- ggplot(
  assigned_summary,
  aes(x = Rank, y = n_assigned_ASVs, colour = Marker, group = Marker)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_colour_brewer(palette = "Dark2") +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(
    x = "Taxonomic rank",
    y = "Assigned ASVs (n)",
    colour = "Marker",
    title = "Number of assigned ASVs across taxonomic ranks"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

p_assigned_combined
```

::: {.cell-output-display}
![](Context_files/figure-html/fig-assigned-asvs-combined-1.png){#fig-assigned-asvs-combined width=960}
:::
:::




::: {.cell}

```{.r .cell-code}
library(microbiome)

my_colors <- c(
  "#FF6666", "#8DD3C7", "#377EB8", "#4DAF4A",
  "#984EA3", "#A6D854", "#FFFF33", "#FF7F00",
  "#F781BF", "#999999", "#66C2A5", "#FC8D62",
  "#8DA0CB", "#E78AC3", "#FFD92F", "#E5C494",
  "#B3B3B3", "#A65628", "#FFFFB3", "#BEBADA",
  "#FB8072", "#80B1D3", "#FDB462"
)

plot_class_composition <- function(ps, marker, top_n = 20) {

  # 1. Keep biological samples only
  ps <- subset_samples(ps, SampleType == "Biological")
  ps <- prune_samples(sample_sums(ps) > 0, ps)
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)

  # 2. Replace missing class assignments
  tax <- as(tax_table(ps), "matrix")

  tax[, "Class"][
    is.na(tax[, "Class"]) |
    trimws(tax[, "Class"]) == "" |
    tolower(trimws(tax[, "Class"])) %in% c("na", "unclassified")
  ] <- "Unassigned"

  tax_table(ps) <- tax_table(tax)

  # 3. Agglomerate to Class
  ps_class <- tax_glom(
    ps,
    taxrank = "Class",
    NArm = FALSE
  )

  # 4. Convert to relative abundance
  ps_rel <- microbiome::transform(
    ps_class,
    "compositional"
  )

  # 5. Select top classes by mean relative abundance
  top_classes <- microbiome::top_taxa(
    ps_rel,
    n = top_n
  )

  # 6. Melt and group remaining classes
  df <- psmelt(ps_rel) %>%
    mutate(
      Class = as.character(Class),
      Class_plot = if_else(
        OTU %in% top_classes | Class == "Unassigned",
        Class,
        "Other"
      ),
      Number = as.integer(
        stringr::str_extract(Sample, "^[0-9]+")
      )
    ) %>%
    group_by(Sample, Number, Class_plot) %>%
    summarise(
      Abundance = sum(Abundance),
      .groups = "drop"
    )

  # 7. Order classes by total relative abundance
  class_order <- df %>%
    filter(!Class_plot %in% c("Other", "Unassigned")) %>%
    group_by(Class_plot) %>%
    summarise(
      Total = sum(Abundance),
      .groups = "drop"
    ) %>%
    arrange(desc(Total)) %>%
    pull(Class_plot)

  class_order <- c(
    class_order,
    intersect("Unassigned", unique(df$Class_plot)),
    intersect("Other", unique(df$Class_plot))
  )

  df$Class_plot <- factor(
    df$Class_plot,
    levels = class_order
  )

  # 8. Assign colours
  plot_colors <- setNames(
    rep(my_colors, length.out = length(class_order)),
    class_order
  )

  if ("Unassigned" %in% class_order) {
    plot_colors["Unassigned"] <- "#666666"
  }

  if ("Other" %in% class_order) {
    plot_colors["Other"] <- "#C0C0C0"
  }

  # 9. Plot
  ggplot(
    df,
    aes(
      x = factor(Number, levels = 1:33),
      y = Abundance,
      fill = Class_plot
    )
  ) +
    geom_col(width = 0.85) +
    scale_fill_manual(
      values = plot_colors,
      drop = FALSE
    ) +
    scale_y_continuous(
      labels = scales::percent,
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      title = paste0(marker, " — Class-level composition"),
      x = "Sample number",
      y = "Relative abundance (%)",
      fill = "Class"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      legend.position = "bottom",
      legend.direction = "horizontal"
    ) +
    guides(
      fill = guide_legend(nrow = 5, byrow = TRUE)
    )
}
```
:::



# Class-level taxonomical composition

## 16S


::: {.cell}
::: {.cell-output-display}
![](Context_files/figure-html/unnamed-chunk-18-1.png){width=1152}
:::
:::


## 18S


::: {.cell}
::: {.cell-output-display}
![](Context_files/figure-html/unnamed-chunk-19-1.png){width=1152}
:::
:::


## COI


::: {.cell}
::: {.cell-output-display}
![](Context_files/figure-html/unnamed-chunk-20-1.png){width=1152}
:::
:::


## MiFish


::: {.cell}
::: {.cell-output-display}
![](Context_files/figure-html/unnamed-chunk-21-1.png){width=1152}
:::
:::


## ITS


::: {.cell}
::: {.cell-output-display}
![](Context_files/figure-html/unnamed-chunk-22-1.png){width=1152}
:::
:::


## rbcL


::: {.cell}
::: {.cell-output-display}
![](Context_files/figure-html/unnamed-chunk-23-1.png){width=1152}
:::
:::



