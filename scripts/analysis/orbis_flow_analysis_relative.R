library(tidyverse)
library(pheatmap)
library(corrplot)
library(RColorBrewer)
library(factoextra)
library(FactoMineR)
library(Polychrome)
library(conflicted)
library(bestNormalize)
library(mice)
library(skimr)
library(qs)
library(SoupX)
library(ICD10gm)
library(datawizard)
library(WRS2)
library(tidymodels)
library(finetune)
options(tidymodels.dark = TRUE)

source("ml_izkf_utils.R")
project <- "relative"

# color palette ------------------------------------------
phmap_colors <- colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100) #nice colors for pheatmap

#large sequential color palette
set.seed(123)
my_cols <- unname(createPalette(100, RColorBrewer::brewer.pal(8, "Set2")))
scales::show_col(my_cols)

# section read in final data for analysis ------------------------------------------
all_data_one_fil <- qs::qread("final_one_rel.qs")

csf_data <- all_data_one_fil$csf
blood_data <- all_data_one_fil$blood

all_data_one_complete <- qs::qread("final_one_rel_complete.qs")
csf_data_complete <- all_data_one_complete$csf
blood_data_complete <- all_data_one_complete$blood

all_data_norm_complete <- qs::qread("final_one_rel_norm_complete.qs")
csf_norm_complete <- all_data_norm_complete$csf
blood_norm_complete <- all_data_norm_complete$blood

combined_complete <- qread("final_one_rel_combined_complete.qs")

combined_norm_complete <- qs::qread("final_one_rel_combined_norm_complete.qs")
combined_umap_full <- qs::qread("final_one_rel_umap_combined.qs")

data_combined_multi <- qs::qread("final_multi_comb_rel.qs")

sum(all_data_one_fil$blood$event_count) #1313 million events
sum(all_data_one_fil$csf$event_count)# 233 million events
# combined over 1,5 billion events after filtering

# section read in processed but unfilter  ------------------------------------------
#all_data duplicate (measure date repeated) removed, but multiple measurements of one patients kept
#all_data_one only one measurement per patient kept
#patient_id - same for each patient, but blood/CSF
#sample_id - same for blood and CSF, distinct for each patient and measurement (patient_id + measure_date)

#all_data <- read_csv("orbis_flow_rel.csv")
all_data_one <- read_csv("orbis_flow_rel_one.csv")
all_data <- read_csv("orbis_flow_rel.csv")

subfolders <- file.path(
  "analysis",
  "relative",
  c("qc", "categories", "correlation", "feature", "heatmap", "umap", "abundance", "top", "models", "interval")
)
lapply(subfolders, dir.create, recursive = TRUE)

# section histogram ------------------------------------------
# visualize data
all_data_one_long <-
    bind_rows(csf_data, blood_data) |>
    select(tissue, granulos:HLA_DR_T, lymphos_basic:lactate, harvest_volume, event_count) |>
    pivot_longer(granulos:event_count, names_to = "variable", values_to = "value")

 ggplot(all_data_one_long, aes(x = value, fill = tissue)) +
#    geom_density(alpha = 0.3) +
    geom_histogram(data = dplyr::filter(all_data_one_long, tissue == "blood"), fill = "red", bins = 100, alpha = 0.2) +
    geom_histogram(data = dplyr::filter(all_data_one_long, tissue == "CSF"), fill = "blue", bins = 100, alpha = 0.2) +
    facet_wrap(vars(variable), scales = "free", ncol = 4) +
    theme_bw()

ggsave(file.path("analysis", "relative", "qc", "histogram.pdf"), width = 10, height = 20)

# section count categories ------------------------------------------
categories <- c("dx_icd_level1", "dx_icd_level2", "dx_biobanklist_level1", "dx_biobanklist_level2", "tx_biobanklist", "dx_andi_level1", "dx_andi_level2", "dx_andi_level3")
sel_categories <- c("dx_icd_level1", "dx_icd_level2")

lapply(sel_categories, count_category, data = combined_complete)

plot_category(data = combined_complete, category = "dx_icd_level1", width = 4, height = 2)
plot_category(data = combined_complete, category = "dx_icd_level2", width = 7, height = 7)

# age sex histograms ------------------------------------------
sex_age_histogram <-
  combined_data |>
  dplyr::filter(!is.na(sex)) |>
  ggplot(aes(x = age, fill = sex)) +
  #    geom_density(alpha = 0.3) +
  ## geom_histogram(data = dplyr::filter(csf_data, sex == "f"), fill = "red", bins = 100, alpha = 0.2) +
  ## geom_histogram(data = dplyr::filter(csf_data, sex == "m"), fill = "blue", bins = 100, alpha = 0.2) +
  geom_histogram(bins = 25) +
  facet_wrap(vars(sex), scales = "free", ncol = 4) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(plot = sex_age_histogram, file.path("analysis", "relative", "basic", "sex_age_histogram.pdf"), width = 7, height = 5)

median(csf_data$age)
dplyr::count(csf_data, sex)


# correlation plot ------------------------------------------
#remove all those with only missing NA
#rename those with two "CSF" in their name, like protein_CSF_CSF
# cor_data <-
#   bind_rows(csf_data, blood_data) |>
#   select(sample_pair_id, tissue, granulos:lactate) |>
#   pivot_wider(names_from = tissue, values_from = granulos:lactate) |>
#   select(where(function(x) !all(is.na(x)))) |>
#   select(-sample_pair_id) |>
#   rename_with(function(x) str_remove(x, "_CSF"), c(protein_CSF_CSF:IgM_ratio_CSF, glucose_CSF_CSF)) |>
#   cor(use = "complete.obs", method = "spearman")

cor_dat <-
  combined_complete |>
  select(granulos_CSF:lactate_CSF) |>
  cor(method = "spearman")

pdf(file.path("analysis", "relative", "correlation", "corplot_spearman.pdf"), width = 8, height = 8)
corrplot(cor_data, order = "hclust", method = "color", col = phmap_colors, tl.col = "black", cl.cex = 0.8, tl.cex = 0.5, hclust.method = "ward.D")
dev.off()

##correlation with age difficult because correlates strongly with diseases


## #compare with andi
## freq_dx_icd_andi <-
##   combined_norm_complete |>
##   dplyr::filter(dx_icd_level2 %in% top_10_dx_icd_level2) |>
##   drop_na(dx_andi_level2) |>
##   count(dx_icd_level2, dx_andi_level2)

## freq_dx_icd_andi_phmap <-
##   freq_dx_icd_andi |>
##   pivot_wider(names_from = dx_andi_level2, values_from = n, values_fill = 0) |>
##   dplyr::filter(!is.na(dx_icd_level2)) |>
##   column_to_rownames(var = "dx_icd_level2") |>
##   pheatmap()


# compare with biobank  ------------------------------------------
top_dx_icd_level2 <-
  combined_norm_complete |>
  count(dx_icd_level2) |>
  arrange(desc(n)) |>
  dplyr::filter(!is.na(dx_icd_level2)) |>
  slice(1:10) |>
  pull(dx_icd_level2)

top_dx_icd_level2_manual <- c("somatoform", "multiple sclerosis", "dementia", "ischemic stroke", "Parkinson’s syndrome", "opticus neuritis", "transient ischemic attack", "viral encephalitis", "bacterial meningitis")


dplyr::count(combined_norm_complete, dx_biobanklist_level2) |>
  arrange(desc(n)) |>
  write_csv(file.path("biobank", "biobank_dx.csv"))

## write_csv(tibble(dx = unique(combined_norm_complete$dx_biobanklist_level2)), "biobank/biobank_dx.csv")

biobank_lookup <- read_csv(file.path("biobank", "biobank_dx_lookup.csv")) |>
  dplyr::filter(group != "remove") |>
  select(-n)

read_csv(file.path("biobank", "biobank_dx_lookup.csv")) |>
  dplyr::filter(group != "remove") |>
  select(n) |>
  sum()

freq_dx_icd_biobank <-
  combined_norm_complete |>
  dplyr::filter(dx_icd_level2 %in% top_dx_icd_level2_manual) |>
  left_join(biobank_lookup, by = c("dx_biobanklist_level2")) |>
  count(dx_icd_level2, group) |>
  drop_na(group)

  freq_dx_icd_biobank |>
  dplyr::filter(dx_icd_level2 == "transient ischemic attack")

  15/(15+2+17)

#plot heatmap, arrange by specific order
freq_dx_icd_biobank_phmap <-
  freq_dx_icd_biobank |>
  arrange(match(group, top_dx_icd_level2_manual)) |>
  pivot_wider(names_from = group, values_from = n, values_fill = 0) |>
  arrange(match(dx_icd_level2, top_dx_icd_level2_manual)) |>
  column_to_rownames(var = "dx_icd_level2")   |>
  pheatmap(
    scale = "row",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    border_color = NA,
    color = phmap_colors,
    cellwidth = 10,
    cellheight = 10,
    )

grDevices::cairo_pdf(
  file.path("analysis", "relative", "categories", "compare_dx_biobank.pdf"),
  width = 10,
  height = 5
)
print(freq_dx_icd_biobank_phmap)
dev.off()




# CSF/blood ratios ------------------------------------------
# use combined data, only keep if blood and CSF are both present
# calculcate albumin ratio manually (otherwise imputed data might be used)
# only keep if both CSF and blood are non-zero
# join with meta data
combined_ratio <-
  bind_rows(csf_data_complete, blood_data_complete) |>
  select(sample_pair_id, tissue, granulos:lactate) |>
  pivot_wider(names_from = tissue, values_from = granulos:lactate) |>
  mutate(albumin_ratio = albumin_CSF_CSF/albumin_serum_CSF) |>
  select(-albumin_ratio_CSF) |>
  select(where(function(x) !all(is.na(x)))) |>
  drop_na() |>
  rename_with(function(x) str_remove(x, "_CSF"), c(protein_CSF_CSF:IgM_ratio_CSF, glucose_CSF_CSF)) |>
  rename_with(.fn = function(x) gsub(x = x, pattern = "serum", replacement = "blood"),
              .cols = matches("_serum$")) |>
  left_join(
    select(csf_data_complete, patient_id, sample_pair_id, dx_icd_level1:lp_interval),
    by = "sample_pair_id")

# sanity check
combined_ratio |>
  dplyr::filter(patient_id == "111301") |>
  select(sample_pair_id, monos_CSF, dx_icd_level2)

csf_data_complete |>
  dplyr::filter(patient_id == "111301") |>
  select(sample_pair_id, monos, dx_icd_level2)

# create vars for ratios (only blood because CSF has some like lymphos basic that do not match)
ratio_vars <- str_subset(colnames(combined_ratio), pattern = "_blood$") |>
  gsub(pattern = "_blood", replacement = "_ratio")

blood_vars <- gsub(x = ratio_vars, pattern = "_ratio", replacement = "_blood")
CSF_vars <- gsub(x = ratio_vars, pattern = "_ratio", replacement = "_CSF")

names_ratio <- paste0(ratio_vars, "_norm")

# create ratios of those vars
combined_ratio[names_ratio] <-
  map2_dfc(
    select(combined_ratio, all_of(CSF_vars)),
    select(combined_ratio, all_of(blood_vars)),
    function(x, y) (x/y)/combined_ratio$albumin_ratio)

#remove infinite values (because divided by 0)
combined_ratio <-
  combined_ratio |>
  dplyr::filter(if_all(.cols = all_of(names_ratio), is.finite))

# section heatmap grouped CSF  with ratios ------------------------------------------
#first normalize then mean
#better results when leaving out step_normalize, especially visuable in individual heatmap
colnames(combined_ratio)

phmap_comb_norm <-
  combined_ratio |>
  ## select(dx_icd_level1, granulos_CSF:HLA_DR_T_ratio) |>
  select(dx_icd_level1, granulos_CSF:lactate_CSF, granulos_ratio_norm:HLA_DR_T_ratio_norm) |>
  recipes::recipe(dx_icd_level1 ~ .) |>
  bestNormalize::step_orderNorm(recipes::all_numeric()) |>
  recipes::prep() |>
  recipes::bake(new_data = NULL)

#histograms
phmap_comb_norm |>
    select(-dx_icd_level1) |>
    pivot_longer(everything(), names_to = "variable", values_to = "value") |>
    ggplot(aes(x=value))+
    geom_histogram(bins = 50) +
    facet_wrap(vars(variable), scales = "free", ncol = 4)
ggsave(file.path("analysis", "relative", "qc", "histogram_comb_ratio_norm.pdf"), width = 10, height = 30)

phmap_comb_group_data <-
  phmap_comb_norm |>
    drop_na(dx_icd_level1) |>
    group_by(dx_icd_level1) |>
    dplyr::summarize(across(c(granulos_CSF:lactate_CSF, granulos_ratio_norm:HLA_DR_T_ratio_norm),
                            function(x) mean(x, na.rm = TRUE))) |>
    column_to_rownames(var = "dx_icd_level1")

phmap_comb_group_data |>
  pivot_longer(everything(), names_to = "variable", values_to = "value") |>
  ggplot(aes(x = value)) +
  geom_histogram(bins = 10) +
  facet_wrap(vars(variable), scales = "free", ncol = 4)
ggsave(file.path("analysis", "relative", "qc", "histogram_comb_norm_mean.pdf"), width = 10, height = 30)

vars_combo_plot <-
  combined_ratio |>
  select(granulos_CSF:lactate_CSF, granulos_ratio_norm:HLA_DR_T_ratio_norm) |>
  names()

debugonce(heatmap_group)

names(combined_ratio)

debugonce(heatmap_group)

heatmap_group(
  category = "dx_icd_level2",
  data = combined_ratio,
  vars = vars_combo_plot,
  label = "combo_ratio",
  cutree_rows = 10,
  height = 15,
  colors = phmap_colors
)

debugonce(heatmap_group)

vars_ratio_plot <-
  combined_ratio |>
  select(granulos_ratio_norm:HLA_DR_T_ratio_norm) |>
  names()

heatmap_group(
  category = "dx_icd_level2",
  data = combined_ratio,
  vars = vars_ratio_plot,
  label = "combo_ratio",
  cutree_rows = 10,
  height = 15,
  colors = phmap_colors
)

qs::qsave(combined_ratio, "combined_ratio.qs")

# experimental: extract icd metadata ------------------------------------------
csf_data |>
  select(patient_id, hd, hd_g) |>

str(csf_data$hd_g[1:10])
test3 <-   ICD10gm::get_icd_labels(csf_data$hd_g[1:10], year = 2022)
test3 <-   ICD10gm::icd_expand(icd_in = data.frame(ICD = csf_data$hd_g[1:10]), year = 2020)

## # import lookup from hpo github repo https://github.com/DiseaseOntology/HumanDiseaseOntology
## # works worse than from private github repo, see below, so don't use
## lookup_icd10_hp <-
##   read_tsv("./ICD10inDO.tsv") |>
##   dplyr::mutate(icd10 = gsub(x = xrefs, pattern = "(ICD10CM):([A-Z0-9])", replacement = "\\2"))

## test4 <-
##   csf_data |>
##   select(patient_id, hd_g) |>
##   mutate(hd_g_short = gsub(x = hd_g, pattern = "\\..*", replacement = "")) |>
##   left_join(lookup_icd10_hp, by = c("hd_g" = "icd10")) |>
##   left_join(lookup_icd10_hp, by = c("hd_g_short" = "icd10")) |>
##   mutate(label = coalesce(label.x, label.y)) |>
##   mutate(id = coalesce(id.x, id.y)) |>
##   mutate(xrefs = coalesce(xrefs.x, xrefs.y)) |>
##   select(-id.x, -id.y, -label.x, -label.y, -xrefs.x, -xrefs.y, -hd_g_short) |>
##   dplyr::distinct(patient_id, .keep_all = TRUE) |>
##   slice(90:100)

## skimr::skim(test4)

# import lookup mesh umsl mapping from https://github.com/kush02/Automated-ICD10-Codes-Assignment/blob/master/MESH_ICD10_Mapping.csv
lookup_icd10_mesh <-
  read_csv("./MESH_ICD10_Mapping.csv")

# join either by entire hd_g (i.e. g or star diagnosis of HD) or by first part of hd_g
# make distinct if multiple matches
csf_ontology_pre <-
  csf_data |>
  select(patient_id, hd_g) |>
  mutate(hd_g_short = gsub(x = hd_g, pattern = "\\..*", replacement = "")) |>
  left_join(lookup_icd10_mesh, by = c("hd_g" = "ICD10CM_CODE"))|>
  left_join(lookup_icd10_mesh, by = c("hd_g_short" = "ICD10CM_CODE")) |>
  mutate(mesh = coalesce(MESH_ID.x, MESH_ID.y)) |>
  mutate(umls = coalesce(UMLS_CUI.x, UMLS_CUI.y)) |>
  mutate(description = coalesce(DRESCP.x, DRESCP.y)) |>
  select(-MESH_ID.x, -MESH_ID.y, -UMLS_CUI.x, -UMLS_CUI.y, -DRESCP.x, -DRESCP.y)|>
  dplyr::distinct(patient_id, .keep_all = TRUE)

  slice(90:100)

skimr::skim(test5)

test5 |>
  dplyr::count(umls, hd_g) |>
  arrange(desc(n))

test5 |>
  dplyr::filter(umls == "C0026769") |>
  dplyr::count(hd_g, mesh)

lookup_mesh_do <-
  read_tsv("MeshinDO.tsv") |>
  dplyr::rename(mesh = xrefs,
                do = id) |>
  dplyr::mutate(mesh = gsub(x = mesh, pattern = "MESH:", replacement = ""))

csf_ontology <-
  csf_ontology_pre |>
  left_join(lookup_mesh_do, by = c("mesh" = "mesh")) |>
  dplyr::distinct(patient_id, .keep_all = TRUE)


names(csf_data)

r_present <-
  csf_data |>
  select(hd:nd_10) |>
  ##create a new column called R_present that is TRUE if any of the columns hd to nd_10 containts the letter R
  dplyr::mutate(R_present = if_any(hd:nd_10, function(x) str_detect(x, "R"))) |>
  print(width = Inf)

sum(r_present$R_present, na.rm = TRUE)/nrow(r_present)


all_data_one_fil$csf |>
  dplyr::filter(age < 1)

  data_combined_multi |>
  dplyr::filter(dx_icd_level2 == "somatoform") |>
  dplyr::group_by(patient_id) |>
    dplyr::filter(n() > 1) |>
    arrange(patient_id, interval) |>
   dplyr::count(patient_id)

# save database for kknms
all_data_one |>
  select(pid, fallnummer, patient_id, sample_pair_id, tissue, measure_date, birthdate, first_name_orbis, last_name_orbis, first_letter_lukas, second_letter_lukas) |>
  write_csv("orbis_flow_kknms.csv")

# save database for time extraction
all_data |>
  dplyr::select(pid, fallnummer, measure_date, aufnahme) |>
  write_csv("orbis_flow_time.csv")

combined_complete |>
  dplyr::select(pid, fallnummer, measure_date) |>
  write_csv("orbis_flow_time_selected.csv")

all_data |>
  dplyr::filter(tissue == "CSF") |>
  dplyr::filter(fallnummer == 50836061) |>
  dplyr::select(fallnummer, pid, measure_date, cell_count:lactate) |>
  write_csv("orbis_flow_example.csv")

all_data |>
  distinct(patient_id)

#add orbis id, measure dateand birthdate to data combined
lookup <- 
  read_csv("orbis_flow_rel.csv") |>
  dplyr::distinct(patient_id, .keep_all = TRUE)

data_combined_multi |>
left_join()