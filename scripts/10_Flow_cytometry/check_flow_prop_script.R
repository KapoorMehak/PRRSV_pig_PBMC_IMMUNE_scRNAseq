
## Flow cytometry for day 84 ###
.libPaths("rstudio/libs/4.4.1")

# 5 populations: GD, B, NK, Myeloid, AB

flow <- read.csv("./PRRSV/Flow_proportions.csv", sep= '')
sc_14 <- read.csv("./PRRSV/scrna_props_14dpi.csv")
sc_84 <- read.csv("./PRRSV/scrna_props_84dpi.csv")


flow14_5 <- flow %>%
  filter(dpi == 14, DataType == "PrcLv") %>%
  mutate(
    pop5 = case_when(
      CellType == "gdTCR"    ~ "GD",
      CellType == "AllBcells"~ "B",
      CellType == "NKcells"  ~ "NK",
      CellType == "CD172sp"  ~ "Myeloid", 
      CellType == "abTCR"    ~ "AB",
      TRUE ~ NA_character_
    ),
    prop_flow = Value / 100   
  ) %>%
  filter(!is.na(pop5)) %>%
  group_by(Sample, Treatment, pop5) %>%
  summarise(prop_flow = sum(prop_flow), .groups = "drop")
flow14_5 <- flow14_5 %>%
  mutate(
    Treatment = case_when(
      Treatment == "neg"     ~ "control",
      Treatment == "persist" ~ "persistent",
      Treatment == "extinct" ~ "extinct",
      TRUE ~ Treatment
    )
  )
valid_samples <- unique(sc_14$Sample)
flow14_5_filtered <- flow14_5 %>%
  filter(Sample %in% valid_samples)

flow84_5 <- flow %>%
  filter(dpi == 84, DataType == "PrcLv") %>%
  mutate(
    pop5 = case_when(
      CellType == "gdTCR"    ~ "GD",
      CellType == "AllBcells"~ "B",
      CellType == "NKcells"  ~ "NK",
      CellType == "CD172sp"  ~ "Myeloid", 
      CellType == "abTCR"    ~ "AB",
      TRUE ~ NA_character_
    ),
    prop_flow = Value / 100   
  ) %>%
  filter(!is.na(pop5)) %>%
  group_by(Sample, Treatment, pop5) %>%
  summarise(prop_flow = sum(prop_flow), .groups = "drop")
flow84_5 <- flow84_5 %>%
  mutate(
    Treatment = case_when(
      Treatment == "neg"     ~ "control",
      Treatment == "persist" ~ "persistent",
      Treatment == "extinct" ~ "extinct",
      TRUE ~ Treatment
    )
  )
sc14_5 <- sc_14 %>%
  mutate(
    pop5 = case_when(
      celltypes %in% c("CD2+ GD T cells" ,"CD2- GD T cells") ~ "GD",
      celltypes %in% c("B cells","ASCs") ~ "B",
      celltypes %in% c("NK cells") ~ "NK",
      celltypes %in% c("Monocytes","pDCs") ~ "Myeloid",
      TRUE ~ "AB"   
    )
  ) %>%
  filter(!is.na(pop5)) %>%  
  group_by(Sample, Treatment, pop5) %>%
  summarise(
    n = sum(n_cells, na.rm = TRUE),
    n_total = max(n_cells_total, na.rm = TRUE),
    prop_sc = n / n_total,
    .groups = "drop"
  )

sc84_5 <- sc_84 %>%
  mutate(
    pop5 = case_when(
      CellTypes %in% c("CD2+ GD T cells" ,"CD2- GD T cells") ~ "GD",
      CellTypes %in% c("B cells","ASC", "Transitional B-like") ~ "B",
      CellTypes %in% c("NK cells") ~ "NK",
      CellTypes %in% c("Monocytes","pDCs","cDCs") ~ "Myeloid",
      CellTypes %in% c("AB T cells") ~ "AB",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(pop5)) %>%          # <<< drop unmapped cell types
  group_by(Sample, Treatment, pop5) %>%
  summarise(
    n       = sum(n_cells, na.rm = TRUE),
    n_total = max(n_cells_total, na.rm = TRUE),
    prop_sc = n / n_total,
    .groups = "drop"
  )


#compare
comp14 <- inner_join(flow14_5, sc14_5,
                     by = c("Sample","Treatment","pop5"))
cor.test(comp14$prop_flow, comp14$prop_sc, method = 'spearman') 
cor_tests <- comp14 %>%
  group_by(pop5) %>%
  summarise(
    rho  = cor(prop_flow, prop_sc, method = "spearman", use = "complete.obs"),
    pval = cor.test(prop_flow, prop_sc, method = "spearman")$p.value,
    n    = n(),
    .groups = "drop"
  )

cor_tests
cor_by_trt <- comp14 %>%
  group_by(Treatment) %>%
  summarise(
    rho  = cor(prop_flow, prop_sc, method = "spearman", use = "complete.obs"),
    pval = cor.test(prop_flow, prop_sc, method = "spearman")$p.value,
    n    = n(),
    .groups = "drop"
  )

cor_by_trt

comp84 <- inner_join(flow84_5, sc84_5,
                     by = c("Sample","Treatment","pop5"))
cor.test(comp84$prop_flow, comp84$prop_sc, method = 'spearman') 
cor_tests_84 <- comp84 %>%
  group_by(pop5) %>%
  summarise(
    rho  = cor(prop_flow, prop_sc, method = "spearman", use = "complete.obs"),
    pval = cor.test(prop_flow, prop_sc, method = "spearman")$p.value,
    n    = n(),
    .groups = "drop"
  )

cor_tests_84
cor_by_trt_84 <- comp84 %>%
  group_by(Treatment) %>%
  summarise(
    rho  = cor(prop_flow, prop_sc, method = "spearman", use = "complete.obs"),
    pval = cor.test(prop_flow, prop_sc, method = "spearman")$p.value,
    n    = n(),
    .groups = "drop"
  )

cor_by_trt_84




###plotting ###
cor_tests_plot <- cor_tests%>%
  mutate(
    pop5 = factor(pop5, levels = c("AB", "B", "GD", "Myeloid", "NK")),
    sig = case_when(
      pval < 0.05  ~ "*"
    )
  )

ggplot(cor_tests_plot, aes(x = pop5, y = rho)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey70") +
  geom_bar(stat = "identity", aes(fill = pop5)) +
  geom_text(aes(label = sprintf("%.3f", rho)), vjust = -0.2, size = 3.5) +
  geom_text(aes(y = rho + 0.05, label = sig), vjust = 0, size = 4) +
  coord_cartesian(ylim = c(0, 1.1)) +
  theme_bw() +
  labs(
    title = "Correlation by celltype - 14 DPI",
    x = "Cellular Phenotype",
    y = "Spearman correlation"
  )+
  theme(
    axis.text.x = element_text( hjust = 1, face = 'bold'),
    axis.text.y = element_text(face = 'bold')
  )
##trt 
cor_trts_plot <- cor_by_trt%>%
  mutate(
    pop5 = factor(Treatment, levels = c("control", "extinct", "persistent")),
    sig = case_when(
      pval < 0.05  ~ "*"
    )
  )

ggplot(cor_trts_plot, aes(x = pop5, y = rho)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey70") +
  geom_bar(stat = "identity", aes(fill = Treatment)) +
  geom_text(aes(label = sprintf("%.3f", rho)), vjust = -0.2, size = 3.5) +
  geom_text(aes(y = rho + 0.05, label = sig), vjust = 0, size = 4) +
  coord_cartesian(ylim = c(0, 1.1)) +
  theme_bw() +
  labs(
    title = "Correlation by Treatment - 14 DPI",
    x = "Treatment",
    y = "Spearman correlation"
  )+
  theme(
    axis.text.x = element_text( hjust = 1, face = 'bold'),
    axis.text.y = element_text(face = 'bold')
  )



summ <- comp14 %>%
  group_by(Treatment, pop5) %>%
  summarise(
    flow_mean = mean(prop_flow, na.rm = TRUE),
    flow_sd   = sd(prop_flow,   na.rm = TRUE),
    sc_mean   = mean(prop_sc,   na.rm = TRUE),
    sc_sd     = sd(prop_sc,     na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(flow_mean, sc_mean, flow_sd, sc_sd),
    names_to = c("Method", ".value"),
    names_pattern = "(flow|sc)_(mean|sd)"
  )

ggplot(summ, aes(x = pop5, y = mean, color = Method)) +
  geom_point(
    position = position_dodge(width = 0.4),
    size = 3
  ) +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd, group = Method),
    position = position_dodge(width = 0.4),
    width = 0.1
  ) +
  facet_wrap(~ Treatment) +
  theme_bw() +
  labs(
    title = "Trt level: 14 DPI",
    x = "Cellular Phenotypes",
    y = "Proportion",
    color = "Method"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = 'bold'),
    axis.text.y = element_text(face = 'bold')
  )


cor_txct <- comp84 %>%
  group_by(Treatment, pop5) %>%
  summarise(
    rho   = cor(prop_flow, prop_sc, method = "spearman", use = "complete.obs"),
    pval  = cor.test(prop_flow, prop_sc, method = "spearman")$p.value,
    MAE   = mean(abs(prop_sc - prop_flow), na.rm = TRUE),
    Ratio = mean(prop_sc / prop_flow,    na.rm = TRUE),
    n     = n(),
    .groups = "drop"
  )
cor_plot <- cor_txct %>%
  mutate(
    pop5 = factor(pop5, levels = c("AB", "B", "GD", "Myeloid", "NK")),
    MAE_cat = case_when(
      MAE < 0.05 ~ "Close (<0.05)",
      MAE < 0.10 ~ "Okay (0.05-0.1)",
      TRUE       ~ "poor (>0.1)"
    ),
    MAE_cat = factor(MAE_cat, levels = c("Close (<0.05)","Okay (0.05-0.1)","poor (>0.1)")),
    p_star = case_when(
      pval < 0.05  ~ "*",
      TRUE         ~ ""
    ),
    ratio_lab = sprintf("R=%.2f", Ratio)
  )

ggplot(cor_plot, aes(x = Treatment, y = pop5)) +
  geom_tile(fill = "grey95", color = "white") +
  geom_point(
    aes(fill = rho, shape = MAE_cat),
    size   = 15,   
    color  = "black",
    stroke = 0.5
  ) +

  geom_text(
    aes(label = p_star),
    color = "white",
    fontface = "bold",
    size = 5,
    vjust = -0.5
  ) +

  geom_text(
    aes(label = ratio_lab),
    color = "white",
    size = 3,
    fontface = "bold",
    vjust = 1.8   
  ) +
  
  scale_fill_gradient2(
    low  = "steelblue",
    mid  = "white",
    high = "firebrick",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman ρ"
  ) +

  scale_shape_manual(
    values = c(
      "Close (<0.05)"    = 21,   # big circle
      "Okay (0.05-0.1)" = 22,   # square
      "poor (>0.1)"   = 24    # triangle
    ),
    name = "MAE level"
  ) +
  
  theme_bw() +
  labs(
    title    = "scRNA vs Flow - Day 84",
    subtitle = "Fill = correlation, Shape = MAE, Star = p-value, Ration = sc/flow ",
    x        = "Treatment",
    y        = "Cell type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid  = element_blank()
  )

comp14$n_total <- NULL; comp14$n <- NULL
comp84$n_total <- NULL; comp84$n <- NULL
