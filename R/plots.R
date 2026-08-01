# ============================================================================
# Figuras del manuscrito — funciones puras de ploteo (Chat 6, Trabajo A)
# ----------------------------------------------------------------------------
# Extraidas de los chunks de figura del .Rmd
# (analysis/integrated_preterm_prediction_workflow.Rmd). Cada funcion RECIBE por
# argumento los objetos YA CALCULADOS (resultados de CV, modelos entrenados,
# metadata, coordenadas de PCA, tablas de frecuencia) y DEVUELVE el objeto de
# figura (ggplot / patchwork / kable). NO re-entrenan modelos, NO recalculan
# metricas de CV, NO leen variables globales y NO imprimen (el .Rmd conserva la
# narracion cat()/print() y toda la preparacion de datos; aqui solo vive la
# construccion del grafico, movida TAL CUAL para preservar identidad byte-a-byte).
#
# Invariante (playbook §7): figure5 usa extract_fit_parsnip(best_rf$final_model)
# en el .Rmd (modelo YA entrenado); estas funciones solo grafican su importancia.
#
# Depends on the plotting/analysis packages being attached (ggplot2, dplyr, tidyr,
# stringr, patchwork, ggrepel, knitr, kableExtra), exactly as in the source .Rmd.
# ============================================================================


# ==========================================================================
# FIGURA 1 — Cohorte y diseño de muestreo longitudinal
# ==========================================================================

#' Violin de edad gestacional al parto por desenlace (Study Population)
#'
#' @param subject_outcomes Tabla a nivel sujeto con `preterm_status` y `sdg_parto`.
#' @param n_subjects,n_preterm,prevalence,n_term Escalares para el subtitulo.
#' @return Un objeto ggplot.
#' @export
build_population_ga_violin <- function(subject_outcomes, n_subjects, n_preterm,
                                       prevalence, n_term) {
  ggplot(subject_outcomes, aes(x = preterm_status, y = sdg_parto, fill = preterm_status)) +
    geom_violin(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 3) +
    geom_hline(yintercept = 37, linetype = "dashed", color = "red", size = 1.2) +
    scale_fill_manual(values = c("Term" = "#619CFF", "Preterm" = "#F8766D")) +
    labs(
      title = "Gestational Age at Delivery by Birth Outcome",
      subtitle = sprintf("n = %d subjects (Preterm: %d [%.1f%%], Term: %d [%.1f%%])",
                         n_subjects, n_preterm, prevalence, n_term, 100-prevalence),
      x = "Birth Outcome",
      y = "Gestational Age at Delivery (weeks)",
      caption = "Dashed red line indicates 37 weeks threshold for preterm birth definition"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(face = "bold")
    )
}

#' Panel A de la Figura 1: flowchart del diseño del estudio
#'
#' @param n_total_subjects,n_preterm,prevalence_subject,n_term Conteos a nivel sujeto.
#' @param n_total_samples,n_preterm_samples,prevalence_sample,n_term_samples Conteos a nivel muestra.
#' @param mean_samples,range_samples Resumen de muestras por sujeto.
#' @param outcome_colors Vector nombrado de colores ("Term"/"Preterm").
#' @return Un objeto ggplot (theme_void).
#' @export
build_figure1_flowchart <- function(n_total_subjects, n_preterm, prevalence_subject,
                                    n_term, n_total_samples, n_preterm_samples,
                                    prevalence_sample, n_term_samples,
                                    mean_samples, range_samples, outcome_colors) {
  ggplot() +
    # Set up coordinate system
    xlim(0, 10) + ylim(0, 10) +

    # Box 1: Total enrollment
    annotate("rect", xmin = 2, xmax = 8, ymin = 8.5, ymax = 9.5,
             fill = "white", color = "black", size = 1) +
    annotate("text", x = 5, y = 9,
             label = sprintf("Total Enrollment\nn = %d pregnant women", n_total_subjects),
             size = 3.5, fontface = "bold") +

    # Arrow 1
    annotate("segment", x = 5, xend = 5, y = 8.5, yend = 7.5,
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             size = 0.8) +

    # Box 2: Study design
    annotate("rect", xmin = 1.5, xmax = 8.5, ymin = 6.5, ymax = 7.5,
             fill = "#E8F4F8", color = "black", size = 1) +
    annotate("text", x = 5, y = 7,
             label = "Nested Case-Control Design\nLongitudinal Vaginal Microbiome Sampling",
             size = 3.2, fontface = "italic") +

    # Arrow 2 (split into two)
    annotate("segment", x = 5, xend = 5, y = 6.5, yend = 6,
             size = 0.8) +
    annotate("segment", x = 5, xend = 2.5, y = 6, yend = 6,
             size = 0.8) +
    annotate("segment", x = 5, xend = 7.5, y = 6, yend = 6,
             size = 0.8) +
    annotate("segment", x = 2.5, xend = 2.5, y = 6, yend = 5.3,
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             size = 0.8) +
    annotate("segment", x = 7.5, xend = 7.5, y = 6, yend = 5.3,
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             size = 0.8) +

    # Box 3a: Preterm births (left)
    annotate("rect", xmin = 0.5, xmax = 4.5, ymin = 3.8, ymax = 5.3,
             fill = alpha(outcome_colors["Preterm"], 0.2), color = outcome_colors["Preterm"], size = 1.2) +
    annotate("text", x = 2.5, y = 5,
             label = "Preterm Birth",
             size = 4, fontface = "bold", color = outcome_colors["Preterm"]) +
    annotate("text", x = 2.5, y = 4.55,
             label = sprintf("n = %d subjects (%.1f%%)", n_preterm, prevalence_subject),
             size = 3.2) +
    annotate("text", x = 2.5, y = 4.2,
             label = "< 37 weeks gestation",
             size = 2.8, fontface = "italic") +

    # Box 3b: Term births (right)
    annotate("rect", xmin = 5.5, xmax = 9.5, ymin = 3.8, ymax = 5.3,
             fill = alpha(outcome_colors["Term"], 0.2), color = outcome_colors["Term"], size = 1.2) +
    annotate("text", x = 7.5, y = 5,
             label = "Term Birth",
             size = 4, fontface = "bold", color = outcome_colors["Term"]) +
    annotate("text", x = 7.5, y = 4.55,
             label = sprintf("n = %d subjects (%.1f%%)", n_term, 100 - prevalence_subject),
             size = 3.2) +
    annotate("text", x = 7.5, y = 4.2,
             label = "≥ 37 weeks gestation",
             size = 2.8, fontface = "italic") +

    # Arrow 3 (convergence)
    annotate("segment", x = 2.5, xend = 2.5, y = 3.8, yend = 3.3,
             size = 0.8) +
    annotate("segment", x = 7.5, xend = 7.5, y = 3.8, yend = 3.3,
             size = 0.8) +
    annotate("segment", x = 2.5, xend = 5, y = 3.3, yend = 3.3,
             size = 0.8) +
    annotate("segment", x = 7.5, xend = 5, y = 3.3, yend = 3.3,
             size = 0.8) +
    annotate("segment", x = 5, xend = 5, y = 3.3, yend = 2.5,
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             size = 0.8) +

    # Box 4: Total samples
    annotate("rect", xmin = 1.5, xmax = 8.5, ymin = 1.3, ymax = 2.5,
             fill = "#FFF9E6", color = "black", size = 1) +
    annotate("text", x = 5, y = 2.2,
             label = sprintf("Total Samples: n = %d", n_total_samples),
             size = 3.8, fontface = "bold") +
    annotate("text", x = 5, y = 1.85,
             label = sprintf("Preterm: %d samples (%.1f%%) | Term: %d samples (%.1f%%)",
                            n_preterm_samples, prevalence_sample,
                            n_term_samples, 100 - prevalence_sample),
             size = 3) +
    annotate("text", x = 5, y = 1.55,
             label = sprintf("Mean %.1f samples/subject (range %s)",
                            mean_samples, range_samples),
             size = 2.8, fontface = "italic") +

    # Box 5: Note about pregnancy loss
    annotate("rect", xmin = 0.3, xmax = 4.2, ymin = 0.3, ymax = 0.9,
             fill = "white", color = "gray40", size = 0.5, linetype = "dashed") +
    annotate("text", x = 2.25, y = 0.6,
             label = "Note: Preterm births include 1 pregnancy\nloss at 20 weeks analyzed as PTB",
             size = 2.3, color = "gray20", hjust = 0.5) +

    # Clean theme
    theme_void() +
    theme(
      plot.margin = margin(10, 10, 10, 10)
    )
}

#' Panel B de la Figura 1: distribución de muestras por trimestre
#'
#' @param trimester_summary Tabla con `trimester`, `outcome`, `n`, `percentage`.
#' @param outcome_colors Vector nombrado de colores ("Term"/"Preterm").
#' @return Un objeto ggplot.
#' @export
build_figure1_trimester <- function(trimester_summary, outcome_colors) {
  ggplot(trimester_summary, aes(x = trimester, y = n, fill = outcome)) +

    # Bars
    geom_col(position = position_dodge(width = 0.8), width = 0.7,
             color = "black", size = 0.5) +

    # Sample counts on bars
    geom_text(aes(label = n),
              position = position_dodge(width = 0.8),
              vjust = -0.5, size = 3.5, fontface = "bold") +

    # Percentage labels inside bars
    geom_text(aes(label = sprintf("(%.1f%%)", percentage)),
              position = position_dodge(width = 0.8),
              vjust = 1.7, size = 3, color = "white") +

    # Color scale
    scale_fill_manual(values = outcome_colors, name = "Birth Outcome") +

    # Y-axis
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.15)),
      breaks = seq(0, 50, by = 10)
    ) +

    # Labels
    labs(
      x = "Trimester at Sample Collection",
      y = "Number of Samples",
      title = NULL
    ) +

    # Theme
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "top",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      axis.title = element_text(face = "bold", size = 11),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.subtitle = element_text(size = 9.5, hjust = 0.5,
                                    margin = margin(b = 10)),
      plot.margin = margin(10, 10, 10, 10)
    )
}

#' Panel C de la Figura 1: swimmer plot de trayectorias longitudinales
#'
#' @param swimmer_data Tabla a nivel sujeto (trayectorias y desenlace).
#' @param sample_points Puntos de muestra individuales.
#' @param n_preterm,n_total_subjects,n_term Conteos para bandas y etiquetas.
#' @param outcome_colors Vector nombrado de colores ("Term"/"Preterm").
#' @return Un objeto ggplot.
#' @export
build_figure1_swimmer <- function(swimmer_data, sample_points, n_preterm,
                                  n_total_subjects, n_term, outcome_colors) {
  ggplot() +

    # Background shading for PTB/Term groups
    annotate("rect",
             xmin = -Inf, xmax = Inf,
             ymin = 0.5, ymax = n_preterm + 0.5,
             fill = alpha(outcome_colors["Preterm"], 0.05)) +
    annotate("rect",
             xmin = -Inf, xmax = Inf,
             ymin = n_preterm + 0.5, ymax = n_total_subjects + 0.5,
             fill = alpha(outcome_colors["Term"], 0.05)) +

    # Vertical line at 37 weeks (PTB threshold)
    geom_vline(xintercept = 37, linetype = "dashed", color = "gray30", size = 0.8) +
    annotate("text", x = 37, y = n_total_subjects + 1,
             label = "37 weeks", size = 3, color = "gray30", hjust = -0.1) +

    # Horizontal lines for each subject (from first sample to delivery)
    geom_segment(data = swimmer_data,
                 aes(x = first_sample_ga, xend = ga_delivery,
                     y = subject_number, yend = subject_number,
                     color = outcome),
                 size = 1.2, alpha = 0.7) +

    # Points for individual samples
    geom_point(data = sample_points,
               aes(x = sdg_visita, y = subject_number, fill = outcome),
               shape = 21, size = 2.5, color = "white", stroke = 0.5, alpha = 0.8) +

    # End point markers (delivery)
    geom_point(data = swimmer_data,
               aes(x = ga_delivery, y = subject_number, fill = outcome),
               shape = 23, size = 3.5, color = "black", stroke = 0.8) +

    # Color scales
    scale_color_manual(values = outcome_colors, name = "Birth Outcome") +
    scale_fill_manual(values = outcome_colors, name = "Birth Outcome") +

    # Axes and labels
    scale_x_continuous(
      breaks = seq(10, 42, by = 4),
      limits = c(8, 43),
      expand = c(0.01, 0)
    ) +
    scale_y_continuous(
      breaks = seq(5, n_total_subjects, by = 5),
      expand = c(0.01, 0)
    ) +

    labs(
      x = "Gestational Age (weeks)",
      y = "Subject Number",
      title = NULL
    ) +

    # Group labels on y-axis
    annotate("text", x = 7, y = n_preterm/2 + 0.5,
             label = sprintf("Preterm\n(n = %d)", n_preterm),
             size = 3.5, fontface = "bold", color = outcome_colors["Preterm"],
             hjust = 0.5) +
    annotate("text", x = 7, y = n_preterm + (n_term/2) + 0.5,
             label = sprintf("Term\n(n = %d)", n_term),
             size = 3.5, fontface = "bold", color = outcome_colors["Term"],
             hjust = 0.5) +

    # Theme
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.title = element_text(face = "bold", size = 11),
      axis.text = element_text(size = 9),
      plot.margin = margin(10, 15, 10, 10)
    )
}

#' Figura 1 combinada (patchwork de los 3 paneles)
#'
#' @param panel_a Flowchart (build_figure1_flowchart()).
#' @param panel_b Barras por trimestre (build_figure1_trimester()).
#' @param panel_c Swimmer plot (build_figure1_swimmer()).
#' @return Un objeto patchwork.
#' @export
build_figure1_combined <- function(panel_a, panel_b, panel_c) {
  figure1_combined <- ((panel_a / panel_b) + plot_layout(heights = c(3, 2))) | panel_c +
    plot_annotation(
      tag_levels = 'A',
      tag_prefix = '',
      tag_suffix = '',
      theme = theme(
        plot.tag = element_text(size = 16, face = "bold"),
        plot.tag.position = c(0, 1)  # Top-left corner
      )
    ) +
    plot_layout(
      widths = c(1.2, 1)  # Adjust relative widths if needed
    )

  # Add main title
  figure1_combined <- figure1_combined +
    plot_annotation(
      title = "Figure 1. Study Cohort and Longitudinal Sampling Design",
      theme = theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0,
                                  margin = margin(b = 10)),
        plot.margin = margin(10, 10, 10, 10)
      )
    )

  figure1_combined
}


# ==========================================================================
# FIGURA 2 — Composición del microbioma (PCA de CLR)
# ==========================================================================

#' Panel A de la Figura 2: scatter de PCA por desenlace
#'
#' @param pca_data Coordenadas de PCA por muestra + `outcome`.
#' @param variance_explained Vector de % de varianza por PC.
#' @param outcome_colors Vector nombrado de colores ("Term"/"Preterm").
#' @return Un objeto ggplot.
#' @export
build_figure2_pca_scatter <- function(pca_data, variance_explained, outcome_colors) {
  ggplot(pca_data, aes(x = PC1, y = PC2, color = outcome, fill = outcome)) +

    # Confidence ellipses (95% level)
    stat_ellipse(geom = "polygon", alpha = 0.15,
                 level = 0.95, type = "norm", show.legend = FALSE) +

    # Individual points
    geom_point(size = 3, alpha = 0.7, stroke = 0.5, shape = 21, color = "white") +

    # Color scales
    scale_fill_manual(values = outcome_colors, name = "Birth Outcome") +
    scale_color_manual(values = outcome_colors, name = "Birth Outcome") +

    # Axes labels with variance explained
    labs(
      x = sprintf("PC1 (%.1f%% variance explained)", variance_explained[1]),
      y = sprintf("PC2 (%.1f%% variance explained)", variance_explained[2]),
      title = NULL
    ) +

    # Theme
    theme_minimal(base_size = 12) +
    theme(
      legend.position = c(0.85, 0.15),
      legend.background = element_rect(fill = "white", color = "gray80", size = 0.5),
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 10),
      axis.title = element_text(face = "bold", size = 11),
      axis.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray60", fill = NA, size = 0.8),
      plot.margin = margin(10, 10, 10, 10)
    )
}

#' Panel B de la Figura 2: scree plot (varianza explicada + acumulada)
#'
#' @param variance_explained Vector de % de varianza por PC.
#' @param cumulative_variance Vector de % de varianza acumulada.
#' @return Un objeto ggplot.
#' @export
build_figure2_scree <- function(variance_explained, cumulative_variance) {
  variance_df <- tibble(
    PC = factor(1:min(20, length(variance_explained)),
               levels = 1:min(20, length(variance_explained))),
    variance = variance_explained[1:min(20, length(variance_explained))],
    cumulative = cumulative_variance[1:min(20, length(variance_explained))]
  )

  ggplot(variance_df, aes(x = PC)) +

    # Individual variance (bars)
    geom_col(aes(y = variance), fill = "steelblue", alpha = 0.7, width = 0.7) +

    # Cumulative variance (line)
    geom_line(aes(y = cumulative, group = 1), color = "darkred", size = 1.2) +
    geom_point(aes(y = cumulative), color = "darkred", size = 2.5) +

    # Horizontal line at 80% and 90%
    geom_hline(yintercept = 80, linetype = "dashed", color = "gray40", size = 0.6) +
    geom_hline(yintercept = 90, linetype = "dashed", color = "gray40", size = 0.6) +

    # Annotations
    annotate("text", x = 18, y = 82, label = "80%",
             size = 3, color = "gray40", hjust = 0) +
    annotate("text", x = 18, y = 92, label = "90%",
             size = 3, color = "gray40", hjust = 0) +

    # Highlight PC1-PC2
    annotate("rect", xmin = 0.5, xmax = 2.5, ymin = 0, ymax = 100,
             alpha = 0.1, fill = "gold") +
    annotate("text", x = 1.5, y = 95,
             label = sprintf("PC1+PC2\n%.1f%%", cumulative_variance[2]),
             size = 3.5, fontface = "bold") +

    # Scales
    scale_y_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, by = 20),
      expand = c(0, 0)
    ) +
    scale_x_discrete(breaks = c(1, 5, 10)) +

    # Labels
    labs(
      x = "Principal Component",
      y = "Variance Explained (%)",
      title = NULL,
      subtitle = "Blue bars: individual variance | Red line: cumulative variance"
    ) +

    # Theme
    theme_minimal(base_size = 11) +
    theme(
      axis.title = element_text(face = "bold", size = 11),
      axis.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.subtitle = element_text(size = 9, hjust = 0.5, margin = margin(b = 10)),
      plot.margin = margin(10, 10, 10, 10)
    )
}

#' Panel C de la Figura 2: biplot con loadings de género
#'
#' @param pca_data Coordenadas de PCA por muestra + `outcome`.
#' @param top_loadings Loadings de los géneros que más contribuyen a PC1/PC2.
#' @param variance_explained Vector de % de varianza por PC.
#' @param outcome_colors Vector nombrado de colores ("Term"/"Preterm").
#' @return Un objeto ggplot.
#' @export
build_figure2_biplot <- function(pca_data, top_loadings, variance_explained, outcome_colors) {
  # Scale loadings for visualization (arrows)
  loading_scale <- 0.69  # Adjust if arrows are too long/short

  biplot_data <- top_loadings %>%
    mutate(
      PC1_arrow = PC1 * loading_scale,
      PC2_arrow = PC2 * loading_scale,
      # Clean genus names for display (remove prefixes if present)
      genus_clean = str_remove(genus, "^[a-z]__")
    )

  ggplot() +

    # Ellipses
    stat_ellipse(data = pca_data,
                 aes(x = PC1, y = PC2, color = outcome, fill = outcome),
                 level = 0.95, type = "norm", alpha = 0.15,
                 geom = "polygon", show.legend = FALSE) +

    # Sample points
    geom_point(data = pca_data,
               aes(x = PC1, y = PC2, fill = outcome),
               size = 3, alpha = 0.7, stroke = 0.5, shape = 21, color = "white") +

    # Loading arrows
    geom_segment(data = biplot_data,
                 aes(x = 0, y = 0, xend = PC1_arrow, yend = PC2_arrow),
                 arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
                 color = "grey20", size = 0.8, alpha = 0.7) +

    # Genus labels
    geom_text_repel(data = biplot_data,
                    aes(x = PC1_arrow, y = PC2_arrow, label = genus_clean),
                    size = 3, color = "grey20", fontface = "bold",
                    box.padding = 0.5, point.padding = 0.3,
                    segment.color = "gray50", segment.size = 0.3,
                    max.overlaps = 20) +

    # Color scales
    scale_fill_manual(values = outcome_colors, name = "Birth Outcome") +
    scale_color_manual(values = outcome_colors, name = "Birth Outcome") +

    # Axes
    labs(
      x = sprintf("PC1 (%.1f%% variance explained)", variance_explained[1]),
      y = sprintf("PC2 (%.1f%% variance explained)", variance_explained[2]),
      title = NULL,
      caption = sprintf("Arrows represent top %d genera contributing to PC1 and PC2 variation",
                       nrow(biplot_data))
    ) +

    # Theme
    theme_minimal(base_size = 12) +
    theme(
      legend.position = c(0.85, 0.15),
      legend.background = element_rect(fill = "white", color = "gray80", size = 0.5),
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 10),
      axis.title = element_text(face = "bold", size = 11),
      axis.text = element_text(size = 10),
      plot.caption = element_text(size = 9, hjust = 0.5, margin = margin(t = 10)),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray60", fill = NA, size = 0.8),
      plot.margin = margin(10, 10, 10, 10)
    )
}


# ==========================================================================
# Approach 3 — Visualización de selección de variables (3 paneles)
# ==========================================================================

#' Frecuencia de selección de variables por fold (Approach 3)
#'
#' @param var_frequency_approach3 Tabla con `variable`, `n_folds_selected`.
#' @return Un objeto ggplot.
#' @export
build_approach3_selection_freq <- function(var_frequency_approach3) {
  var_frequency_approach3 %>%
    mutate(
      selection_category = case_when(
        n_folds_selected == 5 ~ "All 5 folds (100%)",
        n_folds_selected >= 4 ~ "4-5 folds (80-100%)",
        n_folds_selected >= 3 ~ "3 folds (60%)",
        n_folds_selected >= 2 ~ "2 folds (40%)",
        TRUE ~ "1 fold (20%)"
      ),
      selection_category = factor(
        selection_category,
        levels = c("All 5 folds (100%)", "4-5 folds (80-100%)",
                   "3 folds (60%)", "2 folds (40%)", "1 fold (20%)")
      )
    ) %>%
    ggplot(aes(x = reorder(variable, n_folds_selected),
               y = n_folds_selected)) +
    geom_col(aes(fill = selection_category), alpha = 0.8) +
    geom_hline(yintercept = 5, linetype = "dashed",
               color = "red", size = 1) +
    geom_hline(yintercept = 4, linetype = "dashed",
               color = "orange", size = 0.5, alpha = 0.5) +
    scale_fill_manual(
      name = "Selection Stability",
      values = c(
        "All 5 folds (100%)" = "#1B5E20",
        "4-5 folds (80-100%)" = "#388E3C",
        "3 folds (60%)" = "#FBC02D",
        "2 folds (40%)" = "#F57C00",
        "1 fold (20%)" = "#D32F2F"
      )
    ) +
    coord_flip() +
    labs(
      title = "Variable Selection Frequency Across CV Folds (Approach 3)",
      subtitle = "Data-driven univariate screening performed independently in each fold",
      x = "Clinical Variable",
      y = "Number of Folds Selected (out of 5)",
      caption = "Red dashed line: selected in all folds | Orange dashed line: selected in 4+ folds"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 14),
      axis.text.y = element_text(size = 10)
    )
}

#' Significancia y estabilidad de selección (Approach 3)
#'
#' @param var_frequency_approach3 Tabla con frecuencias y p-values por variable.
#' @return Un objeto ggplot.
#' @export
build_approach3_pvalue_stability <- function(var_frequency_approach3) {
  var_frequency_approach3 %>%
    filter(n_folds_selected >= 2) %>%  # Show variables selected in 2+ folds
    ggplot(aes(x = reorder(variable, mean_p_value),
               y = -log10(mean_p_value))) +
    geom_col(aes(fill = n_folds_selected), alpha = 0.7) +
    geom_errorbar(
      aes(ymin = -log10(max_p_value),
          ymax = -log10(min_p_value)),
      width = 0.3, alpha = 0.5
    ) +
    geom_point(aes(size = selection_rate), alpha = 0.8) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "red", size = 1) +
    geom_hline(yintercept = -log10(0.20), linetype = "dashed",
               color = "orange", size = 1) +
    scale_fill_gradient(
      low = "#FFF9C4",
      high = "#1976D2",
      name = "Folds\nSelected",
      breaks = c(2, 3, 4, 5)
    ) +
    scale_size_continuous(
      name = "Selection\nRate (%)",
      range = c(2, 6)
    ) +
    coord_flip() +
    labs(
      title = "Variable Significance and Selection Stability (Approach 3)",
      subtitle = "Variables selected in 2+ folds | Error bars show range of p-values across folds",
      x = "Clinical Variable",
      y = expression(-log[10](mean~p-value)),
      caption = "Red line: p=0.05 | Orange line: p=0.20 (screening threshold)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 14),
      axis.text.y = element_text(size = 10)
    )
}

#' Heatmap de selección de variables por fold (Approach 3)
#'
#' @param selection_matrix Matriz larga fold x variable con `selected` y `p_value`.
#' @param var_order Orden de variables (por estabilidad) para el eje y.
#' @return Un objeto ggplot.
#' @export
build_approach3_selection_heatmap <- function(selection_matrix, var_order) {
  selection_matrix %>%
    mutate(
      variable = factor(variable, levels = var_order),
      fold = factor(fold, levels = 1:5),
      p_value_label = ifelse(selected == 1,
                              sprintf("%.3f", p_value),
                              "")
    ) %>%
    ggplot(aes(x = fold, y = variable, fill = as.factor(selected))) +
    geom_tile(color = "white", size = 1) +
    geom_text(aes(label = p_value_label), size = 3, color = "black") +
    scale_fill_manual(
      values = c("0" = "gray90", "1" = "#1976D2"),
      labels = c("Not selected", "Selected"),
      name = ""
    ) +
    labs(
      title = "Variable Selection Pattern Across CV Folds (Approach 3)",
      subtitle = "P-values shown for selected variables",
      x = "Cross-Validation Fold",
      y = "Clinical Variable",
      caption = "Blue cells: variable selected | Numbers: p-value from univariate screening"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text.x = element_text(size = 11, face = "bold"),
      axis.text.y = element_text(size = 10),
      legend.position = "bottom",
      panel.grid = element_blank()
    )
}


# ==========================================================================
# ANCOM-BC2 — Frecuencia de selección de taxa
# ==========================================================================

#' Frecuencia de selección de taxa ANCOM-BC2 por fold
#'
#' @param taxa_frequency Tabla con `taxon`, `n_folds_selected`.
#' @return Un objeto ggplot.
#' @export
build_ancom_taxa_freq <- function(taxa_frequency) {
  taxa_frequency %>%
    mutate(
      stability = case_when(
        n_folds_selected == 5 ~ "Very Stable (5/5)",
        n_folds_selected >= 4 ~ "Stable (4/5)",
        n_folds_selected >= 3 ~ "Moderate (3/5)",
        n_folds_selected >= 2 ~ "Low-moderate (2/5)",
        TRUE ~ "Low (1/5)"
      )
    ) %>%
    ggplot(aes(x = reorder(taxon, n_folds_selected), y = n_folds_selected)) +
    geom_col(aes(fill = stability), alpha = 0.8) +
    geom_hline(yintercept = 5, linetype = "dashed", color = "red", size = 1) +
    scale_fill_manual(
      values = c(
        "Very Stable (5/5)" = "#0072B2",
        "Stable (4/5)" = "#56B4E9",
        "Moderate (3/5)" = "#F0E442",
        "Low-moderate (2/5)" = "#E69F00",
        "Low (1/5)" = "#D55E00"
      ),
      name = "Selection Stability"
    ) +
    coord_flip() +
    labs(
      title = "ANCOM-BC2 Taxa Selection Frequency Across CV Folds",
      subtitle = "Differentially abundant taxa identified independently in each fold",
      x = "Bacterial Genus",
      y = "Number of Folds Selected (out of 5)",
      caption = "Red dashed line: selected in all folds"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 14)
    )
}


# ==========================================================================
# Tabla de performance del nested CV (presentación, sin cómputo de métricas)
# ==========================================================================

#' Tabla kable de performance de las combinaciones de nested CV
#'
#' Formatea (sin recalcular) los `summary` ya calculados de `cv_results_all` en una
#' tabla HTML ordenada por AUROC. Pura presentación.
#'
#' @param cv_results_all Lista de resultados de [train_with_nested_cv()].
#' @return Un objeto kableExtra (tabla HTML).
#' @export
build_cv_results_table <- function(cv_results_all) {
  cv_performance <- map_dfr(cv_results_all, function(res) {
    tibble(
      Model = res$model_name,
      Approach = res$approach,
      Microbiome = res$microbiome,
      AUROC = sprintf("%.3f ± %.3f",
                      res$summary$AUROC_mean, res$summary$AUROC_sd),
      PRAUC = sprintf("%.3f ± %.3f",
                      res$summary$PRAUC_mean, res$summary$PRAUC_sd),
      Sensitivity = sprintf("%.3f ± %.3f",
                            res$summary$Sensitivity_mean, res$summary$Sensitivity_sd),
      Specificity = sprintf("%.3f ± %.3f",
                            res$summary$Specificity_mean, res$summary$Specificity_sd),
      Accuracy = sprintf("%.3f ± %.3f",
                         res$summary$Accuracy_mean, res$summary$Accuracy_sd),
      Bal_Accuracy = sprintf("%.3f ± %.3f",
                             res$summary$Balanced_Accuracy_mean,
                             res$summary$Balanced_Accuracy_sd),
      Threshold = sprintf("%.3f ± %.3f",
                          res$summary$threshold_mean, res$summary$threshold_sd),
      Youden = sprintf("%.3f ± %.3f",
                       res$summary$Youden_mean, res$summary$Youden_sd),
      AUROC_mean = res$summary$AUROC_mean,  # For sorting
      PRAUC_mean = res$summary$PRAUC_mean   # For secondary sorting
    )
  }) %>%
    arrange(desc(AUROC_mean), desc(PRAUC_mean)) %>%
    select(-AUROC_mean, -PRAUC_mean)

  # Display table
  kable(cv_performance,
        caption = "Model Performance with Nested Cross-Validation (Mean ± SD across 5 folds)",
        format = "html",
        align = c("l", "l", "l", rep("r", 8))) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                  full_width = FALSE,
                  font_size = 10) %>%
    row_spec(0, bold = TRUE, background = "#E3F2FD") %>%
    column_spec(1, bold = TRUE, width = "7em") %>%
    column_spec(2, width = "10em") %>%
    row_spec(1:3, background = "#C8E6C9") %>%
    footnote(
      general = c(
        "Thresholds optimized on inner validation sets using Youden Index.",
        "Performance metrics calculated on independent outer test folds.",
        "This approach avoids data leakage and provides unbiased estimates."
      ),
      general_title = "Notes:",
      footnote_as_chunk = TRUE
    )
}


# ==========================================================================
# Curvas ROC / PR (a partir de curvas ALMACENADAS, sin re-entrenar)
# ==========================================================================

#' Curvas ROC de los top-3 modelos (medias interpoladas + IC 95%)
#'
#' @param roc_summary Curvas ROC medias interpoladas por modelo/rank.
#' @param labels_ordered Etiquetas de modelo en orden de rank.
#' @param colors Colores por rank (1,2,3).
#' @return Un objeto ggplot.
#' @export
build_roc_curves <- function(roc_summary, labels_ordered, colors) {
  ggplot() +
    # Confidence bands
    geom_ribbon(
      data = roc_summary,
      aes(x = fpr_interp,
          ymin = sens_lower,
          ymax = sens_upper,
          fill = factor(model_rank, levels = 1:3)),
      alpha = 0.2
    ) +
    # Lines for ranks 2 and 3 (thinner)
    geom_line(
      data = roc_summary %>% filter(model_rank %in% c(2, 3)),
      aes(x = fpr_interp,
          y = sens_mean,
          color = factor(model_rank, levels = 1:3)),
      linewidth = 1.0,
      alpha = 0.8
    ) +
    # Line for rank 1 (thicker, winner)
    geom_line(
      data = roc_summary %>% filter(model_rank == 1),
      aes(x = fpr_interp,
          y = sens_mean,
          color = factor(model_rank, levels = 1:3)),
      linewidth = 1.5,
      alpha = 1.0
    ) +
    # Diagonal reference
    geom_abline(
      intercept = 0, slope = 1,
      linetype = "dashed",
      color = "gray40",
      linewidth = 0.8
    ) +
    # Color scales with explicit factor levels
    scale_color_manual(
      values = colors,
      labels = labels_ordered,
      name = "",
      breaks = factor(1:3, levels = 1:3)
    ) +
    scale_fill_manual(
      values = colors,
      labels = labels_ordered,
      name = "",
      breaks = factor(1:3, levels = 1:3)
    ) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
      title = "ROC Curves: Top 3 Models",
      subtitle = "5-fold CV with 95% CI",
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 9),
      plot.title = element_text(face = "bold", size = 15),
      panel.grid.minor = element_blank()
    ) +
    guides(
      color = guide_legend(ncol = 1, override.aes = list(linewidth = 2)),
      fill = guide_legend(ncol = 1)
    )
}

#' Curvas Precision-Recall de los top-3 modelos (medias interpoladas + IC 95%)
#'
#' @param pr_summary Curvas PR medias interpoladas por modelo/rank.
#' @param labels_ordered Etiquetas de modelo en orden de rank.
#' @param colors Colores por rank (1,2,3).
#' @param prevalence Prevalencia (línea base de precisión).
#' @return Un objeto ggplot.
#' @export
build_pr_curves <- function(pr_summary, labels_ordered, colors, prevalence) {
  ggplot() +
    # Confidence bands
    geom_ribbon(
      data = pr_summary,
      aes(x = recall_interp,
          ymin = prec_lower,
          ymax = prec_upper,
          fill = factor(model_rank, levels = 1:3)),
      alpha = 0.2
    ) +
    # Lines for ranks 2 and 3
    geom_line(
      data = pr_summary %>% filter(model_rank %in% c(2, 3)),
      aes(x = recall_interp,
          y = prec_mean,
          color = factor(model_rank, levels = 1:3)),
      linewidth = 1.0,
      alpha = 0.8
    ) +
    # Line for rank 1 (winner)
    geom_line(
      data = pr_summary %>% filter(model_rank == 1),
      aes(x = recall_interp,
          y = prec_mean,
          color = factor(model_rank, levels = 1:3)),
      linewidth = 1.5,
      alpha = 1.0
    ) +
    # Baseline
    geom_hline(
      yintercept = prevalence,
      linetype = "dashed",
      color = "gray40",
      linewidth = 0.8
    ) +
    annotate(
      "text",
      x = 0.85, y = prevalence + 0.03,
      label = sprintf("Baseline (prevalence = %.3f)", prevalence),
      size = 3.5, color = "gray40"
    ) +
    # Color scales
    scale_color_manual(
      values = colors,
      labels = labels_ordered,
      name = "",
      breaks = factor(1:3, levels = 1:3)
    ) +
    scale_fill_manual(
      values = colors,
      labels = labels_ordered,
      name = "",
      breaks = factor(1:3, levels = 1:3)
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
      title = "Precision-Recall Curves: Top 3 Models",
      subtitle = "5-fold CV with 95% CI",
      x = "Recall (Sensitivity)",
      y = "Precision (PPV)"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 9),
      plot.title = element_text(face = "bold", size = 15),
      panel.grid.minor = element_blank()
    ) +
    guides(
      color = guide_legend(ncol = 1, override.aes = list(linewidth = 2)),
      fill = guide_legend(ncol = 1)
    )
}


# ==========================================================================
# FIGURA 5 — Importancia de features y estabilidad de selección
# ==========================================================================

#' Panel A de la Figura 5: importancia de features (top 20, Gini)
#'
#' @param top20_importance Top-20 features con `Variable_clean`, `Importance`, `feature_type`.
#' @return Un objeto ggplot.
#' @export
build_figure5_importance <- function(top20_importance) {
  top20_importance %>%
    ggplot(aes(x = reorder(Variable_clean, Importance),
               y = Importance,
               fill = feature_type)) +

    # Bars
    geom_col(alpha = 0.85, width = 0.75, color = "white", linewidth = 0.3) +

    # Flip coordinates
    coord_flip() +

    # Color scale
    scale_fill_manual(
      values = c(
        "Clinical" = "#E74C3C",
        "Microbiome" = "#3498DB",
        "Diversity" = "#2ECC71"
      ),
      name = "Feature Type"
    ) +

    # Labels
    labs(
      x = NULL,
      y = "Variable Importance\n(Mean Decrease in Gini Impurity)",
      title = NULL
    ) +

    # Theme
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 10)),
      plot.caption = element_text(size = 9, color = "gray50", hjust = 0, margin = margin(t = 10)),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    )
}

#' Panel B de la Figura 5: estabilidad de selección clínica (Approach 3)
#'
#' @param var_frequency Tabla con `variable`, `n_folds_selected`, `mean_p_value`.
#' @return Un objeto ggplot.
#' @export
build_figure5_stability <- function(var_frequency) {
  var_frequency_plot <- var_frequency %>%
    mutate(
      # P-value category (SHORTENED LABELS for legend)
      p_category = case_when(
        mean_p_value < 0.01 ~ "p<0.01",
        mean_p_value < 0.05 ~ "p<0.05",
        mean_p_value < 0.10 ~ "p<0.10",
        mean_p_value < 0.20 ~ "p<0.20",
        TRUE ~ "p≥0.20"
      ),
      p_category = factor(p_category, levels = c(
        "p<0.01", "p<0.05", "p<0.10", "p<0.20", "p≥0.20"
      )),
      variable_clean = case_when(
  # Existing (keep)
  variable == "edad_cronologicamujer" ~ "Maternal Age",
  variable == "imc_visita" ~ "BMI at Visit",
  variable == "peso_pregestacional_kg" ~ "Pre-pregnancy Weight",
  variable == "peso_kg" ~ "Weight at Visit",
  variable == "extreme_bmi" ~ "Extreme BMI",
  variable == "underweight" ~ "Underweight",
  variable == "obese" ~ "Obese",
  variable == "hemoglobin_alti_adj" ~ "Hemoglobin (altitude-adj)",
  variable == "hemoglobin_g_dl" ~ "Hemoglobin",
  variable == "low_education" ~ "Low Education",
  variable == "unmarried" ~ "Unmarried Status",
  variable == "ses_risk_score" ~ "SES Risk Score",
  variable == "talla_mujer_cm" ~ "Height",
  variable == "imc_pregestacional" ~ "Pre-pregnancy BMI",
  variable == "rpm_preterm" ~ "Preterm PROM",
  variable == "rpm" ~ "PROM",
  variable == "diabetes_gest" ~ "Gestational Diabetes",
  variable == "complication_count" ~ "Complication Count",
  variable == "any_complication" ~ "Any Complication",
  variable == "anemia_visita" ~ "Anemia at Visit",
  variable == "folic_ac_correg" ~ "Folic Acid Intake",
  variable == "vitaminsup" ~ "Vitamin Supplementation",
  variable == "comp1tribleed" ~ "First Trimester Bleeding",
  variable == "extreme_age" ~ "Extreme Maternal Age",
  variable == "workouthome" ~ "Works Outside Home",
  variable == "maritalstat" ~ "Marital Status",
  variable == "oligohidramnios" ~ "Oligohydramnios",
  variable == "compvaginf" ~ "Cervicovaginal Infection",
  variable == "dietsuppl3mon" ~ "Dietary Supplement Use",
  variable == "sdg_visita" ~ "Gestational Age at Visit",
  variable == "age_risk_category" ~ "Age Risk Category",
  variable == "comppreeclam" ~ "Preeclampsia",
  variable == "compsexualinf" ~ "Sexually Transmitted Infection",
  variable == "rciu" ~ "Intrauterine Growth Restriction",
        TRUE ~ variable
      )
    )

  # Create plot
  var_frequency_plot %>%
    ggplot(aes(x = reorder(variable_clean, n_folds_selected),
               y = n_folds_selected)) +

    # Bars colored by p-value strength
    geom_col(aes(fill = p_category), alpha = 0.85, width = 0.75,
             color = "white", linewidth = 0.3) +

    # Color scale for p-value categories
    scale_fill_manual(
        name = "Univariate Association Strength",
        values = c(
          "p<0.01" = "#0072B2",
          "p<0.05" = "#56B4E9",
          "p<0.10" = "#F0E442",
          "p<0.20" = "#E69F00",
          "p≥0.20" = "#D55E00"
        ),
      labels = c("p<0.01", "p<0.05", "p<0.10", "p<0.20", "p≥0.20"),
      guide = guide_legend(
        direction = "horizontal",
        title.position = "top",
        title.hjust = 0.5,
        label.position = "bottom",
        nrow = 1
      )
    ) +

    # Y-axis
    scale_y_continuous(
      breaks = 1:5,
      limits = c(0, 5.5),
      expand = c(0, 0)
    ) +

    # Flip coordinates
    coord_flip() +

    # Labels
    labs(
      x = NULL,
      y = "Number of CV Folds Selected (out of 5)",
      title = NULL
    ) +

    # Theme with HORIZONTAL LEGEND
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      legend.key.width = unit(2, "cm"),
      legend.key.height = unit(0.5, "cm"),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 10)),
      plot.caption = element_text(size = 8.5, color = "gray50", hjust = 0,
                                   margin = margin(t = 10), lineheight = 1.2),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    )
}

#' Figura 5 combinada (panel A | panel B)
#'
#' @param panel_a Importancia (build_figure5_importance()).
#' @param panel_b Estabilidad (build_figure5_stability()).
#' @param best_rf_info Fila con `auroc_mean`, `auroc_sd` del mejor RF (para el título).
#' @return Un objeto patchwork.
#' @export
build_figure5_combined <- function(panel_a, panel_b, best_rf_info) {
  (panel_a | panel_b) +
    plot_annotation(
      tag_levels = 'A',
      title = sprintf(
        "Feature Importance and Selection Stability in Best-Performing Random Forest Model\n(Data-driven features with full microbiome | AUROC = %.3f ± %.3f)",
        best_rf_info$auroc_mean,
        best_rf_info$auroc_sd
      ),
      theme = theme(
        plot.title = element_text(face = "bold", size = 13, hjust = 0, lineheight = 1.2),
        plot.margin = margin(10, 10, 10, 10)
      )
    ) +
    plot_layout(widths = c(1, 1))
}


# ==========================================================================
# Comparación de performance (boxplots AUROC + trade-off Sens/Spec)
# ==========================================================================

#' Boxplots de AUROC por combinación de modelo
#'
#' @param plot_data_all Métricas por fold de todas las combinaciones.
#' @return Un objeto ggplot.
#' @export
build_performance_auroc <- function(plot_data_all) {
  ggplot(plot_data_all,
         aes(x = reorder(paste(model, approach), -AUROC),
             y = AUROC, fill = microbiome)) +
    geom_boxplot(alpha = 0.7) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
    scale_fill_manual(values = c("ANCOM_Taxa" = "#66C2A5",
                                  "Full_Microbiome" = "#FC8D62")) +
    labs(
      title = "Model Performance: AUROC Across All Combinations",
      subtitle = "12 model combinations with 5-fold nested cross-validation",
      x = "Model + Approach",
      y = "AUROC",
      fill = "Microbiome"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "top"
    )
}

#' Scatter de trade-off Sensibilidad vs Especificidad
#'
#' @param summary_data Resumen por combinación (AUROC/Sens/Spec + sd).
#' @return Un objeto ggplot.
#' @export
build_performance_sens_spec <- function(summary_data) {
  ggplot(summary_data,
         aes(x = Specificity, y = Sensitivity,
             color = AUROC, shape = model)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_errorbar(aes(ymin = Sensitivity - Sens_sd,
                      ymax = Sensitivity + Sens_sd), alpha = 0.3) +
    geom_errorbarh(aes(xmin = Specificity - Spec_sd,
                       xmax = Specificity + Spec_sd), alpha = 0.3) +
    scale_color_gradient2(low = "#D73027", mid = "#FEE090", high = "#1A9850",
                          midpoint = 0.7, limits = c(0.5, 1)) +
    scale_shape_manual(values = c(16, 17)) +
    labs(
      title = "Sensitivity vs Specificity Trade-off",
      subtitle = "Error bars show ±1 SD across folds",
      x = "Specificity",
      y = "Sensitivity",
      color = "AUROC",
      shape = "Model"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold")) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1))
}
