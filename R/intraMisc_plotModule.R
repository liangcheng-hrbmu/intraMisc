###########################################################
#' Visual Representation of Intracellular Microbiome Diversity.
#'
#' @param seurat.object seurat object after running \code{\link{runSeurat}} and \code{\link{calcAlpha}}
#' @param alpha.index Default: \code{"Simpson"},
#' string or vector, alpha diversity index name.
#'
#' @param group.by Default: \code{"orig.ident"},
#' Name of seurat object metadata columns to group (color) cells by.
#'
#' @param reduction Default: \code{NULL}.
#' Which dimensionality reduction to use. Selected from \code{"umap"},  \code{"tsne"}, or \code{"pca"}.
#'
#' @param custom.color Default: \code{NULL}, which will use 20 built-in colors \code{c('#D97894','#7C98CE','#BBCB54', '#A7BED3',
#' '#C9A2E1', '#748977', '#BEBE9E', '#4CD28F',
#' '#CFAB84', '#E3E15C', '#DFC7CF', '#6BB8BF',
#' '#CBCBCB', '#9696eb', '#DED1AD', '#E29F58',
#' '#D8A1CE', '#B16455', '#ad9bf4', '#BEB68D')}.
#' Users could provide a vector of colors, each color corresponds to one group mapped to param \code{"group.by"}
#'
#' @param tax_level Default: \code{"species"} taxonomy level.
#'
#' @import ggplot2
#' @importFrom stringr str_to_title
#'
#' @return ggplot2 format graph
#' @export
#'
#' @examples
#' data("seurat.object")
#' # seurat.object <- runSeurat(seurat.object)
#' seurat.object <- calcAlpha(seurat.object)
#' plotAlpha(seurat.object)
#' # plotAlpha(seurat.object, alpha.index = "Shannon", reduction = "umap")
plotAlpha <- function(seurat.object,
                      alpha.index = "Simpson",
                      tax_level = "species",
                      group.by = "orig.ident",
                      reduction = NULL,
                      custom.color = NULL) {
  if (is.null(custom.color)) {
    custom.color <- c(
      "#D97894", "#7C98CE", "#BBCB54", "#A7BED3", "#C9A2E1",
      "#748977", "#BEBE9E", "#4CD28F", "#CFAB84", "#E3E15C",
      "#DFC7CF", "#6BB8BF", "#CBCBCB", "#9696eb", "#DED1AD",
      "#E29F58", "#D8A1CE", "#B16455", "#ad9bf4", "#BEB68D"
    )
  }
  if (!paste0(alpha.index, "_", tax_level) %in% colnames(seurat.object@meta.data)) {
    stop('Run "seurat.object <- calcAlpha(seurat.object)", then plot alpha diversity.')
  }
  if (is.null(reduction)) {
    p_alpha <-
      ggplot(seurat.object@meta.data, aes(
        seurat.object@meta.data[, group.by],
        seurat.object@meta.data[, paste0(alpha.index, "_", tax_level)]
      )) +
      geom_jitter(aes(color = seurat.object@meta.data[, group.by]), width = 0.1, size = 0.6) +
      geom_boxplot(aes(color = seurat.object@meta.data[, group.by]),
        width = .1,
        position = position_nudge(x = 0.25),
        size = 0.3
      ) +
      labs(title = paste0(alpha.index, " Diversity - ", tax_level) %>% str_to_title(), x = group.by, y = alpha.index) +
      theme_bw() +
      theme(
        text = element_text(family = "serif"),
        aspect.ratio = NULL,
        panel.border = element_rect(fill = NA, linewidth = 0.3, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold")
      ) +
      scale_fill_manual(
        name = str_to_title(gsub("_", " ", group.by)),
        values = custom.color
      ) +
      scale_color_manual(
        name = str_to_title(gsub("_", " ", group.by)),
        values = custom.color
      )
  } else if (reduction %in% c("tsne", "umap", "pca")) {
    p_alpha <- Seurat::FeaturePlot(seurat.object,
      features = paste0(alpha.index, "_", tax_level),
      cols = c("grey", "red"),
      reduction = reduction
    ) +
      theme(aspect.ratio = 1) +
      labs(title = paste0(alpha.index, " - ", tax_level))
  }
  return(p_alpha)
}


#' Visualize Feature Abundance or Expression Across Samples or Groups.
#'
#' @param seurat.object seurat object
#' @param tax_name Default: \code{NULL}. Microbe tax_name from different taxonomy levels
#'   (e.g. [\code{"Escherichia_coli", tax_level = "species"}],
#'   [\code{"Bacillota", tax_level = "phylum"]}.)
#'
#' @param tax_level Default: \code{"species"}. Taxonomy level of param \strong{tax_name}.
#'
#' @param features Default: \code{c("nCount_RNA", "nFeature_RNA",
#' "nCount_Microbe", "nFeature_Microbe")}.
#' Vector of features to plot.
#' Details see \strong{features} in \code{\link[Seurat]{FeaturePlot}}.\cr
#' Features can come from:
#' \itemize{
#'   \item feature (e.g. a gene name - "MS4A1")
#'   \item A column name from meta.data (e.g. mitochondrial percentage - \code{"percent.mito"}, \code{"nCount_Microbe"})
#'   \item A column name from a DimReduc object corresponding to the cell embedding values (e.g. the PC 1 scores - "PC_1")
#' }
#'
#' @param reduction Defalut: \code{"umap"}. Details see \strong{reduction} in \code{\link[Seurat]{FeaturePlot}}.
#' Which dimensionality reduction to use.
#' If not specified, first searches for umap, then tsne, then pca.
#'
#' @param ncol Default: \code{2}.
#' (optional) Number of rows in the plot grid.
#' Details see \strong{ncol} in \code{\link[cowplot]{plot_grid}}.
#'
#' @param nrow Default: \code{2}.
#' (optional) Number of columns in the plot grid.
#' Details see \strong{nrow} in \code{\link[cowplot]{plot_grid}}.
#'
#' @param labels Default: \code{NULL}.
#' (optional) List of labels to be added to the plots.
#' Details see \strong{labels} in \code{\link[cowplot]{plot_grid}}.
#' \describe{
#' \item{\code{"AUTO"}}{auto-generate upper-case labels, "A", "B"...}
#' \item{\code{"auto"}}{auto-generate lower-case labels, "a", "b"...}
#' }
#'
#' @param ... additional params from function \strong{FeaturePlot}.
#' Details params in \code{\link[Seurat]{FeaturePlot}}.
#' @param spt.sc.theme Default: \code{FALSE}.
#' If the data is spatial transcriptomics,
#' should the visualization theme for single-cell sequencing data analysis be used?
#'
#' @return feature plot
#' @export
#'
#' @examples
#' data("seurat.object")
#' seurat.object <- runSeurat(seurat.object)
#' plotFeatures(seurat.object, labels = "auto")
#' # plotFeatures(seurat.object, features = c("PLA2G2A", "CDA", "HSPG2"), labels = "AUTO")
#' # plotFeatures(seurat.object, tax_name = "Escherichia_coli")
#' # plotFeatures(seurat.object, tax_name = "Bacillota", tax_level = "phylum")
#'
plotFeatures <- function(seurat.object,
                         tax_name = NULL,
                         tax_level = "species",
                         features = c("nCount_RNA", "nFeature_RNA", "nCount_Microbe", "nFeature_Microbe"),
                         reduction = "umap",
                         ncol = 2,
                         nrow = 2,
                         labels = NULL,
                         spt.sc.theme = FALSE,
                         ...) {
  if ("RNA" %in% names(seurat.object@assays)) {
    if (!is.null(tax_name)) {
      if (is.null(seurat.object@misc$microbe$phyloseq.object)) {
        microbe.phyloseq.object <- seurat.object@misc$microbe$seurat.object %>%
          seurat2phyloseq()
      } else {
        microbe.phyloseq.object <- seurat.object@misc$microbe$phyloseq.object
      }
      if (tax_level == "species") {
        seurat.object@meta.data$microbe_feature <- microbe.phyloseq.object %>%
          phyloseq::get_sample(tax_name)
      } else {
        seurat.object@meta.data$microbe_feature <- phyloseq_select_tax(microbe.phyloseq.object, tax_name, tax_level)
      }
      p <- seurat.object %>%
        Seurat::FeaturePlot(features = c("microbe_feature"), reduction = reduction, ...) +
        theme(aspect.ratio = 1) +
        labs(fill = paste0(tax_name, " - ", tax_level))
    } else {
      p_list <- list()
      for (i in seq_along(features)) {
        p_list[[i]] <- seurat.object %>%
          Seurat::FeaturePlot(features = features[i], reduction = reduction, ...) +
          theme(aspect.ratio = 1)
      }
      p <- cowplot::plot_grid(plotlist = p_list, ncol = ncol, nrow = nrow, labels = labels)
    }
  } else if ("Spatial" %in% names(seurat.object@assays)) {
    if (spt.sc.theme) {
      spatial_corr <- SeuratObject::GetTissueCoordinates(seurat.object, cols = c("imagerow", "imagecol"))[, c("y", "x")]
      colnames(spatial_corr) <- c("s_1", "s_2")
      spatial_corr <- as.matrix(spatial_corr)
      seurat.object[["spatial"]] <- SeuratObject::CreateDimReducObject(embeddings = spatial_corr, key = "s_", assay = "Spatial")
    }

    if (!is.null(tax_name)) {
      if (is.null(seurat.object@misc$microbe$phyloseq.object)) {
        microbe.phyloseq.object <- seurat.object@misc$microbe$seurat.object %>%
          seurat2phyloseq()
      } else {
        microbe.phyloseq.object <- seurat.object@misc$microbe$phyloseq.object
      }
      if (tax_level == "species") {
        seurat.object@meta.data$microbe_feature <- microbe.phyloseq.object %>%
          phyloseq::get_sample(tax_name)
      } else {
        seurat.object@meta.data$microbe_feature <- phyloseq_select_tax(microbe.phyloseq.object, tax_name, tax_level)
      }
      if (spt.sc.theme) {
        p <- seurat.object %>%
          Seurat::FeaturePlot(
            features = c("microbe_feature"), alpha = 0.8,
            reduction = "spatial"
          ) +
          scale_y_reverse() +
          scale_colour_gradientn(colors = Seurat::CustomPalette(low = "#313695", high = "#A50026", mid = "#FFFFBF", k = 100)) +
          theme_void() +
          theme(
            aspect.ratio = 1,
            legend.position = "right"
          ) +
          labs(colour = paste0(tax_name, " - ", tax_level), title = NULL)
      } else {
        p <- seurat.object %>%
          Seurat::SpatialFeaturePlot(
            features = c("microbe_feature"),
            alpha = c(1, 0.5), pt.size = 2.5
          ) +
          theme(
            aspect.ratio = 1,
            legend.position = "right"
          ) +
          labs(fill = paste0(tax_name, " - ", tax_level), title = NULL)
      }
    } else {
      p_list <- list()
      for (i in seq_along(features)) {
        features[features == "nCount_RNA"] <- "nCount_Spatial"
        features[features == "nFeature_RNA"] <- "nFeature_Spatial"
        if (spt.sc.theme) {
          p_list[[i]] <- seurat.object %>%
            Seurat::FeaturePlot(
              features = features[i], alpha = 0.8,
              reduction = "spatial"
            ) +
            scale_y_reverse() +
            scale_colour_gradientn(colors = Seurat::CustomPalette(low = "#313695", high = "#A50026", mid = "#FFFFBF", k = 100)) +
            theme_void() +
            theme(
              aspect.ratio = 1,
              legend.position = "right"
            ) +
            labs(colour = features[i], title = NULL)
        } else {
          p_list[[i]] <- seurat.object %>%
            Seurat::SpatialFeaturePlot(
              features = features[i],
              alpha = c(1, 0.5), pt.size = 2.5
            ) +
            theme(
              aspect.ratio = 1,
              legend.position = "right"
            ) +
            labs(fill = features[i], title = NULL)
        }
      }
      p <- cowplot::plot_grid(plotlist = p_list, ncol = ncol, nrow = nrow, labels = labels)
    }
  }

  return(p)
}

if (getRversion() >= "2.15.1") utils::globalVariables(c("group", "Abundance", "taxName"))
#' Visualize Taxonomic Composition Across Samples Using Stacked Bar Charts.
#'
#' @param seurat.object seurat object.
#'
#' @param tax_level Default: \code{"species"}. Taxonomy level of param \strong{feature}.
#'
#' @param group.by Default: \code{"orig.ident"}.
#' Name of seurat object metadata columns to group (color) cells by.
#'
#' @param custom.color Default: \code{NULL}. which will use 20 built-in colors \code{c('#D97894','#7C98CE','#BBCB54', '#A7BED3',
#' '#C9A2E1', '#748977', '#BEBE9E', '#4CD28F',
#' '#CFAB84', '#E3E15C', '#DFC7CF', '#6BB8BF',
#' '#CBCBCB', '#9696eb', '#DED1AD', '#E29F58',
#' '#D8A1CE', '#B16455', '#ad9bf4', '#BEB68D')}.
#' Users could provide a vector of colors, each color corresponds to one group mapped to param \code{"group.by"}.
#'
#' @param cluster_group Default: \code{FALSE}. Weather cluster the groups with a tree based on bray distance.
#'
#' @param top_n Default: \code{20}. Numbers of enriched taxonomy to show.
#'
#' @return Stark Bar
#' @export
#'
#' @examples
#' data("seurat.object")
#' plotStackedBar(seurat.object)
#' # seurat.object <- runSeurat(seurat.object)
#' # plotStackedBar(seurat.object, group.by = "seurat_clusters", tax_level = "phylum")
#'
plotStackedBar <- function(seurat.object,
                           tax_level = "species",
                           group.by = "orig.ident",
                           custom.color = NULL,
                           cluster_group = FALSE,
                           top_n = 20) {
  check_package("vegan")
  if (is.null(custom.color)) {
    custom.color <- c(
      "#D97894", "#7C98CE", "#BBCB54", "#A7BED3", "#C9A2E1",
      "#748977", "#BEBE9E", "#4CD28F", "#CFAB84", "#E3E15C",
      "#DFC7CF", "#6BB8BF", "#CBCBCB", "#9696eb", "#DED1AD",
      "#E29F58", "#D8A1CE", "#B16455", "#ad9bf4", "#BEB68D",
      "#B6B6B6"
    )
  }
  if (is.null(seurat.object@misc$microbe$phyloseq.object)) {
    microbe.phyloseq.object <- seurat.object@misc$microbe$seurat.object %>%
      seurat2phyloseq()
  } else {
    microbe.phyloseq.object <- seurat.object@misc$microbe$phyloseq.object
  }

  microbe.phyloseq.object.tax_level <- microbe.phyloseq.object %>%
    phyloseq_select_level(tax_level = tax_level)
  otu_tab.norm <- microbe.phyloseq.object.tax_level@otu_table@.Data %>%
    as.data.frame() %>%
    t() %>%
    stats::aggregate(by = list(seurat.object@meta.data[, group.by]), sum) %>%
    {
      rownames(.) <- .$Group.1
      .[["Group.1"]] <- NULL
      t(.)
    } %>%
    scale(center = FALSE, scale = colSums(.)) %>%
    as.data.frame()
  otu_tab.norm <- otu_tab.norm[, colSums(is.na(otu_tab.norm)) == 0]
  otu_tab.norm <- otu_tab.norm[order(rowSums(otu_tab.norm), decreasing = TRUE), ]
  if (cluster_group) {
    bray_distance <- otu_tab.norm %>%
      t() %>%
      vegan::vegdist(method = "bray")
    tree <- stats::hclust(bray_distance, method = "average")
    # let the tree face right
    tree$height <- -tree$height
    p_tree <- (tree %>% stats::as.dendrogram() %>%
      ggdendro::ggdendrogram()) +
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = margin(t = 0.05, r = 0, b = 0.05, l = 0, unit = "pt")
      ) +
      coord_flip() +
      theme(aspect.ratio = 1)
  }

  if (top_n < nrow(otu_tab.norm)) {
    melData <- otu_tab.norm[1:top_n, ]
    melData["Others", ] <- colSums(otu_tab.norm[(top_n + 1):nrow(otu_tab.norm), ])
  } else {
    melData <- otu_tab.norm
  }
  melData$taxName <- rownames(melData)
  plotData <- reshape2::melt(melData,
    id.vars = c("taxName"),
    measure.vars = colnames(melData)[1:(ncol(melData) - 1)],
    variable.name = "group",
    value.name = "Abundance"
  )
  plotData$taxName <- factor(plotData$taxName, levels = unique(plotData$taxName))
  if (cluster_group) {
    plotData$group <- factor(plotData$group, levels = tree$labels[tree$order])
  }
  p_bar <- ggplot(plotData, aes(x = group)) +
    geom_bar(aes(y = Abundance, fill = taxName),
      stat = "identity", width = 0.85,
      position = position_stack(reverse = T)
    ) +
    labs(title = group.by, x = NULL, y = "Relative Abundance") +
    coord_flip() +
    theme(
      text = element_text(family = "serif"),
      panel.background = element_blank(),
      panel.spacing = unit(0.1, "lines"),
      axis.ticks.x = element_blank(),
      axis.line = element_line(colour = "black", lineend = 0.5),
      # axis.text.x = element_blank(),
      strip.text.x = element_text(size = 10, face = "bold", family = "serif"),
      strip.background.x = element_rect(fill = "#F6F6F6", colour = "#E0E0E0"),
      axis.title.y = element_text(size = 14, face = "bold"),
      # legend.position = "bottom",
      legend.key.size = unit(15, "pt"),
      legend.title = element_text(size = 12, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 0.05, r = 0, b = 0.05, l = 0, unit = "pt"),
      plot.title = element_text(hjust = -1),
      aspect.ratio = 1 / 10 * length(unique(plotData$group))
    ) +
    scale_fill_manual(name = tax_level, values = custom.color) +
    guides(
      fill = guide_legend(
        ncol = 2,
        byrow = TRUE,
        reverse = FALSE
      )
    )
  if (cluster_group) {
    p <- p_tree + p_bar
  } else {
    p <- p_bar
  }
  return(p)
}

#' Core Microbiome Visualization Based on Prevalence and Abundance Thresholds.
#'
#' @param seurat.object seurat object
#' @param tax_level Default: \code{"species"}. Taxonomy level.
#' @param group.by Default: \code{"orig.ident"}.
#' Name of seurat object metadata columns to group (color) cells by.
#'
#' @param color Default: \code{c("white", "red")}.
#' Colors of the core microbe heatmap.
#' User could provide a vector of colours to use for n-colour gradient.
#' E.g: \code{"colorRampPalette(c("white", "red"))(256)"}
#'
#' @return heatmap
#' @export
#'
#' @examples
#' data("seurat.object")
#' plotCore(seurat.object)
#'
plotCore <- function(seurat.object,
                     tax_level = "species",
                     group.by = "orig.ident",
                     color = c("#EFEFEF", "#A50026") # colorRampPalette(c("white", "red"))(256)
) {
  if (is.null(seurat.object@misc$microbe$phyloseq.object)) {
    microbe.phyloseq.object <- seurat.object@misc$microbe$seurat.object %>%
      seurat2phyloseq()
  } else {
    microbe.phyloseq.object <- seurat.object@misc$microbe$phyloseq.object
  }
  microbe.phyloseq.object.tax_level <- microbe.phyloseq.object %>%
    phyloseq_select_level(tax_level = tax_level)

  otu_tab.norm <- microbe.phyloseq.object.tax_level@otu_table@.Data %>%
    t() %>%
    stats::aggregate(by = list(seurat.object@meta.data[, group.by]), sum) %>%
    {
      rownames(.) <- .[, "Group.1"]
      .[, 2:ncol(.)]
    } %>%
    t() %>%
    as.data.frame()
  otu_tab.norm <- otu_tab.norm[rowSums(otu_tab.norm) != 0, ]
  otu_tab.norm <- otu_tab.norm[, colSums(otu_tab.norm) != 0]
  otu_tab.norm <- otu_tab.norm %>% scale(., center = F, scale = colSums(.))

  prevalences <- seq(5e-2, max(otu_tab.norm), 5e-2)
  detections <- round(10^seq(log10(0.01), log10(max(otu_tab.norm)), length = 8), 3)

  # Also define gray color palette
  p_core <- microbiome::plot_core(otu_tab.norm,
    plot.type = "heatmap",
    colours = color,
    prevalences = prevalences,
    detections = detections,
    min.prevalence = 0.50
  ) +
    labs(
      x = "Detection Threshold\n(Relative Abundance (%))",
      y = tax_level,
      title = paste0("Core microbe analysis among ", group.by)
    ) +

    # Adjusts axis text size and legend bar height
    theme(
      axis.text.y = element_text(size = 8, face = "italic"),
      axis.text.x.bottom = element_text(size = 8),
      axis.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 10),
      plot.title = element_text(hjust = 0.5)
    )
  return(p_core)
}

# Temporarily commented out  plotTree function due to orphan package dependencies

# if (getRversion() >= "2.15.1") utils::
#     globalVariables(c("taxon_names", "n_obs"))
# # ' Title
# # '
# # ' @param seurat.object  seurat object
# # '
# # ' @return Tree diagram
# # ' @export
# # '
# # ' @examples
# # ' #
# # '
# plotTree < - function(seurat.object)
# {
#     # check_package("taxize")
#     # read microbe abundance data
#     microbe.seurat.object < - seurat.object @ misc$microbe$seurat.object % > % SeuratObject::JoinLayers()
# abundance < - microbe.seurat.object @ assays$Microbe$counts % > % as.data.frame()
# abundance[, colSums(abundance) != 0] < - abundance[, colSums(abundance) != 0] % > %scale(center=FALSE, scale=colSums(.))
#
#
#
#
# rownames(abundance) < - rownames(abundance) % > %
# stringr::str_replace_all("#", "\\|") % > %
# stringr::str_replace_all("-", "_")
# abundance$Classification < - rownames(abundance)
#
# x < - metacoder::parse_tax_data(abundance,
#                                 class_cols="Classification", class_sep="|",
#                                 class_key=c(tax_rank="taxon_rank", tax_name="taxon_name"),
#                                 class_regex="^(.+)__(.+)$"
#                                 )
#
# # calculate the abundance of taxa.
# x$data$tax_abund < - metacoder::calc_taxon_abund(x, "tax_data",
#                                                  cols=colnames(microbe.seurat.object)
#                                                  )
#
# # # the number of microbial occurrences was counted in groups.
# # x$data$tax_occ <- metacoder::calc_n_samples(x, "tax_abund",
# #                                             cols = colnames(microbe.seurat.object),
# #                                             groups = seurat.object@meta.data[, "orig.ident"])
#
#
# x$data$tax_occ < - metacoder::calc_group_mean(x, "tax_abund",
#                                               cols=colnames(microbe.seurat.object),
#                                               groups=seurat.object @ meta.data[, "orig.ident"]
# )
# group1 < - colnames(x$data$tax_occ)[2]
# group2 < - colnames(x$data$tax_occ)[3]
# x$data$tax_occ$group < - NA
# for (group in colnames(x$data$tax_occ)[2:ncol(x$data$tax_occ)]) {
# x$data$tax_occ$group[x$data$tax_occ[, group1] > x$data$tax_occ[, group2]] < -
# x$data$tax_occ[, group1][x$data$tax_occ[, group1] >= x$data$tax_occ[, group2]]
# x$data$tax_occ$group[x$data$tax_occ[, group1] < x$data$tax_occ[, group2]] < -
# x$data$tax_occ[, group2][x$data$tax_occ[, group1] < x$data$tax_occ[, group2]] * -1
# }
#
#
# # x$data$tax_occ$total <- rowSums(x$data$taxon_counts[, -1]) # -1 = taxon_id column
#
# set.seed(1)  # This makes the plot appear the same each time it is run
# tree < - metacoder::
#     heat_tree(x,
#               node_label=taxon_names,
#               node_size=n_obs,
#               node_color=n_obs,
#               make_node_legend=T,
#               make_edge_legend=F,
#               # node_color_interval = c(-max(x$data$tax_occ[, 2:3])* 2.3e-2,
#               #                         max(x$data$tax_occ[, 2:3])* 2.3e-2),
#               node_color_range=c("#DB9261", "gray", "#5A99CC"),  # The color palette used
#               # node_size_axis_label = "OTU count",
#               # node_color_axis_label = "Cells with reads",
#               layout="davidson-harel",  # The primary layout algorithm
#               initial_layout="reingold-tilford"
#               )  # The layout algorithm that initializes node locations
#
# tree
# return (tree)
# }



#' UpSet Plot for Shared and Unique Microbial Features Across Groups.
#'
#' @param seurat.object seurat object
#' @param tax_level Default: \code{"species"}. Taxonomy level.
#' @param group.by Default: \code{"orig.ident"}.
#' Name of seurat object metadata columns to group (color) cells by.
#'
#' @param custom.color Default: \code{NULL}. which will use 20 built-in colors \code{c('#D97894','#7C98CE','#BBCB54', '#A7BED3',
#' '#C9A2E1', '#748977', '#BEBE9E', '#4CD28F',
#' '#CFAB84', '#E3E15C', '#DFC7CF', '#6BB8BF',
#' '#CBCBCB', '#9696eb', '#DED1AD', '#E29F58',
#' '#D8A1CE', '#B16455', '#ad9bf4', '#BEB68D')}.
#' Users could provide a vector of colors, each color corresponds to one group mapped to param \code{"group.by"}.
#'
#' @return upset plot
#' @export
#'
#' @examples
#' plotUpset(seurat.object)
#'
plotUpset <- function(seurat.object,
                      tax_level = "species",
                      group.by = "orig.ident",
                      custom.color = NULL) {
  check_package("UpSetR")
  if (is.null(custom.color)) {
    custom.color <- c(
      "#D97894", "#7C98CE", "#BBCB54", "#A7BED3", "#C9A2E1",
      "#748977", "#BEBE9E", "#4CD28F", "#CFAB84", "#E3E15C",
      "#DFC7CF", "#6BB8BF", "#CBCBCB", "#9696eb", "#DED1AD",
      "#E29F58", "#D8A1CE", "#B16455", "#ad9bf4", "#BEB68D",
      "#B6B6B6"
    )
  }
  if (is.null(seurat.object@misc$microbe$phyloseq.object)) {
    microbe.phyloseq.object <- seurat.object@misc$microbe$seurat.object %>%
      seurat2phyloseq()
  } else {
    microbe.phyloseq.object <- seurat.object@misc$microbe$phyloseq.object
  }
  microbe.phyloseq.object.tax_level <- microbe.phyloseq.object %>%
    phyloseq_select_level(tax_level = tax_level)

  otu_tab <- microbe.phyloseq.object.tax_level@otu_table %>%
    as.data.frame() %>%
    t() %>%
    stats::aggregate(by = list(seurat.object@meta.data[, group.by]), sum) %>%
    {
      rownames(.) <- .$Group.1
      .[["Group.1"]] <- NULL
      t(.)
    } %>%
    as.data.frame()

  upset_list <- list()
  for (clu in colnames(otu_tab)) {
    upset_list[[clu]] <- rownames(otu_tab)[otu_tab[, clu] != 0]
  }

  p_upset <- UpSetR::upset(UpSetR::fromList(upset_list),
    sets = colnames(otu_tab),
    sets.x.label = paste0(tax_level, " number"),
    mainbar.y.label = paste0("Intersection ", tax_level, " number"),
    order.by = "freq",
    keep.order = TRUE,
    mb.ratio = c(0.7, 0.3),
    sets.bar.color = custom.color[seq_along(upset_list)],
    main.bar.color = "#a9a9a9",
    set_size.show = TRUE
  )

  return(p_upset)
}



#' Trajectory Visualization for Host Transcriptomes and Microbial Communities.
#'
#' @param seurat.object seurat object
#' @param method Default: \code{"monocle3"}.
#' @param feature Default: \code{"pseudotime"}.
#' @param gene_name Default: \code{NULL}.
#' @param tax_name Default: \code{NULL}.
#' @param tax_level Default: \code{"species"}.
#' @param custom.color Default: \code{NULL}.
#' @param ... additional params from function \strong{plot_cells}.
#' Details params in \code{\link[monocle3]{plot_cells}}.
#'
#' @return Trajectory plot.
#' @export
#'
#' @examples
#' # data("seurat.object")
#' # seurat.object <- runSeurat(seurat.object)
#' # seurat.object <- runTrajectory(seurat.object)
#' # plotTrajectory(seurat.object, feature = "State")
#' # plotTrajectory(seurat.object, gene_name = "IGHM")
#' # plotTrajectory(seurat.object, tax_name = "Escherichia", tax_level = "genus")
#'
plotTrajectory <- function(seurat.object,
                           method = "monocle3",
                           feature = "pseudotime",
                           gene_name = NULL,
                           tax_name = NULL,
                           tax_level = "species",
                           custom.color = NULL,
                           ...) {
  if (is.null(seurat.object@misc$host$trajectory$monocle3)) {
    stop('Run "seurat.object <- runTrajectory(seurat.object)" first!')
  }
  cds <- seurat.object@misc$host$trajectory$monocle3
  if (is.null(gene_name) && is.null(tax_name)) {
    p_trajectory <- monocle3::plot_cells(cds,
      color_cells_by = feature,
      label_cell_groups = FALSE, ...
    )
  }
  if (!is.null(gene_name)) {
    p_trajectory <- monocle3::plot_cells(cds,
      genes = c(gene_name),
      label_cell_groups = FALSE, ...
    )
  }
  if (!is.null(tax_name)) {
    if (is.null(seurat.object@misc$microbe$phyloseq.object)) {
      microbe.phyloseq.object <- seurat.object@misc$microbe$seurat.object %>%
        seurat2phyloseq()
    } else {
      microbe.phyloseq.object <- seurat.object@misc$microbe$phyloseq.object
    }
    cds@colData$taxo_abu <- (phyloseq::prune_samples(
      SeuratObject::Cells(cds),
      microbe.phyloseq.object
    ) %>%
      phyloseq_select_tax(tax_name, tax_level)) + 1 %>% log2()
    p_trajectory <- monocle3::plot_cells(cds,
      color_cells_by = "taxo_abu",
      label_cell_groups = FALSE,
      label_branch_points = FALSE, ...
    ) +
      scale_color_gradientn(
        name = paste0(tax_level, " - ", tax_name),
        colours = c("blue", "white", "red")
      )
  }
  return(p_trajectory)
}


#' Visualize Correlations Between Microbial Features and Host Genes.
#'
#' @param seurat.object seurat object
#' @param microbe_list Default: \code{1:10}.
#' @param feature_list Default: \code{1:20}.
#'
#' @return  seurat object
#' @export
#'
#' @examples
#' NULL
#'
plotCorr <- function(seurat.object,
                     microbe_list = 1:10,
                     feature_list = 1:20) {
  if (is.vector(feature_list) && is.integer(feature_list)) {
    feature_list <- rownames(seurat.object@misc$corr$corr_result@correlation_matrix)[unique(feature_list)]
  } else if (is.vector(feature_list) && is.character(feature_list)) {
    feature_list <- intersect(rownames(seurat.object@misc$corr$corr_result@correlation_matrix), unique(feature_list))
  }

  if (is.vector(microbe_list) && is.integer(microbe_list)) {
    microbe_list <- colnames(seurat.object@misc$corr$corr_result@correlation_matrix)[unique(microbe_list)]
  } else if (is.vector(microbe_list) && is.character(microbe_list)) {
    microbe_list <- intersect(colnames(seurat.object@misc$corr$corr_result@correlation_matrix), unique(microbe_list))
  }


  correlation_matrix_selected <- seurat.object@misc$corr$corr_result@correlation_matrix[feature_list, microbe_list] %>% as.matrix()
  p_values_matrix_selected <- seurat.object@misc$corr$corr_result@p_values_matrix[feature_list, microbe_list] %>% as.matrix()
  coul <- c(
    "#00008B", "#1A1A97", "#3535A3", "#5050AF",
    "#6B6BBB", "#8686C8", "#A1A1D4", "#BBBBE0",
    "#D6D6EC", "#F1F1F8", "#F8F1F1", "#ECD6D6",
    "#E0BBBB", "#D4A1A1", "#C88686", "#BB6B6B",
    "#AF5050", "#A33535", "#971A1A", "#8B0000"
  )
  check_package("corrplot")
  corrplot::corrplot(correlation_matrix_selected,
    col = coul,
    tl.cex = 0.7,
    tl.col = "black",
    tl.srt = 45,
    tl.offset = 0.5,
    cl.pos = "b",
    cl.align.text = "l",
    cl.length = 5,
    cl.ratio = 0.1,
    cl.cex = 0.7,
    addgrid.col = "white",
    method = "color",
    p.mat = p_values_matrix_selected,
    insig = "label_sig",
    sig.level = c(0.001, 0.01, 0.05),
    pch.cex = 1.2,
    is.corr = FALSE,
    mar = c(0, 0, 2, 2),
    xpd = T
  )
}


if (getRversion() >= "2.15.1") utils::globalVariables(c("variable", "Gene", "value", "cluster"))
#' Heatmap Visualization of Differentially Expressed Marker Genes.
#'
#' @param seurat.object seurat object
#' @param dtype Default: \code{"host.marker"}.
#' @param row.group.by Default: \code{"seurat_clusters"}.
#' @param col.cluster.by Default: \code{"seurat_clusters"}.
#' @param feature.row Default: \code{3}.
#'
#' @return heatmap plot
#' @export
#'
#' @examples #
plotHeatmap <- function(seurat.object,
                        dtype = "host.marker",
                        row.group.by = "seurat_clusters",
                        col.cluster.by = "seurat_clusters",
                        feature.row = 3) {
  check_package("patchwork")
  check_package("tibble")
  if (dtype == "host.marker") {
    if (is.null(seurat.object@misc$host$markers)) {
      stop("Slot '\033[1;31mseurat.object@misc$host$markers\033[0m' not found! Run 'runSeurat' function (which include 'Seurat::FindAllMarkers' process) or run 'Seurat::FindAllMarkers' first.")
    } else {
      markers <- seurat.object@misc$host$markers
      gene <- data.frame()
      clu <- 1
      for (clu in unique(markers$cluster)) {
        marker.clu <- markers[markers$cluster == clu, ]
        marker.clu <- marker.clu[!marker.clu$gene %in% gene$gene, ]
        if (nrow(marker.clu) > 15) {
          gene <- rbind(gene, marker.clu[1:15, c("gene", "cluster")])
        } else {
          gene <- rbind(gene, marker.clu[, c("gene", "cluster")])
        }
      }
      rownames(gene) <- gene$gene

      clu <- gene$cluster[1]
      gene$x <- 0
      gene$y <- 0
      x <- 1
      y <- 1

      for (g in seq_along(gene$gene)) {
        if (gene$cluster[g] == clu) {
          if (x > 3) {
            x <- 1
            y <- y + 1
          }
        } else {
          clu <- gene$cluster[g]
          x <- 1
          y <- y + 1
        }

        gene$x[g] <- x
        gene$y[g] <- y
        x <- x + 1
      }

      gene$y <- max(gene$y) - gene$y + min(gene$y)

      # Ensure that gene$gene is a factor with the correct levels
      gene$gene <- factor(gene$gene, levels = gene$gene)

      # Calculate average expression
      vag_exp <- Seurat::AggregateExpression(seurat.object,
        assays = "RNA",
        features = gene$gene,
        group.by = "singleR_cluster",
        slot = "data"
      )

      # Get unique cell types
      celltype <- unique(Seurat::FetchData(seurat.object, vars = c("singleR_cluster")))

      # Set celltype clusters as factors
      celltype$singleR_cluster <- factor(celltype$singleR_cluster,
        levels = unique(celltype$singleR_cluster)
      )

      # Remove duplicate genes
      gene <- gene[!duplicated(gene$gene), ]

      # Normalize the expression data (z-score normalization)
      dat <- t(apply(vag_exp$RNA, 1, function(x) (x - mean(x)) / stats::sd(x)))
      dat <- as.data.frame(dat)

      # Subset dat to the order of genes in gene$gene
      dat <- dat[gene$gene %>% as.character(), ]

      # dat.sum <- dat %>% stats::aggregate(by = list(as.character(gene$cluster)), sum) %>%
      #   {
      #     rownames(.) <- .$Group.1
      #     .[["Group.1"]] <- NULL
      #     t(.)
      #   } %>% t()
      #
      # dat.sum <- dat.sum %>% as.data.frame() %>%
      #   rownames_to_column('Cluster') %>%
      #   reshape2::melt()

      # Create a list to store the rows for dat.new
      dat.list <- list()

      # Loop through unique clusters
      for (clu in unique(gene$cluster)) {
        # Extract data for the current cluster
        dat.clu <- dat[gene$gene[gene$cluster == clu], ]

        # Check if the number of rows is divisible by 3
        if (nrow(dat.clu) %% 3 == 0) {
          dat.list[[length(dat.list) + 1]] <- dat.clu
        } else {
          # Calculate the number of rows to add (padding)
          pad_rows <- 3 - nrow(dat.clu) %% 3
          matrix_na <- matrix(NA, nrow = pad_rows, ncol = ncol(dat))

          # Add NA rows with appropriate column names
          colnames(matrix_na) <- colnames(dat)

          # Append the data for the current cluster with padding
          dat.list[[length(dat.list) + 1]] <- rbind(dat.clu, matrix_na)
        }
      }

      # Combine all the rows into a single data frame
      dat.new <- do.call(rbind, dat.list)
      # Final assignment to dat
      dat <- dat.new
      dat <- dat %>%
        tibble::rownames_to_column("Gene") %>%
        reshape2::melt()

      dat$Gene <- factor(dat$Gene, levels = rev(unique(dat$Gene)))
      dat$variable <- factor(dat$variable, levels = unique(dat$variable))
      dat <- dat[order(dat$Gene), ]



      heatmap_color <- c(
        "#67001F", "#B2182B", "#D6604D", "#F4A582",
        "#FDDBC7", "#F7F7F7", "#D1E5F0", "#92C5DE",
        "#4393C3", "#2166AC", "#053061"
      )

      pal <- rev(grDevices::colorRampPalette(heatmap_color)(500))

      label <- levels(dat$variable)

      dat[dat$value > 2 & !is.na(dat$value), "value"] <- 2
      dat[dat$value < -2 & !is.na(dat$value), "value"] <- -2


      custom.color <- c(
        "#D97894", "#7C98CE", "#BBCB54", "#A7BED3", "#C9A2E1",
        "#748977", "#BEBE9E", "#4CD28F", "#CFAB84", "#E3E15C",
        "#DFC7CF", "#6BB8BF", "#CBCBCB", "#9696eb", "#DED1AD",
        "#E29F58", "#D8A1CE", "#B16455", "#ad9bf4", "#BEB68D"
      )

      p1 <- ggplot(dat, aes(as.numeric(variable), Gene, fill = value)) +
        geom_tile() +
        scale_y_discrete(expand = c(0, 0)) +
        scale_x_continuous(
          expand = c(0, 0),
          breaks = seq(1, length(levels(celltype$singleR_cluster)), by = 1),
          labels = label,
          # sec.axis = dup_axis(
          #   breaks = seq(1, length(levels(celltype$singleR_cluster)), by = 1),
          #   labels = label)
        ) +
        scale_fill_gradientn(colors = pal, limits = c(-2, 2), name = "Z Score", na.value = "#c0c0c0") +
        geom_hline(yintercept = as.numeric(cumsum(rev(ceiling(table(gene$cluster)[-1] / 3) * 3))) + .5, color = "white") +
        geom_vline(xintercept = as.numeric(cumsum(table(celltype$singleR_cluster)) + .5), linetype = 2) +
        theme(
          text = element_text(face = "bold"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = .5),
          axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = .5)
        )

      p2 <- ggplot(gene, aes(x, y, fill = cluster)) +
        geom_tile() +
        geom_text(aes(label = gene), family = "serif", fontface = "italic") +
        scale_y_continuous(expand = c(0, 0)) +
        scale_fill_manual(values = custom.color) +
        theme(
          text = element_text(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank()
        ) +
        scale_x_continuous(expand = c(0, 0)) +
        geom_hline(yintercept = as.numeric(cumsum(rev(ceiling(table(gene$cluster)[-1] / 3) * 3) / 3)) + .5, color = "white") +
        guides(
          fill = guide_legend(
            ncol = 2,
            byrow = TRUE
          )
        )
      pic.heatmap <- p2 + p1 + patchwork::plot_layout(ncol = 2, widths = c(3, 3), guides = "collect")
      p <- pic.heatmap & theme(plot.margin = margin(0, 0, 0, 0))
      return(p)
    }
  }
}
