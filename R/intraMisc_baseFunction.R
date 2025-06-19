#' Convert mpa format object to taxonomy table
#'
#' @param taxonomy_list String format vector
#'
#' @importFrom stringr str_to_title
#'
#' @return taxonomy table
#'
#' @export
#'
#'
mpa2tax <- function(taxonomy_list) {
  taxonomy_list <- c(taxonomy_list)
  tax_df <- data.frame(matrix(ncol = 7, nrow = length(taxonomy_list)))
  colnames(tax_df) <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

  # Split each tax string and fill the data frame tax_df.
  for (i in seq_along(taxonomy_list)) {
    tax_levels <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
    tax_list <- unlist(strsplit(taxonomy_list[i], "\\|"))
    for (level in tax_list) {
      prefix <- substr(level, 1, 2)
      value <- substr(level, 4, nchar(level))
      column_name <- switch(prefix,
        "k_" = "kingdom",
        "p_" = "phylum",
        "c_" = "class",
        "o_" = "order",
        "f_" = "family",
        "g_" = "genus",
        "s_" = "species",
        ""
      )
      tax_df[i, column_name] <- value
    }
  }
  for (i in 7:2) {
    count <- 1
    for (t in unique(tax_df[, tax_levels[i]])) {
      if (is.na(unique(tax_df[tax_df[, tax_levels[i]] == t, tax_levels[i - 1]]))) {
        tax_df[tax_df[, tax_levels[i]] == t, tax_levels[i - 1]] <- paste0("Unknown", str_to_title(tax_levels[i - 1]), count)
        count <- count + 1
      }
    }
  }
  rownames(tax_df) <- tax_df$species
  return(tax_df)
}


#' Convert mpa format abundance matrix to phyloseq object
#'
#' @param microbe.counts mpa format abundance matrix
#'
#' @return phyloseq object
#' @export
#'
#'
constructPhyloseq <- function(microbe.counts) {
  taxTab <- mpa2tax(rownames(microbe.counts))

  rownames(microbe.counts) <- taxTab$species

  phyloseq.object <- phyloseq::phyloseq(
    phyloseq::otu_table(microbe.counts, taxa_are_rows = TRUE),
    phyloseq::tax_table(taxTab %>% as.matrix())
  )

  phyloseq::sample_data(phyloseq.object) <- data.frame(
    row.names = phyloseq.object@otu_table@.Data %>% colnames(),
    orig.ident = (phyloseq.object@otu_table@.Data %>%
      colnames() %>%
      stringr::str_split("_", simplify = TRUE))[, 1],
    nCount_Microbe = colSums(phyloseq.object@otu_table@.Data),
    nFeature_Microbe = colSums(phyloseq.object@otu_table@.Data != 0)
  )
  return(phyloseq.object)
}




#' Select Microbial Features by Taxonomic Level from a phyloseq Object.
#'
#' @param microbe.phyloseq.object phyloseq object
#' @param tax_level Default: \code{"genus"}
#'
#' @return microbe phyloseq object
#' @export
#'
phyloseq_select_level <- function(microbe.phyloseq.object,
                                  tax_level = "genus") {
  if (tax_level == "species") {
    return(microbe.phyloseq.object)
  }
  # tax_levels <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  # tax.tab <- phyloseq::tax_table(microbe.phyloseq.object) %>% as.data.frame()
  # otu.tab <- phyloseq::otu_table(microbe.phyloseq.object) %>% as.data.frame()
  #
  # otu.tab <- otu.tab %>% stats::aggregate(by = list(tax.tab[, tax_level]), sum)
  #
  # rownames(otu.tab) <- otu.tab$Group.1
  # otu.tab$Group.1 <- NULL
  #
  #
  #
  # for (i in c(which(tax_levels == tax_level) + 1):length(tax_levels)) {
  #   tax.tab[, tax_levels[i]] <- NULL
  # }
  # tax.tab <- tax.tab %>% dplyr::filter(!duplicated(.))
  # rownames(tax.tab) <- tax.tab[, tax_level]
  # tax.tab <- tax.tab[rownames(otu.tab), ]
  #
  # microbe.phyloseq.object.tax_level <- phyloseq::phyloseq(
  #   phyloseq::otu_table(otu.tab, taxa_are_rows = T),
  #   phyloseq::tax_table(tax.tab %>% as.matrix()),
  #   phyloseq::sample_data(microbe.phyloseq.object@sam_data)
  # )
  microbe.phyloseq.object.tax_level <- microbe.phyloseq.object %>% phyloseq::tax_glom(taxrank = tax_level)
  rownames(microbe.phyloseq.object.tax_level@otu_table@.Data) <- microbe.phyloseq.object.tax_level@tax_table@.Data %>%
    as.data.frame() %>%
    .[, tax_level]
  rownames(microbe.phyloseq.object.tax_level@tax_table@.Data) <- microbe.phyloseq.object.tax_level@tax_table@.Data %>%
    as.data.frame() %>%
    .[, tax_level]

  return(microbe.phyloseq.object.tax_level)
}

# select abundance of tax_name at tax_level level.
phyloseq_select_tax <- function(phyloseq.object, tax_name = "Bacillus", tax_level = "genus") {
  tax_tab <- phyloseq.object@tax_table %>% as.data.frame()
  selected_tax_species <- tax_tab[tax_tab[, tax_level] == tax_name, ]$species
  return(colSums(phyloseq::prune_taxa(selected_tax_species, phyloseq.object)@otu_table))
}

# select abundance of tax_name at tax_level level.
#' select abundance of tax_name at tax_level level.
#'
#' @param phyloseq.object phyloseq object
#' @param tax_name Default: \code{"Bacillus"}
#' @param tax_level Default: \code{"genus"}
#'
#' @return phyloseq object at selected tax level.
#' @export
#'
#' @examples
#' NULL
phyloseq_select_tax <- function(phyloseq.object,
                                tax_name = "Bacillus",
                                tax_level = "genus") {
  tax_tab <- phyloseq.object@tax_table %>% as.data.frame()
  selected_tax_species <- tax_tab[tax_tab[, tax_level] == tax_name, ]$species
  return(colSums(phyloseq::prune_taxa(selected_tax_species, phyloseq.object)@otu_table))
}


# Convert microbe seurat object to microbe phyloseq object.
#' Convert microbe seurat object to microbe phyloseq object.
#'
#' @param microbe.seurat.object Microbe seurat object
#'
#' @return Microbe phyloseq object
#' @export
#'
#' @examples
#' NULL
#'
seurat2phyloseq <- function(microbe.seurat.object) {
  if (microbe.seurat.object@active.assay == "Microbe") {
    microbe.phyloseq.object <- microbe.seurat.object %>%
      SeuratObject::JoinLayers() %>%
      {
        as.data.frame(as.matrix(.@assays$Microbe@layers$counts))
      } %>%
      {
        rownames(.) <- microbe.seurat.object %>%
          SeuratObject::Features() %>%
          stringr::str_replace_all("#", "\\|") %>%
          stringr::str_replace_all("-", "_")
        colnames(.) <- microbe.seurat.object %>% SeuratObject::Cells()
        .
      } %>%
      constructPhyloseq()
    return(microbe.phyloseq.object)
  } else {
    stop("ERROR: Not microbe seurat object!")
  }
}

# check suggest packages
check_package <- function(pack_name) {
  if (!requireNamespace(pack_name, quietly = TRUE)) {
    stop(paste0("The package {", pack_name, "} is required for this function"))
  }
}

# rewrite the sbuset method
#' reconstructed sbuset method
#'
#' @param object S4 object.
#' @param ... params of \code{subset} function from \strong{\code{SummarizedExperiment}}
#'
#' @return S4 object.
#' @export
#'
#' @examples
#' NULL
#'
subset <- function(object, ...) {
  check_package("SummarizedExperiment")
  object.subset <- SummarizedExperiment::subset(object, ...)
  if (typeof(object) == "S4") {
    if (!is.null(object@misc$microbe$seurat.object)) {
      cells <- SeuratObject::Cells(object.subset)
      object.subset@misc$microbe$seurat.object <- object@misc$microbe$seurat.object %>%
        SummarizedExperiment::subset(cells = cells)
      object.subset$nCount_Microbe <- object.subset@misc$microbe$seurat.object$nCount_Microbe
      object.subset$nFeature_Microbe <- object.subset@misc$microbe$seurat.object$nFeature_Microbe
    }
    if (!is.null(object@misc$microbe$phyloseq.object)) {
      cells <- SeuratObject::Cells(object.subset)
      object.subset@misc$microbe$phyloseq.object <- object@misc$microbe$phyloseq.object %>%
        phyloseq::prune_samples(cells, .)
      object.subset$nCount_Microbe <- colSums(object.subset@misc$microbe$phyloseq.object@otu_table@.Data)
      object.subset$nFeature_Microbe <- colSums(object.subset@misc$microbe$phyloseq.object@otu_table@.Data != 0)
    }
  }
  return(object.subset)
}
