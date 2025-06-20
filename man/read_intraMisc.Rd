% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/intraMisc_readData.R
\name{read_intraMisc}
\alias{read_intraMisc}
\title{Load and Integrate Host Transcriptome and Intracellular Microbiome}
\usage{
read_intraMisc(
  named.outs.dirs,
  data.type = "sc",
  align.method = "kraken2",
  max.dim = 30000
)
}
\arguments{
\item{named.outs.dirs}{String or string vector, representing the paths to the outs directories of the output results from Cell Ranger or Spacer Ranger.}

\item{data.type}{Default: \code{"sc"}. Selected from \code{"sc"} and \code{"spt"}, sc represents single cell transcriptome data;
spt represents spatial transcriptome sequencing data.}

\item{align.method}{Default: \code{"kraken2"}. intraAlignment alignment methods, selected from ["kraken2", "blast"]}

\item{max.dim}{Default: \code{30000}.
Return phyloseq object if cell numbers <= \strong{\code{max.dim}}, else return seurat object.}
}
\value{
seurat object with microbe data. "seurat.object@misc$microbe"
}
\description{
Load and Integrate Host Transcriptome and Intracellular Microbiome
}
\examples{
named.outs.dirs <- c(
  system.file("extdata",
    "sc/SAMN20336351/outs/",
    package = "intraMisc"
  ),
  system.file("extdata",
    "sc/SAMN20336352/outs/",
    package = "intraMisc"
  )
)
names(named.outs.dirs) <- c("SAMN20336351", "SAMN20336352")
seurat.object <- read_intraMisc(named.outs.dirs)

}
