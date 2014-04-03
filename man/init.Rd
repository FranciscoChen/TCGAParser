\name{init}
\alias{init}
\title{Initialise variables}
\usage{
  init(cancer, array)
}
\arguments{
  \item{cancer}{The cancer type we intend to work with,
  abbreviated as in the TCGA.}

  \item{array}{Array from which the data was processed.}
}
\description{
  Declares and assigns variables like URLs and regular
  expressions.
}
\seealso{
  https://tcga-data.nci.nih.gov/datareports/codeTablesReport.htm
  Currently fully supported arrays: humanmethylation450,
  illuminahiseq_rnaseqv2.
}

