\name{analyzeNAs}
\alias{analyzeNAs}
\title{Analyze the NA distribution}
\usage{
  analyzeNAs(cancer, array)
}
\description{
  Analyze the NA distribution in the filtered TCGA samples,
  per probe and sample (uses the files from
  filterBarcodes). Also, get the standard deviation for
  every probe.
}
\details{
  It will write two text files to disk, one containing the
  number of NAs and standard deviation of every probe and
  the other containing the number of NAs for every sample.
}

