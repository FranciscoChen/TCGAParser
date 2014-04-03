\name{pairBarcodes}
\alias{pairBarcodes}
\title{Pair TCGA barcodes}
\usage{
  pairBarcodes (cancer, array1, array2)
}
\description{
  Pairs barcodes between two arrays for a specific cancer
  type.
}
\details{
  Filters for the first 16 TCGA barcode characters, which
  refer to participant, sample type and vial. Duplicates
  are removed, along with normal and control tissue
  samples. Finally, it will write into disk the data frames
  for each array.
}

