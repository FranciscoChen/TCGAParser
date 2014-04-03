\name{generatePairingTable}
\alias{generatePairingTable}
\title{Generate a TCGA barcode pairing table}
\usage{
  generatePairingTable()
}
\description{
  For each tumor's humanmethylation450/ get the number of
  matching barcodes in agilentg4502a_07_3/,
  illuminaga_rnaseq/, illuminaga_rnaseqv2/,
  illuminahiseq_rnaseq/ and illuminahiseq_rnaseqv2/.
}
\details{
  Filters for the first 16 TCGA barcode characters, which
  refer to participant, sample type and vial. Duplicates
  are removed, along with normal and control tissue
  samples. Finally, it will write the result in the
  specified path.
}

