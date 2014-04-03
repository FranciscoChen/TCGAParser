\name{fixDataframe}
\alias{fixDataframe}
\title{Fix a data frame}
\usage{
  fixDataframe(dataframe, array)
}
\arguments{
  \item{dataframe}{Data frame to fix.}

  \item{array}{Array from which the data was processed.
  Currently supported arrays: humanmethylation450,
  illuminahiseq_rnaseqv2.}
}
\value{
  Returns the fixed data frame or nothing.
}
\description{
  Format a data frame to a common standard. Varies from
  array to array.
}

