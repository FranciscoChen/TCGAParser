\name{getDataframe}
\alias{getDataframe}
\title{Get a data frame from a URL}
\usage{
  getDataframe(URL, header = TRUE, sep = "\t", quote =
    "\"'", dec = ".", row.names, col.names, as.is =
    !stringsAsFactors, na.strings = "NA", colClasses = NA,
    nrows = -1, skip = 0, check.names = TRUE, fill =
    !blank.lines.skip, strip.white = FALSE,
    blank.lines.skip = TRUE, comment.char = "#",
    allowEscapes = FALSE, flush = FALSE, stringsAsFactors =
    default.stringsAsFactors(), fileEncoding = "", encoding
    = "unknown", text)
}
\arguments{
  \item{URL}{Link to a text file.}
}
\value{
  Returns a data frame or nothing.
}
\description{
  Download a file as a temporary .txt file in the working
  directory and read it into R as a data frame.
}
\details{
  The file is read with tab as separators and expecting a
  header. All other arguments are as in read.table.
}
\seealso{
  read.table
}

