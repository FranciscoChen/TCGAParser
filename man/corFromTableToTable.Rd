\name{corFromTableToTable}
\alias{corFromTableToTable}
\title{Do PostgreSQL table correlations}
\usage{
  corFromTableToTable(drv, ..., from.schema = NULL,
    from.table = NULL, from.join = NULL, from.condition =
    NULL, and.schema = NULL, and.table = NULL, and.join =
    NULL, and.condition = NULL, to.schema = NULL, to.table
    = NULL, stdev.threshold.from = 0, stdev.threshold.and =
    0, pval.threshold = 1, corr.threshold = 0.50, nthreads
    = 1)
}
\arguments{
  \item{drv}{A character string specifying the database
  management system driver, "PostgreSQL".}

  \item{...}{Arguments needed to connect to the database,
  such as user, password, dbname, host, port, etc.}

  \item{from.schema}{PostgreSQL schema where the first
  table is.}

  \item{from.table}{PostgreSQL input table name.}

  \item{from.join}{Optional: a complementary PostgreSQL
  schema.table for the input.}

  \item{from.condition}{Everything that comes after the
  PostgreSQL WHERE (at input).}

  \item{and.schema}{Optional: a second PostgreSQL schema to
  correlate with.}

  \item{and.table}{Optional: a second PostgreSQL table
  name.}

  \item{and.join}{Optional: a complementary PostgreSQL
  schema.table for the second input.}

  \item{and.condition}{Everything that comes after the
  PostgreSQL WHERE (for the second input).}

  \item{to.schema}{PostgreSQL schema that contains the
  output table.}

  \item{to.table}{PostgreSQL output table name.}

  \item{stdev.threshold.from}{Filters out all rows
  containing values with lower standard deviation than the
  specified threshold (in the first input).}

  \item{stdev.threshold.and}{Filters out all rows
  containing values with lower standard deviation than the
  specified threshold (in the second input).}

  \item{pval.threshold}{Correlations with a p-value equal
  or lower than this threshold will be registered in the
  output.}

  \item{corr.threshold}{Correlations with an absolute
  correlation value equal or higher than this threshold
  will be registered in the output.}

  \item{nthreads}{Number of threads this function will use.
  Allows hyperthreading.}
}
\description{
  Do correlations of tables stored in PostgreSQL and store
  significant correlations in another PostrgreSQL table.
}
\details{
  Requires parallel package.
}
\seealso{
  http://www.postgresql.org/docs/
}

