\name{corFromTableToTable}
\alias{corFromTableToTable}
\title{Do PostgreSQL table correlations}
\usage{
  corFromTableToTable(drv, ..., from.query = NULL,
    and.query = NULL, to.table = NULL, stdev.threshold.from
    = 0, stdev.threshold.and = 0, pval.threshold = 1,
    nthreads = 1)
}
\description{
  Do correlations of tables stored in PostgreSQL and store
  significant correlations in another PostrgreSQL table.
}
\details{
  Requires parallel package.  drv A character string
  specifying the database management system driver.  ...
  Arguments needed to connect to the database, such as
  user, password, dbname, host, port, etc.  from.query
  PostgreSQL statement to select or filter the data in a
  specified PostgreSQL table for correlations.  and.query
  PostgreSQL statement to select or filter the data in a
  second specified PostgreSQL table for correlations.
  to.table PostgreSQL table where the significant
  correlations will be written.  stdev.threshold.from
  Filters out all rows containing values with lower
  standard deviation than the specified threshold (in
  from.query).  stdev.threshold.and Filters out all rows
  containing values with lower standard deviation than the
  specified threshold (in and.query).  pval.threshold
  Correlations with a p-value higher than this threshold
  will be considered significant and registered in
  to.table.  nthreads Number of threads this functions will
  use.
}

