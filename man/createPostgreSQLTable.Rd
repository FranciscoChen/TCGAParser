\name{createPostgreSQLTable}
\alias{createPostgreSQLTable}
\title{Create PostgreSQL table}
\usage{
  createPostgreSQLTable(cancer, array, drv, ...)
}
\description{
  Create the table in the PostgreSQL database (needs the
  sampleinfo file from filterBarcodes).
}
\details{
  Requires RPostgreSQL package and a connection to a
  PostgreSQL database.  drv A character string specifying
  the database management system driver.  ...  Arguments
  needed to connect to the database, such as user,
  password, dbname, host, port, etc.
}

