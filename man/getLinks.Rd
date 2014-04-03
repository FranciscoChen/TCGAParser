\name{getLinks}
\alias{getLinks}
\title{Get matching links in a URL}
\usage{
  getLinks(URL, REGEXPR)
}
\arguments{
  \item{URL}{Link to a website to scrape.}

  \item{REGEXPR}{Regular expression to match to the scraped
  contents.}
}
\value{
  Returns a list[[i]][[1|2]] with i matches and two
  sublevels for every match: 1 is the link and 2 is the
  last modified date (if any). If the given URL does not
  exist or there are no matches, returns nothing.
}
\description{
  Given a URL and a regular expression, returns all links
  in the URL matching the regular expression along its last
  modified date.
}
\details{
  Requires RCurl package to getURL from https.
}

