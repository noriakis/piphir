% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/profileMetagenome.R
\name{profileMetagenome}
\alias{profileMetagenome}
\title{profileMetagenome}
\usage{
profileMetagenome(
  taxTable,
  copyNumTable,
  KOTable,
  blastRes,
  full = FALSE,
  fullOutput = "file",
  fullTemp = "temporary_full",
  fullID = NULL,
  forceCalc = TRUE
)
}
\arguments{
\item{taxTable}{taxonomy abundance table}

\item{copyNumTable}{16S copy number table (row.names corresponds to reference tax ID)}

\item{KOTable}{precalculated KO table}

\item{blastRes}{output of alignSequences function (vsearch)
The second column must contain representative taxon ID linked to KO table.}

\item{full}{return the full output, default to FALSE}

\item{fullOutput}{"file" or "inmemory"}

\item{fullTemp}{if fullOutput is "file", {{fileTemp}}.txt will be created.}

\item{fullID}{subset to this ID in full mode, when the data size is large.}

\item{forceCalc}{force calculation when the index of normalized hit table and
KO copy number table does not match, by removing blast hits that not in CN table.}
}
\value{
KO profile table
}
\description{
the main workflow
}
