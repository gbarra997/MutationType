% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generalStat.R
\name{generalStat}
\alias{generalStat}
\title{General Statistics on Mutations}
\usage{
generalStat(mutations, plotMutation = FALSE)
}
\arguments{
\item{mutations}{A data frame containing mutation information.
The data frame should be the result of the function seekMutAndSeq have a 'ref' column.}

\item{plotMutation}{Logical. If TRUE, a bar plot of mutation frequencies is
generated. Default is FALSE.}
}
\value{
A list with two elements:
- Element 1: A dataframe with the mutation frequencies.
- Element 2: If plotMutation is TRUE, a bar plot of mutation frequencies.
}
\description{
This function takes a data frame of mutations and calculates the frequency of
each mutation type (only two types in this case). It can optionally generate a bar plot of the mutation
frequencies.
}
\details{
Through this function one can gets the 'ref' column from the input dataframe and calculates
the frequency of each mutation type. If plotMutation is TRUE, it generates a
bar plot using the ggplot2 package, displaying the mutation frequencies with
specific labels. The output is a list containing the table and the plot
}
\examples{
df <- data.frame(
  id = c("rs7410291", "rs1234567"),
  chr = c("chr22", "chrX"),
  pos = c("50300078", "100200300"),
  start = c("0", "0"),
  end = c("0", "0"),
  ref = c("A", "C"),
  alt = c("G", "T"),
  sequences = c("CCGCGCCCT", "TCCCGGGGA"),
  mutation_type = c("C[T>C]C", "C[C>T]G")
)
generalStat(df)
}
