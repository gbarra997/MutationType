#' General Statistics on Mutations
#'
#' This function takes a data frame of mutations and calculates the frequency of
#' each mutation type (only two types in this case). It can optionally generate a bar plot of the mutation
#' frequencies.
#'
#' @param mutations A data frame containing mutation information.
#'                  The data frame should be the result of the function seekMutAndSeq have a 'ref' column.
#' @param plotMutation Logical. If TRUE, a bar plot of mutation frequencies is
#'                      generated. Default is FALSE.
#'
#' @return A list with two elements:
#'         - Element 1: A dataframe with the mutation frequencies.
#'         - Element 2: If plotMutation is TRUE, a bar plot of mutation frequencies.
#'
#' @details
#' Through this function one can gets the 'ref' column from the input dataframe and calculates
#' the frequency of each mutation type. If plotMutation is TRUE, it generates a
#' bar plot using the ggplot2 package, displaying the mutation frequencies with
#' specific labels. The output is a list containing the table and the plot
#' @import ggplot2
#' @import methods
#' @examples
#' df <- data.frame(
#'   id = c("rs7410291", "rs1234567"),
#'   chr = c("chr22", "chrX"),
#'   pos = c("50300078", "100200300"),
#'   start = c("0", "0"),
#'   end = c("0", "0"),
#'   ref = c("A", "C"),
#'   alt = c("G", "T"),
#'   sequences = c("CCGCGCCCT", "TCCCGGGGA"),
#'   mutation_type = c("C[T>C]C", "C[C>T]G")
#' )
#' generalStat(df)
#' @export
generalStat <- function(mutations, plotMutation = FALSE) {
  # Takes in input the data frame
  # One column is enough to get the freq. of the mutations
  refMutation <- mutations$ref
  # Get the freq of the two bases
  countTable <- as.data.frame(table(refMutation))
  colnames(countTable) <- c("refMutation", "Freq")
  # Make a bar plot and set specific labels
  g <- ggplot(mutations, aes(x = refMutation, fill = refMutation)) +
    geom_bar() +
    scale_x_discrete(labels = c("[C>*]", "[T>*]")) +
    labs(x = "", fill = "Mutations Type") +
    theme(legend.position = "none")

  if (plotMutation == TRUE) {
    show(g)
  }
  mutSeq <- mutations$mutation_type
  countTable1 <- as.data.frame(table(mutSeq))
  colnames(countTable1) <- c("Type", "Freq")
  g1 <- ggplot(countTable1, aes(x = countTable1$Type, y = countTable1$Freq)) +
    geom_bar(stat = "identity", position = "identity") +
    labs(x = "Mutation Type", y = "Frequency", fill = "blue") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  if (plotMutation == TRUE) {
    show(g1)
  }
  # Return a list where the first element is the countTable with the freq. while the second is the plot
  return(list(countTable, g,countTable1, g1))
}




