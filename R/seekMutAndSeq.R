#' Seek Mutations and Sequences around Single Nucleotide Variants (SNVs)
#'
#' This function takes a dataframe of SNVs and extracts sequences downward and upward
#' for each SNV based on the provided context length.
#' It replaces 'A's and 'G's in 'ref'
#' and 'alt' with 'T' and 'C', retrieves genomic ranges, and fetches sequences
#' from the reference genome. The resulting data frame includes mutation type
#' information and the extracted sequences.
#'
#' @param df A data frame containing SNV information, including 'chr', 'pos', 'ref',
#'           'alt', and 'stra' columns.
#' @param contextLength The length of the sequence around each SNV. It must be an
#'                     odd number greater than or equal to 3. Default is 3.
#' @param refGen Reference genome object.
#'
#' @return A data frame with information about mutations and the sequences around
#'         each SNV.
#'
#' @details
#' The function checks the validity of the input parameters and processes the
#' data frame to extract sequences around each SNV. It ensures that the sequences
#' do not go beyond the chromosome range. The resulting data frame includes
#' mutation type information, where 'A's and 'G's in 'ref' and 'alt' are replaced
#' with 'T' and 'C', respectively.
#'
#' @importFrom BSgenome getSeq
#' @importFrom GenomicRanges GRanges
#' @import GenomeInfoDb
#' @import VariantAnnotation
#' @import methods
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @examples
#' df <- data.frame(
#' id = c("rs7410291", "rs1234567"),
#' chr = c("chr22", "chrX"),
#' pos = c("50300078", "100200300"),
#' start = c("0", "0"),
#' end = c("0", "0"),
#' ref = c("A", "C"),
#' alt = c("G", "T"),
#' stra = c("*", "*")
#' )
#' refgenome <- BSgenome.Hsapiens.UCSC.hg38
#' seekMutAndSeq(df, contextLength=3, refgenome)
#'
#' @export


seekMutAndSeq <- function(df, contextLength = 3, refGen) {
  #Chek for the input number
  if (contextLength < 3 || contextLength %% 2 == 0) {
    #Stop the execution
    stop("Error: The contextLength must be an odd number greater than or equal to 3.")
  }
  #Get the length of the sequence that are DOWN and UP the SNV
  halfContext <- (contextLength - 1) / 2
  #Ensure the fact that it is not possible to go over the range (i.e. length over the chr)
  #Getting all the chromosome
  unique_chromosomes <- unique(df$chr)
  seq_lengths<- data.frame()
  #And checking that all the rows are within the range otherwise skip that row
  seq_lengths <- sapply(unique_chromosomes, function(chr) {
    chr_seq_length <- GenomeInfoDb::seqlengths(refGen)[chr]
    df_chr <- df[df$chr == chr, ]
    df_chr$start <- as.character(as.numeric(df_chr$pos) - halfContext)
    df_chr$end <- as.character(as.numeric(df_chr$pos) + halfContext)
    df_chr_filtered <- df_chr[df_chr$start >= 0 & df_chr$end <= chr_seq_length, ]
    return(df_chr_filtered)
  }, simplify = FALSE)
  #Combine the data frames for different chromosomes
  mutations <- do.call(rbind, seq_lengths)
  #Replace the bases A G in C T
  mutations$ref <- gsub("A", "T", mutations$ref)
  mutations$ref <- gsub("G", "C", mutations$ref)
  mutations$alt <- gsub("A", "T", mutations$alt)
  mutations$alt <- gsub("G", "C", mutations$alt)
  #Get the vector with the strings will be used to get the ranges
  genomicRanges <- apply(mutations, 1, function(row) {
    paste0(row["chr"], ":", row["start"], "-", row["end"], ":", row["stra"])
  })
  #Get the GRanges
  genomicLocation <- GRanges(genomicRanges, seqinfo = seqinfo(refGen))
  #Get the sequences  with the given coordinates
  sequences <- as.character(getSeq(refGen, genomicLocation))
  #Save the sequencies
  mutations$sequences <- sequences
  #Append a new column with the string that indicate how changed the base with the UP DOWN seq
  mutations$mutation_type <- mapply(function(seq, ref, alt) {
    left_part <- substr(seq, 1, halfContext)
    right_part <- substr(seq, nchar(seq) - halfContext + 1, nchar(seq))
    paste0(left_part, "[", ref, ">", alt, "]", right_part)
  }, mutations$sequences, mutations$ref, mutations$alt)
  return(mutations)
}

