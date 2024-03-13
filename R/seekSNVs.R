#' Seek Single Nucleotide Variants (SNVs) from a VCF file
#'
#' This function extracts SNVs from a given VCF file and returns a dataframe
#' with necessary information for further analysis.
#'
#' @param vcfFile Path to the VCF file or an object of class "CollapsedVCF".
#' @param referenceGenome Character string specifying the reference genome.
#'
#' @return A dataframe with SNV information, including ID, chromosome, position,
#'         reference base, alternative base, and strand information.
#'
#' @details
#' The function first checks the input parameters and then processes the VCF file
#' to extract SNVs. SNVs are defined as variants where both the reference and
#' alternative alleles have a length of 1 according to the Variant Classification
#' definition (https://genome.sph.umich.edu/wiki/Variant_classification).
#' @import VariantAnnotation
#' @import methods
#' @importFrom MatrixGenerics rowRanges
#' @importFrom GenomicRanges width
#' @examples
#' vcfdata <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
#' results <- seekSNVs(vcfdata, referenceGenome="hg19")
#'
#' @export


seekSNVs <- function(vcfFile, referenceGenome) {

# Check if reference_genome is a character
  if (!is.character(referenceGenome)) {
    stop("Error: The 'reference_genome' parameter must be a character.")
  }
# Check if vcffile is an object of class "CollapsedVCF"
  if (is(vcfFile, "CollapsedVCF")) {
#If so get the GRanges
    setOfMutations <- MatrixGenerics::rowRanges(vcfFile)
  } else {
#Otherwise Check if vcffile is a path to a file with extension vcf.gz or vcf
    if (grepl("\\.vcf\\.gz$|\\.vcf$", vcfFile)) {
#If so read the vcf file and get again the GRanges
      vcfFile <- VariantAnnotation::readVcf(vcfFile, referenceGenome)
      setOfMutations <- MatrixGenerics::rowRanges(vcfFile)
    } else {
#Stop the exectution
      stop("The file does not have a valid extension (vcf.gz or vcf) or is not a collapsed VCF file of VariantAnnotation package.")
    }
  }
#Find only SNVs according to the definition: #https://genome.sph.umich.edu/wiki/Variant_classification [lenght of REF and ALT == 1]
  refWidths <- GenomicRanges::width(setOfMutations$REF)
  SNV <- subset(setOfMutations, GenomicRanges::width(setOfMutations$REF) == 1)
  SNV <- subset(SNV, width(unlist(SNV$ALT)) == 1)

#Create a dataframe with the SNVs general necessary information for next steps
  mutation_df <- data.frame(
    id = as.vector(SNV@ranges@NAMES),
    chr = as.character(paste("chr", SNV@seqnames@values, sep = "")),
    pos = as.character(SNV@ranges@start),
    start = as.character("0"),
    end = as.character("0"),
    ref = as.character(SNV$REF),
    alt = as.character(unlist(SNV$ALT)),
    stra = as.character(SNV@strand)
  )
  return(mutation_df)
}


