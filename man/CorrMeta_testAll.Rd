\name{CorrMeta_testAll}
\alias{CorrMeta_testAll}
\title{
	Main function for fast eQTL analysis of correlated samples in CorrMeta package.
}
\description{
	\code{CorrMeta_testAll} function tests association of every gene and every SNP in the input data.
	The input data must be split in two so that the samples in each group are independent.
	The first data set is defined by (\code{snps1}, \code{gene1}, and \code{cvrt1}) and the second one by (\code{snps2}, \code{gene2}, and \code{cvrt2})
	The input data sets must be provided in the form of \code{\link[MatrixEQTL]{SlicedData}} objects.
}
\usage{
CorrMeta_testAll(
    snps1,
    gene1,
    cvrt1,
    snps2,
    gene2,
    cvrt2,
    pvThreshold)
}
\arguments{
	\item{snps1, snps2}{
		\code{\linkS4class{SlicedData}} objects with genotype information. 
		The slice size should be set to a relatively large number (i.e. 50,000) [see \code{\linkS4class{SlicedData}}].
	}
	\item{gene1, gene2}{
		\code{\linkS4class{SlicedData}} objects with gene expression information. 
		Must have columns matching those of \code{snps1} and \code{snps2} respectively.
		The slice size should be set to a relatively small number (i.e. 50 or 100) if the number of SNPs is large [see \code{\linkS4class{SlicedData}}].
	}
	\item{cvrt1, cvrt2}{
		\code{\linkS4class{SlicedData}} objects with additional covariates. 
		Can be an empty \code{SlicedData} object in case of no covariates.
		The constant is always included in the model and would cause an error if included in \code{cvrt}.
		The order of columns must match those in \code{snps1} and \code{snps2} respectively.
	}
	\item{pvThreshold}{
		\code{numeric}. Significance threshold. Only associations significant at this level will be returned.
	}
}
\details{
	Note that the columns of \code{gene1}, \code{snps1}, and \code{cvrt1} must match. Same for the second data set.
	If they do not match in the input files, use \code{ColumnSubsample} method to subset and/or reorder them.
}
\value{
	Returns data frame with the significant associations.
}
\references{
	The package website: \url{http://www.bios.unc.edu/research/genomic_software/CorrMeta/}
}
\author{
	Andrey Shabalin \email{ashabalin@vcu.edu}
}
\seealso{
	The code below is the sample code for CorrMeta analysis NOT using gene/SNP locations.	
}
\examples{
library(CorrMeta)

pvThreshold = 1e-4;

# Location of data files

data.path = paste0(find.package('CorrMeta'),'/data/');

# Joint data set file names

snps.file = paste0(data.path,'snps.txt');
gene.file = paste0(data.path,'gene.txt');
cvrt.file = paste0(data.path,'cvrt.txt');

# Load joint data sets

snps = SlicedData$new()$LoadFile(snps.file, sliceSize = 50000);
gene = SlicedData$new()$LoadFile(gene.file, sliceSize = 100);
cvrt = SlicedData$new()$LoadFile(cvrt.file);

# Indicator for samples in set 1

set1 = ((1:ncol(gene)) \%\% 2) == 1;

# Create set 1 of data (odd samples), preserving original sets

snps1 = snps$Clone()$ColumnSubsample(set1);
gene1 = gene$Clone()$ColumnSubsample(set1);
cvrt1 = cvrt$Clone()$ColumnSubsample(set1);

# Create set 2 of data (even samples), destroying original sets

snps2 = snps$ColumnSubsample(!set1);
gene2 = gene$ColumnSubsample(!set1);
cvrt2 = cvrt$ColumnSubsample(!set1);
rm(cvrt, gene, snps);

# Run the analysis

result = CorrMeta_testAll(
    snps1, gene1, cvrt1,
    snps2, gene2, cvrt2,
    pvThreshold)

# Display the results

head(result)
}
