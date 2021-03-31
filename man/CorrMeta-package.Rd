\name{CorrMeta-package}
\alias{CorrMeta}
\docType{package}
\title{
   CorrMeta: Fast Association Analysis for eQTL and GWAS
   Data with Related Samples and Correlated Phenotypes
}
\description{
Conducts eQTL analysis on correlated data by first splitting
samples into two independent groups on which multiple linear regression analysis
is performed separately, followed by an appropriate meta-analysis
to combine the two non-independent test results.
}
\details{
	\tabular{ll}{
		Package: \tab CorrMeta\cr
		Type: \tab Package\cr
		Version: \tab 1.0\cr
		Date: \tab 2015-09-01\cr
		License: \tab LGPL-3 \cr
		LazyLoad: \tab yes\cr
		Depends: \tab methods, utils, stats\cr
	}
}
\author{
	Andrey Shabalin \email{ashabalin@vcu.edu}
	
	Maintainer: Andrey Shabalin <ashabalin@vcu.edu>
}
\references{
	The package website: \url{http://www.bios.unc.edu/research/genomic_software/CorrMeta/}
}
\keyword{ package }
\keyword{ CorrMeta }
\seealso{
	See \code{\link{CorrMeta_testAll}} for reference and sample code.
}
