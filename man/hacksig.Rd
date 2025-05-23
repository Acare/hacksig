% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hacksig-package.R
\docType{package}
\name{hacksig}
\alias{hacksig}
\alias{hacksig-package}
\title{hacksig: A Tidy Framework to Hack Gene Expression Signatures}
\description{
The \code{hacksig} package has been designed for the purpose of simplifying the way
in which gene expression signatures scores are computed. It is a manually
curated collection of gene expression signatures found in literature and makes
use of three different single sample score calculation methods. Moreover,
parallel computation is supported through the \code{future} framework.
}
\section{Get gene signatures scores in different ways}{

The main function of the package is \code{\link[=hack_sig]{hack_sig()}} and it can be used to:
\itemize{
\item obtain single sample scores with one of three methods (z-score, ssGSEA,
singscore) for a custom list of gene signatures;
\item obtain single sample scores for a number of manually curated gene signatures
either with the original publication method or with one of the three single
sample methods.
}

Once single sample scores are obtained, you can assign your samples into
signature classes with \code{\link[=stratify_sig]{stratify_sig()}}.

In addition, other more complex methods are implemented through:
\itemize{
\item \code{\link[=hack_cinsarc]{hack_cinsarc()}}, for the CINSARC classification;
\item \code{\link[=hack_estimate]{hack_estimate()}}, for the ESTIMATE method;
\item \code{\link[=hack_immunophenoscore]{hack_immunophenoscore()}}, for the Immunophenoscore.
}

Information about implemented signatures can be obtained with \code{\link[=get_sig_info]{get_sig_info()}}.
}

\section{Check if gene signatures are applicable to your data}{

Sometimes your gene expression matrix can miss some genes due to some prior
filtering procedure. The function \code{\link[=check_sig]{check_sig()}} can be used to check how many
genes your expression matrix contain for every input signature.
}

\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/Acare/hacksig}
  \item \url{https://acare.github.io/hacksig/}
  \item Report bugs at \url{https://github.com/Acare/hacksig/issues}
}

}
\author{
Andrea Carenzo, Loris De Cecco, Federico Pistore
}
\keyword{internal}
