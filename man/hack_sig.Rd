% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hack_sig.R
\name{hack_sig}
\alias{hack_sig}
\title{Score samples by gene signatures}
\usage{
hack_sig(
  expr_data,
  signatures = "all",
  method = "original",
  direction = "none",
  sample_norm = "raw",
  rank_norm = "none",
  alpha = 0.25
)
}
\arguments{
\item{expr_data}{A normalized gene expression matrix (or data frame) with
gene symbols as row names and samples as columns.}

\item{signatures}{It can be a list of signatures or a character vector indicating
keywords for a group of signatures. The default (\code{"all"}) will cause the
function to compute single sample scores for all the signatures implemented
in \code{hacksig}.}

\item{method}{A character string specifying which method to use for computing
the single sample score for each signature. You can choose one of:
\itemize{
\item \code{"original"}, the original method used by the authors of the signature;
\item \code{"zscore"}, the combined z-score method;
\item \code{"ssgsea"}, the single sample GSEA method;
\item \code{"singscore"}, the singscore method;
}}

\item{direction}{A character string specifying the \strong{singscore} computation
method depending on the direction of the signatures. Can be on of:
\itemize{
\item \code{"none"} (default), undirected signatures, that is you do not know whether
the genes are up- or down-regulated;
\item \code{"up"}, all genes in the signature are supposed to be up-regulated;
\item \code{"down"}, all genes in the signature are supposed to be down-regulated;
}}

\item{sample_norm}{A character string specifying the type of normalization
affecting the \strong{single sample GSEA} scores. Can be one of:
\itemize{
\item \code{"raw"} (default), obtain raw scores;
\item \code{"separate"}, normalize raw scores in \eqn{[0, 1]} across samples for
each signature separately.
\item \code{"all"}, normalize raw scores both across samples and signatures.
}}

\item{rank_norm}{A character string specifying how gene expression ranks should
be normalized in the \strong{single sample GSEA} procedure. Valid choices are:
\itemize{
\item \code{"none"} (default), no rank normalization;
\item \code{"rank"}, ranks are multiplied by \code{10000 / nrow(expr_data)};
\item \code{"logrank"}, normalized ranks are logged.
}}

\item{alpha}{A numeric scalar. Exponent in the running sum of the \strong{single sample GSEA}
score calculation which weighs the gene ranks. Defaults to \eqn{\alpha = 0.25}.}
}
\value{
A tibble with one row for each sample in \code{expr_data}, a column \code{sample_id}
indicating sample identifiers and one column for each input signature giving
single sample scores.
}
\description{
Compute gene signature single sample scores in one of different ways.
You can choose to apply either the \emph{original} procedure or one of three single
sample scoring methods: the combined z-score (\emph{Lee et al., 2008}), the single
sample GSEA (\emph{Barbie et al., 2009}) or the singscore method (\emph{Foroutan et al., 2018}).
}
\details{
For \code{"original"} method, it is intended the procedure used in the original
publication by the authors for computing the signature score.
\code{hack_sig()} can compute signature scores with the original method only if
this is a relatively simple procedure (e.g weighted sum of fitted model
coefficients and expression values).
For more complex methods, such as CINSARC, ESTIMATE and Immunophenoscore,
use the dedicated functions.

If \code{signatures} is a custom list of gene signatures, then the \code{"ssgsea"}
method will be applied by default.
}
\section{Algorithm}{

This section gives a brief explanation of how single sample scores are obtained
from different methods.
\subsection{Combined z-score}{

Gene expression values are centered by their mean value and scaled by their
standard deviation across samples for each gene (z-scores). Then, for each
sample and signature, corresponding z-scores are added up and divided by the
square root of the signature size (i.e. the number of genes composing a signature).

The combined z-score method is also implemented in the R package \code{GSVA}
(\emph{Hänzelmann et al., 2013}).
}

\subsection{Single sample GSEA}{

For each sample, genes are ranked by expression value in increasing order and
rank normalization may follow (see argument \code{rank_norm}). Then, two probability-like
vectors are computed for each sample and signature:
\itemize{
\item \eqn{P_{in}}, the cumulative sum of weighted ranks divided by their total
sum for genes in the signature;
\item \eqn{P_{out}}, the cumulative sum of ones (indicating genes not in the signature)
divided by the number of genes not in the signature.
}

The single sample GSEA score is obtained by adding up the elements of the
vector difference \eqn{P_{in} - P_{out}}.
Finally, single sample scores could be normalized either across samples or across
gene signatures and samples.

The single sample GSEA method is also implemented in the R package \code{GSVA}
(\emph{Hänzelmann et al., 2013}).
}

\subsection{Singscore}{

For signatures whose genes are supposed to be up- or down-regulated, genes
are ranked by expression value in increasing or decreasing order, respectively.
For signatures whose direction is unknown, genes are ranked by absolute expression
in increasing order and are median-centered.
Enrichment scores are then computed for each sample and signature by averaging
gene ranks for genes in the signature.
Finally, normalized scores are obtained by subtracting the theoretical minimum
mean rank from the score and dividing by the difference between the theoretical
maximum and minimum mean ranks.

The \code{hacksig} implementation of this method works only with unidirectional (i.e.
all genes up- or down-regulated) and undirected gene signatures.
If you want to get single sample scores for bidirectional gene signatures (i.e.
signatures composed of both up- and down-regulated genes), please use the R
package \code{singscore} (\emph{Foroutan et al., 2018}).
}
}

\examples{
# Raw ssGSEA scores for all implemented signatures can be obtained with:
\donttest{hack_sig(test_expr, method = "ssgsea")}

# To obtain 0-1 normalized ssGSEA scores, use:
\donttest{hack_sig(test_expr, method = "ssgsea", sample_norm = "separate")}

# You can also change the exponent of the ssGSEA running sum with:
\donttest{hack_sig(test_expr, method = "ssgsea", sample_norm = "separate", alpha = 0.5)}

# To obtain combined z-scores for custom gene signatures, use:
custom_list <- list(rand_sig1 = rownames(test_expr)[1:5],
                    rand_sig2 = c(rownames(test_expr)[6:8], "RANDOMGENE"))
hack_sig(test_expr, custom_list, method = "zscore")
}
\references{
Barbie, D. A., Tamayo, P., Boehm, J. S., Kim, S. Y., Moody, S. E., Dunn, I. F.,
Schinzel, A. C., Sandy, P., Meylan, E., Scholl, C., Fröhling, S., Chan, E. M.,
Sos, M. L., Michel, K., Mermel, C., Silver, S. J., Weir, B. A., Reiling, J. H.,
Sheng, Q., Gupta, P. B., … Hahn, W. C. (2009). Systematic RNA interference
reveals that oncogenic KRAS-driven cancers require TBK1. \emph{Nature}, 462(7269),
108–112. \doi{10.1038/nature08460}.

Foroutan, M., Bhuva, D. D., Lyu, R., Horan, K., Cursons, J., &
Davis, M. J. (2018). Single sample scoring of molecular phenotypes.
\emph{BMC bioinformatics}, 19(1), 404. \doi{10.1186/s12859-018-2435-4}.

Hänzelmann, S., Castelo, R., & Guinney, J. (2013). GSVA: gene set variation
analysis for microarray and RNA-seq data. \emph{BMC bioinformatics}, 14, 7.
\doi{10.1186/1471-2105-14-7}.

Lee, E., Chuang, H. Y., Kim, J. W., Ideker, T., & Lee, D. (2008). Inferring
pathway activity toward precise disease classification.
\emph{PLoS computational biology}, 4(11), e1000217.
\doi{10.1371/journal.pcbi.1000217}.
}
\seealso{
\code{\link[=get_sig_info]{get_sig_info()}} to get information about all implemented signatures.

\code{\link[=check_sig]{check_sig()}} to check if signatures are applicable to your data.

\code{\link[=hack_cinsarc]{hack_cinsarc()}} to apply the original CINSARC procedure.

\code{\link[=hack_estimate]{hack_estimate()}} to obtain the original ESTIMATE scores.

\code{\link[=hack_immunophenoscore]{hack_immunophenoscore()}} to apply the original Immunophenoscore procedure.
}
