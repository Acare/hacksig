% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hack_immunophenoscore.R
\name{hack_immunophenoscore}
\alias{hack_immunophenoscore}
\title{Hack the Immunophenoscore}
\source{
\href{https://github.com/icbi-lab/Immunophenogram}{github.com/icbi-lab/Immunophenogram}
}
\usage{
hack_immunophenoscore(expr_data, extract = "ips")
}
\arguments{
\item{expr_data}{A normalized gene expression matrix (or data frame) with
gene symbols as row names and samples as columns.}

\item{extract}{A string controlling which type of biomarker scores you want
to obtain. Possible choices are:
\itemize{
\item \code{"ips"} (default), only raw and discrete IPS scores;
\item \code{"class"}, IPS scores together with the four summary class scores;
\item \code{"all"}, all possible biomarker scores.
}}
}
\value{
A tibble with one row for each sample in \code{expr_data}, a column \code{sample_id}
indicating sample identifiers and a number of additional columns depending
on the choice of \code{extract}.
}
\description{
Obtain various immune biomarkers scores, which combined together give the
immunophenoscore (\emph{Charoentong et al., 2017}).
}
\details{
The immunophenoscore is conceived as a quantification of tumor immunogenicity.
It is obtained by aggregating multiple immune biomarkers scores, which are
grouped into four major classes:
\itemize{
\item \emph{MHC molecules} (\strong{MHC}), expression of MHC class I, class II, and non-classical molecules;
\item \emph{Immunomodulators} (\strong{CP}), expression of certain co-inhibitory and co-stimulatory molecules;
\item \emph{Effector cells} (\strong{EC}), infiltration of activated CD8+/CD4+ T cells and Tem (effector memory) CD8+/CD4+ cells;
\item \emph{Suppressor cells} (\strong{SC}), infiltration of immunosuppressive cells (Tregs and MDSCs).
}

The table below shows in detail the 26 immune biomarkers and cell types grouped
by class together with the number of genes which represent them:\tabular{lcc}{
   \strong{Class} | \tab \strong{Biomarker/cell type} | \tab \strong{No. genes} \cr
   MHC \tab B2M \tab 1 \cr
   MHC \tab HLA-A \tab 1 \cr
   MHC \tab HLA-B \tab 1 \cr
   MHC \tab HLA-C \tab 1 \cr
   MHC \tab HLA-DPA1 \tab 1 \cr
   MHC \tab HLA-DPB1 \tab 1 \cr
   MHC \tab HLA-E \tab 1 \cr
   MHC \tab HLA-F \tab 1 \cr
   MHC \tab TAP1 \tab 1 \cr
   MHC \tab TAP2 \tab 1 \cr
   CP \tab CD27 \tab 1 \cr
   CP \tab CTLA-4 \tab 1 \cr
   CP \tab ICOS \tab 1 \cr
   CP \tab IDO1 \tab 1 \cr
   CP \tab LAG3 \tab 1 \cr
   CP \tab PD1 \tab 1 \cr
   CP \tab PD-L1 \tab 1 \cr
   CP \tab PD-L2 \tab 1 \cr
   CP \tab TIGIT \tab 1 \cr
   CP \tab TIM3 \tab 1 \cr
   EC \tab Act CD4 \tab 24 \cr
   EC \tab Act CD8 \tab 26 \cr
   EC \tab Tem CD4 \tab 27 \cr
   EC \tab Tem CD8 \tab 25 \cr
   SC \tab MDSC \tab 20 \cr
   SC \tab Treg \tab 20 \cr
}
}
\section{Algorithm}{

Samplewise gene expression z-scores are obtained for each of 26 immune cell
types and biomarkers. Then, weighted averaged z-scores are computed for each
class and the raw immunophenoscore (\eqn{IPS_{raw}}) results as the sum of the
four class scores. Finally, the immunophenoscore (\eqn{IPS}) is given as an
integer value between 0 and 10 in the following way:
\itemize{
\item \eqn{IPS = 0}, if \eqn{IPS_{raw} \le 0};
\item \eqn{IPS = [10 * (IPS_{raw} / 3)]}, if \eqn{0 < IPS_{raw} < 3};
\item \eqn{IPS = 10}, if \eqn{IPS_{raw} \ge 3}.
}
}

\examples{
hack_immunophenoscore(test_expr)
hack_immunophenoscore(test_expr, extract = "class")
}
\references{
Charoentong, P., Finotello, F., Angelova, M., Mayer, C.,
Efremova, M., Rieder, D., Hackl, H., & Trajanoski, Z. (2017). Pan-cancer
Immunogenomic Analyses Reveal Genotype-Immunophenotype Relationships and
Predictors of Response to Checkpoint Blockade. \emph{Cell reports}, 18(1), 248–262.
\doi{10.1016/j.celrep.2016.12.019}.
}
\seealso{
\code{\link[=hack_sig]{hack_sig()}} to compute Immunophenoscore biomarkers in different
ways (e.g. use \code{signatures = "ips"} and \code{method = "singscore"}).

\code{\link[=check_sig]{check_sig()}} to check if all/most of the Immunophenoscore biomarkers are
present in your expression matrix (use \code{signatures = "ips"}).
}
