#' Flow cytometry data
#'
#' The Raf network is implied in the regulation of cellular proliferation.
#' The data were collected by (Sachs et al., 2005).
#' Flow cytometry measurements consist in sending unique cells suspended in a fluid through a laser beam, and measuring
#' parameters of interest by collecting the light re-emitted by the cell by diffusion or fluores- cence.
#' In this study, the parameters of interest are the activation level of 11 proteins and phospholipids involved in the Raf pathway,
#'  and are measured by flow cytometry across 100 different cells.
#'
#' @docType data
#'
#' @usage data(cyto)
#'
#' @format A dataframe \code{X.raf} with 100 lines (cells) and 11 columns (proteins)
#'  and an adjacency matrix \code{A.raf} describing the network.
#'
#' @keywords datasets
#'
#' @references Sachs, K., O. Perez, D. Pe’er, D. A. Lauffenburger, and G. P. Nolan (2005).
#' Causal protein- signaling networks derived from multiparameter single-cell data. Science 308, 523–529.
#'
#' @examples
#' data(cyto)
#' \donttest{plot.igraph(graph_from_adjacency_matrix(A.raf),edge.arrow.size=0, vertex.color="gold", vertex.size=15, vertex.frame.color="gray", vertex.label.color="black", vertex.label.cex=0.8,vertex.label.dist=0)}
#' @name cyto
#' @author Christophe Ambroise \email{christophe.ambroise@@univ-evry.fr}
#' @keywords data
NULL
