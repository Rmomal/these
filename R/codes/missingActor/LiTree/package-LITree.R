#' LITree: A package for infering Gaussian Graphical Models with unobserved variables using latent tree agreggation
#'
#'
#' The LITree package provides four classes  graphModel, GGMmodel, GGMfit,  GGMexperiment and one function em.litree()
#'
#' The function em.litree implements the algorithm of Gaussian Graphical Model Inference with missing variable described
#' in Robin et. al (2017). The underlying model is based on the aggregation of spanning trees, and the estimation procedure
#' on the Expectation-Maximization algorithm. We treat the graph structure and the unobserved nodes as missing variables
#' and compute posterior probabilities of edge appearance. To provide a complete methodology,
#' we also propose three model selection criteria to estimate the number of missing nodes.
#'
#' @author Geneviève Robin \email{genevieve.robin@@polytechnique.edu}
#' @author Christophe Ambroise \email{christophe.ambroise@@univ-evry.fr}
#' @author Stéphane Robin \email{stéphane.robin@@agroparistech.edu}
#' @references Geneviève Robin, Christophe Ambroise, Stéphane Robin (Submitted on 26 May 2017).
#' Graphical model inference with unobserved variable via latent tree aggregation.
#' \url{Arxiv Paper. https://arxiv.org/abs/1705.09464}


#' @docType package
#' @name LITree
NULL
