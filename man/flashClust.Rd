\name{flashClust}
\alias{flashClust}
\alias{hclust}
\title{ Faster alternative to hclust }
\description{
  This function implements optimal hierarchical clustering with the same interface as
\code{\link{hclust}}. 
}
\usage{
hclust(d, method = "complete", members=NULL)
flashClust(d, method = "complete", members=NULL)
}
\arguments{
  \item{d}{ a dissimilarity structure as produced by 'dist'.}
  \item{method}{ the agglomeration method to be used. This should be (an
          unambiguous abbreviation of) one of \code{"ward"}, \code{"single"},
          \code{"complete"}, \code{"average"}, \code{"mcquitty"}, \code{"median"} or
          \code{"centroid"}. }
   \item{members}{\code{NULL} or a vector with length size of
    \code{d}. See the \sQuote{Details} section.}

}
\details{
  See the description of \code{\link{hclust}} for details on available clustering methods.

  If \code{members!=NULL}, then \code{d} is taken to be a
  dissimilarity matrix between clusters instead of dissimilarities
  between singletons and \code{members} gives the number of observations
  per cluster.  This way the hierarchical cluster algorithm can be
  \sQuote{started in the middle of the dendrogram}, e.g., in order to
  reconstruct the part of the tree above a cut (see examples).
  Dissimilarities between clusters can be efficiently computed (i.e.,
  without \code{hclust} itself) only for a limited number of
  distance/linkage combinations, the simplest one being squared
  Euclidean distance and centroid linkage.  In this case the
  dissimilarities between the clusters are the squared Euclidean
  distances between cluster means.

  \code{flashClust} is a wrapper for compatibility with older code.
}
\value{
  Returned value is the same as that of \code{\link{hclust}}: 
  An object of class \bold{hclust} which describes the
  tree produced by the clustering process.
  The object is a list with components:
  \item{merge}{an \eqn{n-1} by 2 matrix.
    Row \eqn{i} of \code{merge} describes the merging of clusters
    at step \eqn{i} of the clustering.
    If an element \eqn{j} in the row is negative,
    then observation \eqn{-j} was merged at this stage.
    If \eqn{j} is positive then the merge
    was with the cluster formed at the (earlier) stage \eqn{j}
    of the algorithm.
    Thus negative entries in \code{merge} indicate agglomerations
    of singletons, and positive entries indicate agglomerations
    of non-singletons.}

  \item{height}{a set of \eqn{n-1} non-decreasing real values.
    The clustering \emph{height}: that is, the value of
    the criterion associated with the clustering
    \code{method} for the particular agglomeration.}

  \item{order}{a vector giving the permutation of the original
    observations suitable for plotting, in the sense that a cluster
    plot using this ordering and matrix \code{merge} will not have
    crossings of the branches.}

  \item{labels}{labels for each of the objects being clustered.}

  \item{call}{the call which produced the result.}

  \item{method}{the cluster method that has been used.}

  \item{dist.method}{the distance that has been used to create \code{d}
    (only returned if the distance object has a \code{"method"}
    attribute).}

}
\references{ F.Murtagh's software web site, http://astro.u-strasbg.fr/~fmurtagh/mda-sw/splus/ 

  Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988)
  \emph{The New S Language}.
  Wadsworth \& Brooks/Cole. (S version.)

  Everitt, B. (1974).
  \emph{Cluster Analysis}.
  London: Heinemann Educ. Books.

  Hartigan, J. A. (1975).
  \emph{Clustering  Algorithms}.
  New York: Wiley.

  Sneath, P. H. A. and R. R. Sokal (1973).
  \emph{Numerical Taxonomy}.
  San Francisco: Freeman.

  Anderberg, M. R. (1973).
  \emph{Cluster Analysis for Applications}.
  Academic Press: New York.

  Gordon, A. D. (1999).
  \emph{Classification}. Second Edition.
  London: Chapman and Hall / CRC

  Murtagh, F. (1985).
  \dQuote{Multidimensional Clustering Algorithms}, in
  \emph{COMPSTAT Lectures 4}.
  Wuerzburg: Physica-Verlag
  (for algorithmic details of algorithms used).

  McQuitty, L.L. (1966).
  Similarity Analysis by Reciprocal Pairs for Discrete and Continuous
  Data.
  \emph{Educational and Psychological Measurement}, \bold{26}, 825--831.
 }
\author{ Fionn Murtagh, adapted and packaged by Peter Langfelder}
\seealso{ \code{\link{hclust}} }
\examples{

  # generate some data to cluster
  set.seed(1);
  nNodes = 2000;

  # Random "distance" matrix
  dst = matrix(runif(n = nNodes^2, min = 0, max = 1), nNodes, nNodes);

  # Time the flashClust clustering
  system.time( {
     h1 = hclust(as.dist(dst), method= "average");
    } );

  # Time the standard R clustering
  system.time( {
     h2 = stats::hclust(as.dist(dst), method = "average");
    } );

  all.equal(h1, h2)
  # What is different:

  h1[[6]]
  h2[[6]]

  # Everything but the 'call' component is the same; in particular, the trees are exactly equal.

}

\keyword{multivariate}
\keyword{cluster}

