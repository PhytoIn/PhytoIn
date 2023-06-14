\name{BAplot}
\alias{BAplot}
\title{Plot basal areas on a map of quadrats}
\description{
Plot basal areas of trees on a map of quadrats. If individual tree coordinates are not known, the coordinates inside the quadrats are randomly defined, according to the uniform distribution.
}
\usage{
BAplot(formula, data, taxon="taxon", circumference=TRUE, quadrat.size, dead="dead", rm.dead=FALSE, origin=c(0, 0), col="grey40", alpha=1, cex.radius=1, ind.coord=FALSE, legend=TRUE, long=FALSE)
}
\arguments{
\item{formula}{
A model formula indicating the trunk measure (circumference [default] or diameter) in centimeters and the xy coordinates of each quadrat the tree belongs (‘ind.coord=FALSE’) or, the actual tree coordinates (‘ind.coord=TRUE’). Example: ‘measure ~ x + y’. See Details.
}
\item{data}{
A data frame containing the community sample data. See Details.
}
\item{taxon}{
Name of the column representing the sampled taxa. Default is ‘taxon’.
}
\item{circumference}{
Logical. If TRUE (the default) the function assumes that the circumference at breast height was measured.
}
\item{quadrat.size}{
A vector indicating the side lengths (in meters) of the x and y quadrat sides (ex. c(x, y)). It can be indicated as unique value if the quadrat is a square.
}
\item{dead}{
String used to identify the dead individuals. Default is ‘dead’.
}
\item{rm.dead}{
Logical. If FALSE (the default) basal areas of dead individuals are plotted.
}
\item{origin}{
A numeric vector indicating the map origin coordinates. Default is ‘c(0, 0)’.
}
\item{col}{
Circle color of represented basal areas. This argument has value only if legend=FALSE.
}
\item{alpha}{
Value of transparency factor: zero (100% transparent) – 1 (no transparency). Default is 1.
}
\item{cex.radius}{
A numerical value giving the amount by which the tree radius should be magnified relative to the actual measure. Default is 1.
}
\item{ind.coord}{
Logical indicating whether the individual coordinates are given. If FALSE (the default) the tree coordinates inside the quadrats are randomly defined. If TRUE, the actual tree coordinates are plotted.
}
\item{legend}{
Logical. If TRUE (default) the taxon is used as color legend, otherwise circle color will be defined by the argument ‘col’ (default is “blue”).
}
\item{long}{
Logical. If FALSE (default) the function does not return the result data frame, which contains xy coordinates, radius and taxon name of each sampled tree.
}
}
\details{
‘BAplot’ uses the function ‘circleRepelLayout()’ from ‘packcircles’ package to rearrange circle coordinates to avoid overlapping. The minimum distance allowed among trees is 1 meter.
The packages ‘ggforce’ and ‘ggplot2’ are used to draw the map.
The data frame passed to ‘data’ argument must include two columns indicating x and y coordinates of each quadrat that the tree belongs or the actual tree coordinates. If actual coordinates are supplied, the ‘ind.coord’ argument must be set TRUE.
Circumference/diameter measures accept the traditional notation for multiple trunks, e.g. ‘17.1+8+5.7+6.8’. The plus sign indicates the delimitator for each trunk measure. Decimal delimitator can be point or comma and, spaces after or before “+” are ignored by the function.
}
\value{
A plot representing the quadrat map and tree basal areas. If ‘long=TRUE’, the function returns a dataframe containing taxon name, xy coordinates and radius of each sampled tree.
}
\references{
\code{packcircles} package \url{https://github.com/mbedward/packcircles}
\code{ggforce} package \url{https://github.com/thomasp85/ggforce}
}
\examples{
# Using random coordinates
data <- quadrat1.df
BAplot (formula= pap ~ x + y, data=data, taxon="especie", circumference=TRUE, quadrat.size=5, dead="Morta", rm.dead=FALSE, col="blue", alpha=0.4, cex.radius=2, legend=FALSE, long=FALSE, ind.coord=FALSE)

# Using actual coordinates
data2 <- quadrat2.df
BAplot (formula= pap ~ x + y, data=data, taxon="especie", circumference=TRUE, quadrat.size=5, dead="Morta", rm.dead=FALSE, col="blue", alpha=0.4, cex.radius=2, legend=FALSE, long=FALSE, ind.coord=TRUE)
}
\author{Rodrigo A. S. Pereira}
\keyword{BAplot}
