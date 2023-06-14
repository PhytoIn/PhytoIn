\name{phytoparam}
\alias{phytoparam}
\title{Estimate the Phytosociological Parameters and Diversity Indices}
\description{
Estimate the phytosociological parameters and the Shannon-Wiener, Pielou and Simpson’s diversity indices, using the quadrat or the point-centered quarter methods.
}
\usage{
phytoparam(x, measure.label, h ="h", taxon="taxon", family="family", dead="dead", circumference=TRUE, su="quadrat", height=FALSE, quadrat=TRUE, su.size, d="distance", shape.factor=1, rm.dead=FALSE, check.spelling=TRUE)
}
\arguments{
  \item{x}{A data frame containing the community sample data.}
  \item{measure.label}{Name of the column representing the circumference/diameter measured at breast height. If omitted the function assumes the default names ‘cbh’ or ‘dbh’ for circumference or diameter at breast height, respectively. See ‘circumference’ argument.}
  \item{h}{Name of the column representing the trunk height. Default is ‘h’.}
  \item{taxon}{Name of the column representing the sampled taxons. Default is ‘taxon’. Use "UTF-8" character format, accents and special characters are not allowed.}
  \item{family}{Name of the column representing the family names of the sampled taxons. Default is ‘family’. This argument is used to calculate the number of individuals and number of species per family. If you do not want these parameters, set family = NA. Use "UTF-8" character format, accents and special characters are not allowed.}
  \item{dead}{String used to identify the dead individuals. Default is ‘dead’.}
  \item{circumference}{Logical. If TRUE (the default) the function assumes that the circumference at breast height was measured.}
  \item{su}{Name of the column representing the sample units identification. Default is ‘quadrat’ for quadrat method and ‘point’ for point-centered quarter method.}
  \item{height}{Logical. If FALSE (the default) trunk volume is not calculated.}
  \item{quadrat}{Logical. If TRUE (the default) the function assumes that data was sampled using the quadrat method. If FALSE, it assumes that data was sampled using the centered-point quarter method.}
  \item{su.size}{Numeric of length 1 indicating the quadrat area in m², necessary only if quadrat= TRUE.}
  \item{d}{Name of the column representing the point to tree distance, necessary only if quadrat=FALSE. Default is ‘distance’.}
  \item{shape.factor}{A number 0 < shape.factor <= 1 indicating the trunk shape. Factor = 1 considers the trunk has a perfect cylindrical form.}
  \item{rm.dead}{Logical. If FALSE (the default) the phytosociological parameters for the dead individuals are calculated.}
  \item{check.spelling}{Logical. If TRUE (the default) taxon names are checked for misspelling.}
}
\details{
The function estimates the phytosociological parameters for plant tree communities sampled by quadrat or point-centered quarter methods (‘quadrat=T’ or ‘quadrat=F’ respectively).
For quadrat method, data frame (x) must contain individual columns (variables) representing the sample unit labels, taxon names, cbh or dbh measures for each sampled tree. Additionally, trunk height and family can be included to estimate volume and family parameters.
For point-centered quarter method, data frame (x) must contain (in addition to the mandatory columns for the quadrat method) a column for the distance from the point to the individual.
The column for cbh/dbh accepts the traditional notation for multiple trunks, e.g. ‘17.1+8+5.7+6.8’. The plus sign indicates the delimitator for each trunk measure. Decimal delimitator can be point or comma and, spaces after or before “+” are ignored by the function.
Column names of the data frame (x) are coerced to low cases when the function runs, so it is case insensitive. If ‘x’ contains the default column names, the arguments ‘h’, ‘taxon’, ‘family’, ‘dead’, ‘su’ and ‘d’ can be omitted.
Unbiased absolute density in the point-centered quarter method is estimated according to Pollard (1971) and Seber (1982).
Measure units: Individual cbh/dbh must be in centimeters. Trunk height and point to individual distance must be in meters.
}
\value{
The function returns a ‘param’ object containing two or four data frames. ‘$global’ includes the total parameters and the diversity indices. Sampled area is in hectare (ha), total density is in individual/ha, total dominance is in m³/ha and Shannon-Wiener (H’) is in nats/individual (natural log).
‘$param’ includes the taxon identification, number of sampled individuals/taxon (N), absolute and relative density (Ade and RDe), absolute and relative frequency (AFr and RFr), absolute and relative dominance (ADo and RDo), absolute and relative volume (AVol and RVol), index of importance value (IV) and index of cover value (CV). Absolute parameters are in individual/ha and relative parameters are in percentage.
If ‘family’ != NA, two more data frames are returned, listing the number individuals (‘ind.fam’) and number of species (‘spp.fam’) per family.

}
\references{
Pollard, J. H. 1971. On distance estimators of density in randomly distributed forests. Biometrics. 27: 991-1002.
Seber, G. A. F. 1982. The Estimation of Animal Abundance and Related Parameters. New York, MacMillan. 41-45.
}
\examples{
\dontrun{
# Quadrat method
names(quadrat.df)
[1] "parcela" "Especie" "PAP"     "h"

phytoparam(x= quadrat.df, measure.label="pap", taxon="especie", dead="morta", circumference=TRUE, su="parcela", height=T, su.size=25)  # using customized column names

quadrat.default<-quadrat.df
colnames(quadrat.default)<-c("quadrat", "taxon", "cbh", "h")
phytoparam(x=quadrat.default, dead="morta", circumference=TRUE, height=T, su.size=25)  # using default column names

# Point-centered quarter method
names(point.df)
[1] "Ponto"     "Especie"   "Distancia" "PAP"       "h"

phytoparam(x= point.df, measure.label="pap", taxon="especie", dead="morta", circumference=T, su="ponto", height=T, quadrat=F, d="distancia")  # using customized column names
point.default<- point.df
colnames(point.default)<-c("point", "taxon", "distance", "cbh", "h")

phytoparam(x= point.default, dead="morta", circumference=T, height=T, quadrat=F)  
}
}
\author{Rodrigo A. S. Pereira}
\keyword{phytosociology}
