\name{AGB}
\alias{AGB}
\title{Estimate the Above Ground Biomass}
\description{
Estimate the above ground biomass, carbon and CO2 equivalent of trees.
}
\usage{
AGB(x, measure.label, h, taxon="taxon", dead="dead", circumference=TRUE, su="quadrat", area, coord, rm.dead=TRUE, check.spelling=FALSE, correct.taxon=TRUE, sort=TRUE, decreasing=TRUE, long=FALSE)
}
\arguments{
  \item{x}{A data frame containing the community sample data.}
  \item{measure.label}{Name of the column representing the circumference/diameter in centimeters (cm) measured at breast height.}
  \item{h}{Name of the column representing the tree height in meters (m). If omitted in x, the tree height is estimated by a model based on the sample coordinates. See ‘coord’.}
  \item{taxon}{Name of the column representing the sampled taxons. Default is ‘taxon’. Use "UTF-8" character format, accents and special characters are not allowed.}
  \item{dead}{String used to identify the dead individuals. Default is ‘dead’.}
  \item{circumference}{Logical. If TRUE (the default) the function assumes that the circumference at breast height was measured.}
  \item{su}{Name of the column representing the sample unit’s identification. Default is ‘quadrat’.}
  \item{area}{Numeric value (vector of length 1) indicating the total sampled area in hectare (ha).}
  \item{coord}{A vector (ex. c(longitude, latitude)) or a matrix/dataframe with two columns (ex. cbind(longitude, latitude)) giving coordinates, in decimal degrees, of the sites. The coordinates are used to estimate the height-diameter relationship, according to Chave et al. (2014). ‘coord’ must to be given if ‘h’ is missing in ‘x’.}
  \item{rm.dead}{Logical. If TRUE (the default) the biomass of dead individuals is not calculated.}
  \item{check.spelling}{Logical. If TRUE (the default) the spelling of taxon names are compared, and close names are pointed.}
  \item{correct.taxon}{Logical. If TRUE (the default) the taxon names are corrected name using the Taxonomic Name Resolution Service (Boyle et al., 2013).}
  \item{sort}{Logical. If TRUE (the default) the result table is sorted by AGB values.}
  \item{decreasing}{Logical. If TRUE (the default) the result table is sorted in decrescent order.}
  \item{long}{Logical. If FALSE (the default) the ‘tree’ table is not returned. See ‘Value’.}
}
\details{
AGB function is a wrap of the functions getWoodDensity, computeAGB, and correctTaxo from BIOMASS package (Réjou‐Méchain et al., 2017). Tree biomasses are calculated according to the allometric model by Chave et al. (2014).
It is expected that taxon names are composed by two words (genus species). The function split the taxon name in two columns, i.e., ‘genus’ and ‘species’, which are used to get the wood density data. If taxon name is a single word (ex. Indet, Dead etc), species column is <NA>.
Wood density (in g/cm3) is extracted from a reference database of ~ 16,500 tree species. If a species name is not in the database, the mean value for the genus is used. If the genus is not in the database, the mean value of the sample unit (su) is used.
The data frame (x) must contain individual columns (variables) representing the sample unit labels, taxon names, cbh or dbh measures and tree height for each sampled tree. If tree height is not given, ‘coord’ argument is mandatory.
The column for cbh/dbh accepts the traditional notation for multiple trunks, e.g. ‘17.1+8+5.7+6.8’. The plus sign indicates the delimitator for each trunk measure. Decimal delimitator can be point or comma and, spaces after or before “+” are ignored by the function.
Column names of the data frame (x) are coerced to low cases when the function runs, so it is case insensitive. If ‘x’ contains the default column names, the arguments ‘taxon’, ‘dead’ and ‘su’ can be omitted.
Measure units: Individual cbh/dbh must be in centimeters. Tree height must be in meters.
}
\value{
The function returns a ‘biomass’ object with four data frames. ‘$tree’ gives the individual data of taxon name, wood density (wd), level at which wood density has been calculated (i.e., species, genus, or average of sample unit), biomass above ground (AGB), carbon content (C) and CO2 equivalent (CO2e). ‘$taxon’ includes the above ground biomass (AGB), carbon content (C) and CO2 equivalent (CO2e) per taxon in Mg (or ton) per ha. ‘$total’ gives the total AGB, C and CO2e per ha (irrespective of the taxa). 'DW.level' gives the percentage of wood density calculated at species, genus and sample unit levels.
}
\references{
Boyle, B. et al. (2013). The taxonomic name resolution service: An online tool for automated standardization of plant names. BMC bioinformatics, 14, 1. doi:10.1186/1471-2105-14-16
Chave et al. 2014. Improved allometric models to estimate the aboveground biomass of tropical trees, Global Change Biology, 20 (10), 3177-3190
Réjou‐Méchain, M., Tanguy, A., Piponiot, C., Chave, J., Hérault, B., 2017. biomass: an  package for estimating above‐ground biomass and its uncertainty in tropical forests. Methods Ecol Evol 8, 1163–1167. https://doi.org/10.1111/2041-210X.12753
Zanne, A. E. et al. 2009. Global wood density database. Dryad. Identifier: http://hdl. handle. net/10255/dryad 235.
}
\examples{
\dontrun{
data <- quadrat.df

resul1<-AGB(data, measure.label="cbh", h="h", taxon="taxon", dead="Dead", circumference=TRUE, su="quadrat", area=0.0625, rm.dead=TRUE, check.spelling=F, correct.taxon=T, sort=T, decreasing=T, long=T)

head(resul1$tree)
resul1$taxon
resul1$total
resul1$WD.level

# Using default arguments
Resul2<-AGB(data, measure.label="cbh", circumference=TRUE, h="pap", area=0.0625)

# Missing tree height
Resul3<-AGB(data, measure.label="pap", taxon="especie", dead="Morta", circumference=TRUE, su="parcela", area=0.0625, coord=c(-47.85, -21.17))
}
}
\author{Rodrigo A. S. Pereira}
\keyword{agb}
