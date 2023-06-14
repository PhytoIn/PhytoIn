AGB <- function(x, measure.label, h, taxon="taxon", dead="dead", circumference=TRUE, su="quadrat", area, coord, 
                rm.dead=TRUE, check.spelling=FALSE, correct.taxon=TRUE, sort=TRUE, decreasing=TRUE, long=FALSE)
  
{
  if(all((.packages())!= "BIOMASS")) require(BIOMASS)
  measure.label <- tolower(measure.label)
  taxon <- tolower(taxon)
  if(!missing(h)) h<- tolower(h)
  dead<- tolower(dead)
  su<- tolower(su)
  
  colnames(x) <- tolower(colnames(x)) # renomeia as colunas do data.frame usando apenas letras minusculas
  y<-is.na(x) | x == "" # Teste logico para buscar todas as celulas vazias ou com NAs
  x <- x[!apply(y == TRUE, 1, all), !apply(y == TRUE, 2, all)]  # Remove todas as linhas/colunas vazias/NAs
  
  # Check taxon spelling
  if(check.spelling){
    initial.list<-unique(x[[taxon]])
    spp.list<-initial.list
    for(k in 1:length(spp.list)){
      close.strings<-agrep(initial.list[k], spp.list, value=TRUE)
      if(length(close.strings)>1)
      {
        opt<-menu(close.strings, title="Choose the number corresponding the correct name. Type 0 to skip.")
        correct.name<-close.strings[opt]
        incorrect<-close.strings[-opt]
        if(opt!=0) for(i in 1:length(incorrect)) {
          spp<-gsub(incorrect[i], correct.name, x[[taxon]])
          spp.list<-unique(gsub(incorrect[i], correct.name, spp.list))
        }
      }
    }
  }
  # Remove dead individuals
  if (rm.dead) # se a opcao for remover as plantas mortas do data frame 
  {
    filter = tolower(x[[taxon]]) == tolower(dead) # testa quais individuos sao "mortos"
    x <- x[!filter, ] # remocao das linhas contendo plantas mortas do data frame
    cat ("\n", sum(filter), "dead individuals removed from the dataset \n") # exibe na tela uma mensagem com o numero de "mortas" removidas do data frame
  }
  
  # Split taxon in genus and species
  x[[taxon]] <- gsub(" +", " ", x[[taxon]])
#  list<-strsplit(x[[taxon]], split=" ", fixed=T)
#  taxon.names<-as.data.frame(matrix(unlist(list), ncol=2, byrow = T))
#  if(dim(taxon.names)[1]!=dim(x)[1]) stop("Check multiple spaces in the taxon column")
#  colnames(taxon.names)<-c("genus","species")
  
  list <- lapply(x[[taxon]], function(z) {
    split_result <- strsplit(z, split = " ", fixed = TRUE)[[1]]
    if (length(split_result) > 1) {
      split_result
    } else {
      c(split_result, NA)
    }
  })
  taxon.names <- data.frame(genus = sapply(list, "[", 1), species = sapply(list, "[", 2))
  
  # Correct taxon
  if(correct.taxon) {
    taxon.corrected<-correctTaxo(genus=taxon.names$genus, species=taxon.names$species, useCache=T)
    taxon.names <- data.frame(genus=taxon.corrected$genusCorrected, species=taxon.corrected$speciesCorrected)
      }
  #####################
  # Individual diameter
  
  measure <- x[[measure.label]]  # extrai a coluna de valores das medidas
  measure <- gsub(",", ".", measure) # substitui eventuais vigulas por ponto
  
  # Separa os fustes multiplos em colunas
  my_list<-strsplit(x=measure, split="+", fixed=TRUE)
  n.obs <- sapply(my_list, length)
  seq.max <- seq_len(max(n.obs))
  mat <- t(sapply(my_list, "[", i = seq.max))
  measure.table<-matrix(as.numeric(mat), ncol = max(seq.max))    # Convert to numeric matrix
  measure.table[is.na(measure.table)] <- 0
  if (circumference) AB <- (measure.table)^2/(4*pi) else 
    AB <- (pi*(measure.table)^2)/4  # calcula a �rea basal em m� de cada fuste a partir dos valores de cap ou dap
  ABi <- apply(AB, 1, sum)        # Calcula a area basal individual
  d <- sqrt(ABi/pi)*2   # Calcula o diametro individual
  
  # Wood density
  WDdata<-getWoodDensity(genus=taxon.names$genus, species=taxon.names$species, stand=x[[su]])
  
  # Biomass
  if(missing(h)) AGBtree<-computeAGB(D=d, WD=WDdata$meanWD, coord=coord) else 
    AGBtree<-computeAGB(D=d, WD=WDdata$meanWD, H=x[[h]])
  res.tree <- data.frame(su=x[[su]], taxon=x[[taxon]], wd=WDdata$meanWD, level=WDdata$levelWD, 
                         AGB=AGBtree, C=AGBtree*0.471, CO2e=AGBtree*0.471*3.6705069)
  AGBsp <- aggregate(cbind(AGB=AGB/area, C=C/area, CO2e=CO2e/area) ~ taxon, data = res.tree, sum)
  if(sort) AGBsp <- AGBsp[order(AGBsp$AGB, decreasing=decreasing), ]
  row.names(AGBsp) <- 1:dim(AGBsp)[1]
  AGBtot <- cbind(AGBtot=sum(AGBsp$AGB), Ctot=sum(AGBsp$C), CO2etot=sum(AGBsp$CO2e))
  # WD level
  n <- dim(res.tree)[1] # number of trees
  sp <- sum(res.tree$level=="species")
  gn <- sum(res.tree$level=="genus")
  s.u <- n - sp - gn
  WD.level <- data.frame(percentage=c(sp, gn, s.u)/n*100)
  row.names(WD.level) <- c("species", "genus", "sample unit")
  if(long) res <- list(tree = res.tree, taxon = AGBsp, total = AGBtot, WD.level=WD.level)
  else res <- list(taxon = AGBsp, total = AGBtot, WD.level=WD.level)
  class(res) <- "biomass"
  invisible(res)
}
