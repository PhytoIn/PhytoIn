phytoparam <- function(x, measure.label, h="h", taxon="taxon", family="family", dead="dead", circumference=TRUE, su="quadrat", height=TRUE, 
                        quadrat=TRUE, su.size, d="distance", shape.factor=1, rm.dead=FALSE, check.spelling=TRUE)
{
  if(!missing(measure.label)) measure.label <- tolower(measure.label)
  taxon <- tolower(taxon)
  family <- tolower(family)
  h <- tolower(h)
  dead <- tolower(dead)
  su <- tolower(su)
  d <- tolower(d)

  
  colnames(x) <- tolower(colnames(x)) # renomeia as colunas do data.frame usando apenas letras minusculas
  if (!family %in% names(x) & !is.na(family)) stop("family is missing or misspelled")
  if (!h %in% names(x) & height) stop("h is missing or misspelled")
  if (!measure.label %in% names(x)) stop("measure.label is missing or misspelled")
  if (!taxon %in% names(x)) stop("taxon is missing or misspelled")
  if (!su %in% names(x)) stop("su is missing or misspelled")
  if (!d %in% names(x) & !quadrat) stop("d is missing or misspelled")

  y<-is.na(x) | x == "" # busca todas as celulas vazias ou com NAs
  x <- x[!apply(y == TRUE, 1, all), !apply(y == TRUE, 2, all)]  # Remove todas as linhas/colunas vazias/NAs
  
  filter <- tolower(x[[taxon]]) == tolower(dead) # testa quais individuos sao "mortos"
  
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
  N.total <- aggregate(x[[taxon]]!="" ~ x[[taxon]], FUN = sum)  # calcula o N por especie incluindo as mortas
  
  if (rm.dead) # se a opcao for remover as plantas mortas do data frame 
  {
    x <- x[!filter, ] # remove linhas contendo plantas mortas do data frame
    cat ("\n", sum(filter), "dead individuals removed from the dataset \n") # exibe na tela uma mensagem com o n�mero de "mortas" removidas do data frame
  }
  x[[su]] <- factor(x[[su]]) # transforma a coluna parcelas em fator
  x[[taxon]] <- factor(x[[taxon]]) # transforma a coluna taxon em fator
  
  # *** QUADRAT METHOD ***       
  
  if(quadrat)
  {
    area <- length(unique(x[[su]])) * su.size / 10000  # sampled area
    
    #  *** DENSITY QUADRAT ***
    
    N <- aggregate(x[[taxon]]!="" ~ x[[taxon]], FUN = sum)  # calcula o N por esp�cie
    ADe <- N[,2]/area                                       # calcula a densidade absoluta
    RDe <- ADe/sum(ADe) * 100
    N.spp <- dim(N[tolower(N[, 1]) != tolower(dead), ])[1]

    #  *** FREQUENCY QUADRAT ***
    
    frequency.table <- table(x[[taxon]], x[[su]]) > 0
    AFr <- rowSums(frequency.table) / dim(frequency.table)[2] * 100
    RFr <- AFr/sum(AFr) * 100
    
    #  *** DOMINANCE QUADRAT ***
    
    if(missing(measure.label) & circumference==TRUE) measure.label = "cbh" # testa se a medida foi em CAP e atribui valor default para a variavel 
    if(missing(measure.label) & circumference==FALSE) measure.label = "dbh" # testa se a medida foi em DAP e atribui valor default para a variavel
    measure <- x[[measure.label]]  # extrai a coluna de valores das medidas
    measure <- gsub(",", ".", measure) # substitui eventuais v�gulas por ponto
    
    # Separa os fustes multiplos em colunas
    my_list<-strsplit(x=measure, split="+", fixed=TRUE)
#print(my_list)    
    n.obs <- sapply(my_list, length)
#print(n.obs)    
    seq.max <- seq_len(max(n.obs))
#print(seq.max)    
    mat <- t(sapply(my_list, "[", i = seq.max))
#print(mat)    
    measure.table<-matrix(as.numeric(mat), ncol = max(seq.max))    # Convert to numeric matrix
    measure.table[is.na(measure.table)] <- 0
#print(measure.table)    
    if (circumference) AB <- (measure.table/100)^2/(4*pi) else AB <- (pi*(measure.table/100)^2)/4  # calcula a �rea basal em m� de cada fuste a partir dos valores de cap ou dap
#print(AB)    
    
    x$ABi <- apply(AB, 1, sum) # calcula a area basal individual e junta � matriz de dados
#print(x$ABi)    
    G <- aggregate(x$ABi ~ x[[taxon]], FUN = sum)[,2] # calcula a area basal por especie
    ADo <- G/area         # Dominancia absolula
    tADo <- sum(ADo)      # Dominancia total
    RDo <- ADo/tADo*100   # Dominancia relativa
    
    table <- data.frame(Taxon = N[,1], N = N[,2], ADe = round(ADe, 2), RDe = round(RDe, 2), AFr= round(AFr, 2), RFr= round(RFr, 2), ADo = round(ADo, 2), RDo = round(RDo, 2)) # cria a tabela de par�metros
    
    #  *** VOLUME QUADRAT ***
    
    if(height)
    {
      Vi <- x$ABi * x[[h]] * shape.factor   # calcula o volume individual
      volume <- aggregate(Vi ~ x[[taxon]], FUN = sum)   # calcula o volume por esp�cie
      AVol <- volume[, 2]/area           # calcula o volume absoluto em m�/ha
      RVol <- AVol/sum(AVol) * 100       # calcula o volume relativo
      table <- cbind(table, AVol= round(AVol, 2), RVol= round(RVol, 2)) # Adiciona as colunas AVol e RVol
    }
    
    ############################################################################
    
  } else {   # fecha o bloco de parcelas e comeca o de ponto-quadrante
    
    if(su=="quadrat") su="point"
    if(missing(measure.label) & circumference==TRUE) measure.label = "cbh" # testa se a medida foi em CAP e atribui valor default para a variavel 
    if(missing(measure.label) & circumference==FALSE) measure.label = "dbh" # testa se a medida foi em DAP e atribui valor default para a variavel
    measure <- x[[measure.label]]  # extrai a coluna de valores das medidas
    measure <- gsub(",", ".", measure) # substitui eventuais v�gulas por ponto
    
    # Separa os fustes multiplos em colunas
    my_list<-strsplit(x=measure, split="+", fixed=TRUE)
    n.obs <- sapply(my_list, length)
    seq.max <- seq_len(max(n.obs))
    mat <- t(sapply(my_list, "[", i = seq.max))
    measure.table<-matrix(as.numeric(mat), ncol = ncol(mat))    # Convert to numeric matrix
    measure.table[is.na(measure.table)] <- 0

    if (circumference) AB <- (measure.table/100)^2/(4*pi) else AB <- (pi*(measure.table/100)^2)/4  # calcula a �rea basal em m� de cada fuste a partir dos valores de cap ou dap
    
    x$ABi <- apply(AB, 1, sum) # calcula a area basal individual e junta � matriz de dados

    N <- aggregate(x[[taxon]]!="" ~ x[[taxon]], FUN = sum)  # calcula o N por esp�cie
    
    N.spp <- dim(N[tolower(N[, 1]) != tolower(dead), ])[1]
    
    G <- aggregate(x$ABi ~ x[[taxon]], FUN = sum)[,2] # calcula a area basal por especie
    
    DC <- x[[d]]  + sqrt(x$ABi/pi)  # Dist�ncia corrigida em metros

    Ai <- DC^2 * pi / 4             # �rea ocupada por cada indiv�duo
    
    Ai.su <- aggregate(Ai ~ x[[su]], FUN = sum)  # Area occupied by a point
    
    AM <- sum(Ai)/(length(Ai) - 1)  # �rea m�dia ocupada por um individuo

    area <- sum(Ai)/10000  # �rea total amostrada em ha
    

    #  *** DENSITY POINT ***
    
    DT <- 10000/AM   # Densidade total em ha
    ADe <- N[, 2] * DT/sum(N[, 2])             # calcula a densidade absoluta
    RDe <- ADe/sum(ADe) * 100
    
    #  *** FREQUENCY POINT ***
    
    frequency.table <- table(x[[taxon]], x[[su]]) > 0
    AFr <- rowSums(frequency.table) / dim(frequency.table)[2] * 100
    RFr <- AFr/sum(AFr) * 100
    
    #  *** DOMINANCE POINT ***
    
    ABM <- mean(x$ABi)
    tADo<- DT * ABM
    ADo <- G * tADo / sum(x$ABi)   # Dominancia absolula
    tADo <- sum(ADo)               # Dominancia total
    RDo <- ADo/tADo*100            # Dominancia relativa
    
    table <- data.frame(Taxon = N[,1], N = N[,2], ADe = round(ADe, 2), RDe = round(RDe, 2), AFr = round(AFr, 2), RFr= round(RFr, 2), ADo = round(ADo, 2), RDo = round(RDo, 2)) # cria a tabela de par�metros
    
    #  *** VOLUME POINT ***
    
    if(height)
    {
      Vi <- x$ABi * x[[h]] * shape.factor   # calcula o volume individual
      volume <- aggregate(Vi ~ x[[taxon]], FUN = sum)   # calcula o volume por esp�cie
      AVol <- volume[, 2]/area           # calcula o volume absoluto em m�/ha
      RVol <- AVol/sum(AVol) * 100       # calcula o volume relativo
      table <- cbind(table, AVol= round(AVol, 2), RVol= round(RVol, 2)) # Adiciona as colunas AVol e RVol
    }
  } # fecha o bloco de quadrantes
  
  if(quadrat) vars <- list(taxon, h, dead, su, su.size) else vars <- list(taxon, h, dead, su, Ai.su)
  
  # *** Estatisticas por familias ***
  
  if(!is.na(family)) {
    if(rm.dead == F)     xf <- x[!filter, ] else xf <- x
    ind.fam<-table(xf[[family]])
    spp.fam <- table(xf[[family]], xf[[taxon]]) > 0
    spp.error <- colnames(spp.fam)[colSums(spp.fam)>1]
    if(sum(colSums(spp.fam)>1)) stop(paste("\n", spp.error, "is assigned to more than one family"))
    table.fam <- data.frame(Family=row.names(ind.fam), Indivuals=as.vector(ind.fam), Species=rowSums(spp.fam)) # Individuals per family
    row.names(table.fam) <- NULL
    N.fam<-dim(table.fam)[1]    # numero de familias
  }
  
  #  *** IV and CV ***
  
  IV <- RDe + RFr + RDo  # calcula o �ndice de valor de import�ncia (IVI)
  CV <- RDe + RDo        # calcula o �ndice de valor de cobertura (IVC)
  
  table <- cbind(table, IV = round(IV, 2), CV = round(CV, 2)) # Adiciona as coluna IV e CV
  N.taxon<-dim(table)[1]  # Numero de especies
  
  #  *** Shannon, Simpson and Pielou indeces ***
  
  Pe.dead <- N.total[, 2]/sum(N.total[,2]) # densidade relativa incluindo as mortas
  H.dead <- -sum(Pe.dead*log(Pe.dead))     # H' incluindo as mortas
  J.dead <- H.dead/log(dim(N.total)[1])    # J incluindo as mortas
  C.dead <- 1 - (sum(N.total[, 2] * (N.total[, 2] - 1))/(sum(N.total[,2])*(sum(N.total[,2]) - 1)))  # indice de Simpson incluindo as mortas
  
  if(!rm.dead){                                         # Teste para remover os indiv�duos mortos
    #    filter.dead = tolower(x[[taxon]]) == tolower(dead)       # transforma para letra min�scula todas as palavras "morta" no objeto dados 
    x1 <- x[!filter, ]
    N1 <- aggregate(x1[[taxon]]!="" ~ x1[[taxon]], FUN = sum)  # calcula o N por esp�cie
  } else {
    x1 <- x  # remove linhas contendo plantas mortas
    N1 <- N
  }
  
  Pe <- N1[, 2]/sum(N1[,2]) # densidade relativa excluindo as mortas
  H <- -sum(Pe*log(Pe))     # H' excluindo as mortas
  J <- H/log(dim(N1)[1])    # J excluindo as mortas
  C <- 1 - (sum(N1[, 2] * (N1[, 2] - 1))/(sum(N1[,2])*(sum(N1[,2]) - 1)))  # indice de Simpson excluindo as mortas
  
  if(quadrat){  # prepara a tabela com os parametros totais
    if(!is.na(family)) {
      global.var <- c("N. of individuals", "N. of samples", "Sampled area", "Total density", "Total dominace", "N. of families", "N. of species",
                      "Shannon-Wiener with dead", "Shannon-Wiener excluding dead", "Pielou with dead", "Pielou excluding dead", 
                      "Simpson with dead", "Simpson excluding dead")
      global.values <- round(c(sum(N[, 2]), dim(frequency.table)[2], area, sum(ADe), sum(ADo), N.fam, N.spp, H.dead, H,  J.dead, J, C.dead, C), 4)
      global <- data.frame(Parameter = global.var, Value = global.values)
    }
    else
    {
      global.var <- c("N. of individuals", "N. of samples", "Sampled area", "Total density", "Total dominace", "N. of species",
                      "Shannon-Wiener with dead", "Shannon-Wiener excluding dead", "Pielou with dead", "Pielou excluding dead", 
                      "Simpson with dead", "Simpson excluding dead")
      global.values <- round(c(sum(N[, 2]), dim(frequency.table)[2], area, sum(ADe), sum(ADo), N.spp, H.dead, H,  J.dead, J, C.dead, C), 4)
      global <- data.frame(Parameter = global.var, Value = global.values)
    }
    
  } else 
  {
    if(!is.na(family)) {
      global.var <- c("N. of individuals", "N. of samples", "Sampled area", "Unbiased total density", "Total dominace", "N. of families", 
                      "N. of species", "Shannon-Wiener with dead", "Shannon-Wiener excluding dead", "Pielou with dead", "Pielou excluding dead",
                      "Simpson with dead", "Simpson excluding dead")
      global.values <- round(c(sum(N[, 2]), dim(frequency.table)[2], area, DT, tADo, N.fam, N.spp, H.dead, H,  J.dead, J, C.dead, C), 4)
      global <- data.frame(Parameter = global.var, Value = global.values)
    }
    else
    {
      global.var <- c("N. of individuals", "N. of samples", "Sampled area", "Unbiased total density", "Total dominace",  
                      "N. of species", "Shannon-Wiener with dead", "Shannon-Wiener excluding dead", "Pielou with dead", "Pielou excluding dead",
                      "Simpson with dead", "Simpson excluding dead")
      global.values <- round(c(sum(N[, 2]), dim(frequency.table)[2], area, DT, tADo, N.spp, H.dead, H,  J.dead, J, C.dead, C), 4)
      global <- data.frame(Parameter = global.var, Value = global.values) 
    }
  }
  table <- table[order(table$IV, decreasing=T), ]
  rownames(table) <- seq(1, dim(table)[1])
  if(!is.na(family)){
    result <- list(vars = vars, data = x, global = global, family = table.fam, param = table)
    class(result) <- "param"
    invisible(result)
  }
  else
  {
    result <- list(vars = vars, data = x, global = global, param = table)
    class(result) <- "param"
    invisible(result)
  }
}


summary.param <- function(obj) {
  rownames(obj$global)<-NULL
  if("N. of families" %in% obj$global[,1]) print(obj$global[1:7,])
  else print(obj$global[1:6,])
}
plot.param<-function(obj){
  param<-obj$param
  param <- param[order(param$IV, decreasing=F), ]
  mar<-max(nchar(as.character(param$Taxon)))*0.46
  par(mar=c(6, mar, 4, 2)+0.1)
  barplot(cbind(param$RDe, param$RFr, param$RDo) ~ Taxon, data=param,
          border=NA, las=1,
          #args.legend=c(x="bottomright", bty="n"),
          #legend.text=c("Rel. Density", "Rel. Frequency", "Rel. Dominance"),
          ylab=NA, xlab = "IV", beside=F, 
          horiz=T)
  grid(lty = 1, col = "gray", lwd = 1)
  barplot(cbind(param$RDe, param$RFr, param$RDo) ~ Taxon, data=param,
          border=NA, las=1,
          args.legend=c(x="bottomright", bty="n", border=NA),
          legend.text=c("Rel. Density", "Rel. Frequency", "Rel. Dominance"),
          ylab=NA, beside=F, 
          horiz=T, add=T)
}

