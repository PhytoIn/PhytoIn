BAplot <- function(formula, data, taxon="taxon", circumference=TRUE, quadrat.size, dead="dead", rm.dead=FALSE, 
                    origin=c(0, 0), col="grey40", alpha=1, cex.radius=1, ind.coord=FALSE, legend=TRUE, long=FALSE)
{
  if (!taxon %in% names(data)) stop("taxon is missing or misspelled")
  if(missing(quadrat.size)) stop("Quadrat size is mandatory")
  y<-is.na(data) | data == "" # search empty and NAs cells
  data <- data[!apply(y == TRUE, 1, all), !apply(y == TRUE, 2, all)]
  if(all((.packages())!= "ggplot2")) require(ggplot2)
  if(all((.packages())!= "ggforce")) require(ggforce)
  if(all((.packages())!= "packcircles")) require(packcircles)
  
  if (rm.dead) # remove dead trees
  {
    filter = tolower(data[[taxon]]) == tolower(dead) 
    data <- data[!filter, ] 
    cat ("\n", sum(filter), "dead individuals removed from the dataset \n")
  }
  fm <- as.formula(formula) # convert 'formula' in a formula object
  vars <- all.vars(fm)
  vars <- c(vars[1], sort(vars[2:3]))
  
  # Individual radius
  measure <- data[[vars[1]]]  # extract the measure column
  measure <- gsub(",", ".", measure)
  
  # split multiple trunks
  my_list<-strsplit(x=measure, split="+", fixed=TRUE)
  n.obs <- sapply(my_list, length)
  seq.max <- seq_len(max(n.obs))
  mat <- t(sapply(my_list, "[", i = seq.max))
  measure.table<-matrix(as.numeric(mat), ncol = max(seq.max))   # Convert to numeric matrix
  measure.table[is.na(measure.table)] <- 0
  
  if (circumference) AB <- (measure.table/100)^2/(4*pi) else AB <- (pi*(measure.table/100)^2)/4 
  ABi <- apply(AB, 1, sum)   # basal area of each tree
  r <- sqrt(ABi/pi)
  if(length (quadrat.size)==1) quadrat.size <- c(quadrat.size, quadrat.size)
  
  # Test whether xy is present in 'data'
  if(ind.coord) {
    # Grid
    grid.x <- seq(origin[1], ceiling(max(data[[vars[2]]])/quadrat.size[1]) * quadrat.size[1], by=quadrat.size[1])
    grid.y <- seq(origin[2], ceiling(max(data[[vars[3]]])/quadrat.size[2]) * quadrat.size[2], by=quadrat.size[2])
    res.plot <- data.frame(taxon=data[[taxon]], x=data[[vars[2]]], y=data[[vars[3]]], radius=r)    # Junta a coluna com o nome das especies
  }
  else {
    quadrat.label<-with(data, paste(x,y))
    quadrat.index<-unique(quadrat.label)
    data<-cbind(quadrat.label, data)
    
    r <- r + 0.5   # Individual radius + 0.5 to avoid circle overlap
    xi <- NULL
    yi <- NULL
    
    for(i in 1:dim(data)[1])
    {
      xi[i]<-runif(1, min = data[[vars[2]]] [i] + 0.3, max = data[[vars[2]]] [i] + quadrat.size[1] - 0.3) # sample a random position within the quadrat
      yi[i]<-runif(1, min = data[[vars[3]]] [i] + 0.3, max = data[[vars[3]]] [i] + quadrat.size[2] - 0.3) 
    }
    # repel the coordinates to avoid circle overlap
    data <- cbind(xi, yi, r, data)
    res.plot <- NULL   
    for(j in quadrat.index){
      data.sub <- subset(data, subset=quadrat.label == j) 
      res <- circleRepelLayout(x=data.sub, xlim = c(data.sub[[vars[2]]] [1], data.sub[[vars[2]]] [1] + quadrat.size[1]),
                               ylim = c(data.sub[[vars[3]]] [1], data.sub[[vars[3]]] [1] + quadrat.size[2]), xysizecols = 1:3, sizetype="radius")$layout
      res <- cbind(res, taxon=data.sub[[taxon]])      # bind the column with taxon names
      res.plot <- rbind(res.plot, res)                # add data of quadrat j
    }
    res.plot$radius <- res.plot$radius - 0.5        # Subtract 0.5 from the radius
    
    grid.x <- unique(data[[vars[2]]])
    grid.y <- unique(data[[vars[3]]])
    
    grid.x <- c(grid.x, tail(grid.x, n=1) + quadrat.size[1])
    grid.y <- c(grid.y, tail(grid.y, n=1) + quadrat.size[1])
    
  } # End of 'if()'
  
  # plot the map
  
  if(legend){
    map<-ggplot() +
      theme(panel.grid.minor = element_blank()) +
      geom_circle(aes(x0 = x, y0 = y, r = radius*cex.radius, fill = taxon, color=taxon), alpha=alpha, data = res.plot) +
      scale_x_continuous(breaks = grid.x) +
      scale_y_continuous(breaks = grid.y) +
      coord_fixed()
    print(map)
  } else {
    map<-ggplot() +
      theme(panel.grid.minor = element_blank()) +
      geom_circle(aes(x0 = x, y0 = y, r = radius*cex.radius), fill = col, color=NA, alpha=alpha, show.legend=FALSE, data = res.plot) +
      scale_x_continuous(breaks = grid.x) +
      scale_y_continuous(breaks = grid.y) +
      coord_fixed()
    print(map)
  }
  if(long) return(res.plot)
}