sf_prod <- function(n = NULL, m = NULL, cpi = NULL, single = NULL){
  
  sfDf = data.frame(m = m,
                    N = n)
  sfDf = sfDf[order(m),]
  
  sfTab <- data.frame(matrix(NA_real_, nrow(sfDf)+1,7))
  
  names(sfTab) <- c('massClass', #1
                    'abundance', #2
                    'N.loss', #3
                    'biomass', #4
                    'mass.at.loss', #5
                    'biomass.loss', #6
                    'multiply.by.no.sizes') #7
  
  # set 1) & 2) as mass and abundance
  sfTab[,1] <- c(sfDf$m, NA)
  sfTab[,2] <- c(sfDf$N, NA)

  # 3) add in the density loss column
  sfTab[,3] <- c(NA, (diff(sfDf[,2]) * -1), sfTab[(nrow(sfTab)-1),2])
 
    # 4) add in the biomass column
  sfTab[,4] <- sfTab[['massClass']] * sfTab[['abundance']]
  
  # 5) add in the mass at loss column
  sfTab[,5] <- exp(stats::filter(log(sfTab[['massClass']]), c(1,1)/2, sides = 1))
  sfTab[nrow(sfTab),5] <- sfTab[(nrow(sfTab)-1),1]
  
  # 6) add in biomass loss
  sfTab[,6] <- sfTab[['mass.at.loss']] * sfTab[['N.loss']]
  
  if(!single){
  # 7) multiply by # of size classes
  sfTab[,7] <- sfTab[['biomass.loss']] * (length(sfDf$m))
  } else{
    sfTab[,7] <- sfTab[['biomass.loss']]
  }
  # if first biomass.loss x size classes is negative, set to 0
  if(sfTab[2,'biomass.loss'] < 0) sfTab[2,'biomass.loss'] <-0
  
  P.uncorr.samp <- sum(sfTab[,7], na.rm = TRUE)
  
  P.ann.samp <- P.uncorr.samp * (365/cpi)
  
  return(P.ann.samp)
  
}