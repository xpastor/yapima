extractCoords <- function(xx) 
{
  coords <- sapply(xx, strsplit, '[:-]')
  coords <- as.data.frame(do.call(rbind, coords), stringsAsFactors=F)
  colnames(coords) <- c('chrom', 'chromStart', 'chromEnd')
  return(coords)
}
