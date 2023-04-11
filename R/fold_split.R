#' Function for creating K-folds for cross-validation
#' 
#' Internal function. Do not use externally
#' 
#' @param K number of folds
#' @index id for observations

fold_split <- function(K = 3, 
                       index = NULL){
  data_indices <- data.frame(cbind(id=1:length(index),
                        zip_id=index))
  ziplist <- unique(index)
  sapply(ziplist, function(x){
    temp_data <- data_indices[!is.na(match(data_indices[,2],x)),1]
    n <- length(temp_data)
    kn <- floor(n/K)
    setn <- temp_data
    id <- list()
    i <- 1
    while(i<K){
      idx <- sample(setn,kn, replace=F)
      id[[i]] <- idx
      new.setn = setn[! setn %in% idx]
      setn = new.setn
      i <- i+1
    }
    id[[K]] <- setn
    id
  },simplify = F)
}

