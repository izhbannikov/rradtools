#rradtools.R
build_radsites <- function(fnames,keys,astable=T,merge_sites=T) 
{
    ans <- .Call("BuildRadSites",fnames,keys,merge_sites)
    
    
    if(astable & merge_sites) {
	    out <- list()
	    for(i in 1:length(ans)) {
		      ans1 <- matrix(unlist(ans[[i]]), nrow = length(ans[[i]]), byrow = TRUE)
      		rownames(ans1) <- names(ans[[i]])
      		colnames(ans1) <- names(ans[[i]][[1]])
      		out[[i]] <- ans1
	    }
	    out
    } else if(astable & !merge_sites){
        ans1 <- matrix(unlist(ans), nrow = length(ans), byrow = TRUE)
        rownames(ans1) <- names(ans)
        colnames(ans1) <- names(ans[[1]])
        ans1
    }else if(!astable & merge_sites){
      ans
    } else if(!astable & !merge_sites) {
      ans
    }
   
}
