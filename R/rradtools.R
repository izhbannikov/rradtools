#rradtools.R
build_radsites <- function(fnames,keys,astable=T) 
{
    ans <- .Call("BuildRadSites",fnames,keys)
    
    
    if(astable) {
	out <- list()
	for(i in 1:length(ans)) {
		ans1 <- matrix(unlist(ans[[i]]), nrow = length(ans[[i]]), byrow = TRUE)
      		rownames(ans1) <- names(ans[[i]])
      		colnames(ans1) <- names(ans[[i]][[1]])
      		out[[i]] <- ans1
	}
	out
    } else {
      ans
    }
   
}
