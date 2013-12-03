#rradtools.R
build_radsites <- function(fnames,keys,astable=T) 
{
    ans <- .Call("BuildRadSites",fnames,keys)
    
    
    if(astable) {
      ans1 <- matrix(unlist(ans), nrow = length(ans), byrow = TRUE)
      rownames(ans1) <- names(ans)
      colnames(ans1) <- names(ans[[1]])
      ans1
    } else {
      ans
    }
   
}