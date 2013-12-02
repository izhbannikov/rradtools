#rradtools.R
build_radsites <- function(fnames,keys) 
{
   .Call("BuildRadSites",fnames,keys)
}