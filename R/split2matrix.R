split2matrix=function(x, split){
  #Split a character vector x using the defined split pattern and convert to a matrix
  #with number of rows equal to the length of x and number of columns equal to the
  #longest vector resulting from the split
  if(!inherits(x, "character"))
    stop("x must be a character array")
  
  s=strsplit(x, split=split)
  n=sapply(s, length)
  res=matrix(data=as.character(NA), nc=max(n), nr=length(s))
  index1=rep(1:length(x), n)
  index2=unlist(sapply(s, function(y){return(1:length(y))}, simplify=F))
  
  s=unlist(s)
  res[cbind(index1, index2)]=s
  return(res)
  
}