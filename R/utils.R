
last <- plant:::last

category_to_logical <- function(x, trueval) {
  x[x == ""] <- NA
  x == trueval
}

change_column_names <- function(dat, table) {
  i <- match(names(dat), table$var_in)
  j <- !is.na(i) & !is.na(table$var_out[i])
  names(dat)[j] <- table$var_out[i[j]]
  dat[,j, drop=FALSE]
}

change_column_names_file <- function(dat, table_file) {
  change_column_names(dat, read.csv(table_file, stringsAsFactors=FALSE))
}

render_md_as_html <- function(filename) {
  rmarkdown::render(filename, "html_document", quiet=TRUE)
}


# Returns up to 80 unique, nice colors, generated using
# http://tools.medialab.sciences-po.fr/iwanthue/
# Starts repeating after 80
nice_colors<-function(n=80){
  cols<-rep(c("#75954F","#D455E9","#E34423","#4CAAE1","#5DE737","#DC9B94",
    "#DC3788","#E0A732","#67D4C1","#5F75E2","#1A3125","#65E689","#A8313C",
    "#8D6F96","#5F3819","#D8CFE4","#BDE640","#DAD799","#D981DD","#61AD34",
    "#B8784B","#892870","#445662","#493670","#3CA374","#E56C7F","#5F978F",
    "#BAE684","#DB732A","#7148A8","#867927","#918C68","#98A730","#DDA5D2",
    "#456C9C","#2B5024","#E4D742","#D3CAB6","#946661","#9B66E3","#AA3BA2",
    "#A98FE1","#9AD3E8","#5F8FE0","#DF3565","#D5AC81","#6AE4AE","#652326",
    "#575640","#2D6659","#26294A","#DA66AB","#E24849","#4A58A3","#9F3A59",
    "#71E764","#CF7A99","#3B7A24","#AA9FA9","#DD39C0","#604458","#C7C568",
    "#98A6DA","#DDAB5F","#96341B","#AED9A8","#55DBE7","#57B15C","#B9E0D5",
    "#638294","#D16F5E","#504E1A","#342724","#64916A","#975EA8","#9D641E",
    "#59A2BB","#7A3660","#64C32A", "#451431"),
            ceiling(n/80))
  cols[seq_len(n)]
}

make.transparent <- function(col, opacity=0.5) {
  tmp <- col2rgb(col)/255
  rgb(tmp[1,], tmp[2,], tmp[3,], alpha=opacity)
}

## Position label at a fractional x/y position on a plot
label <- function(px, py, lab, ..., adj=c(0, 1), log.y=FALSE, log.x=FALSE) {
  u <- par("usr")
  x <- u[1] + px*(u[2]-u[1])
  y <- u[3] + py*(u[4]-u[3])

  if(log.x)
    x <- 10^x

  if(log.y)
    y <- 10^y

  text(x,y,lab, adj=adj,...)
}
