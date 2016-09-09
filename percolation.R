

library(osmar)
library(rgdal)
library(rgeos)
library(igraph)
library(ggplot2)

# load data, convert to spatial polygons
start=proc.time()[3]
buildings = as_sp(get_osm(corner_bbox(0.0, 45.0,20.0, 50.0),source=osmsource_file('buildings_lux.osm')))$polygons
gc(verbose=TRUE)
show(start-proc.time()[3])
# for luxembourg -> 13703s cpu time

##
# Function to construct topological network from spatial polygons
#   returns igraph object
constructNetwork <- function(polygons){
  # get building centroids
  centroids = gCentroid(polygons,byid=TRUE)
  #plot(centroids,col="red",add=TRUE)
  # non-unique centroids/polygons sometimes (2 for lux)
  centroids = SpatialPoints(unique(centroids@coords),proj4string = centroids@proj4string)
  
  # delaunay triangulation
  delaunay=gDelaunayTriangulation(centroids)
  #plot(delaunay,col="blue",add=T)
  
  # construction of the graph
  # let do it dirtily taking coordinates as vertices id
  from=c();to=c();elengths=c();vids=c();vxcors=c();vycors=c()
  n=length(delaunay@polygons)
  for(p in 1:n){
    if(p%%100==0){show(100*p/n)}
    polygon = delaunay@polygons[[p]]@Polygons[[1]]
    for(i in 2:nrow(polygon@coords)){
      coords=polygon@coords
      # polygons are assumed closed, can add vertices only once
      v1=paste0(coords[i-1,1],"-",coords[i-1,2])
      v2=paste0(coords[i,1],"-",coords[i,2])
      # take euclidian dist for length, not assumed to be used on long distances
      #elength = sqrt((coords[i-1,1]-coords[i,1])^2+(coords[i-1,2]-coords[i,2])^2)
      # better use spDists to have in kms
      elength=spDistsN1(matrix(polygon@coords[i,],ncol=2),polygon@coords[i-1,],longlat = TRUE)
      elengths=append(elengths,elength)
      from=append(from,v1);to=append(to,v2)
      vids=append(vids,v2);vxcors=append(vxcors,coords[i,1]);vycors=append(vycors,coords[i,2])
    }
  }
  edf=data.frame(from=from,to=to,length=elengths);vdf=unique(data.frame(id=vids,x=vxcors,y=vycors))
  
  g = graph_from_data_frame(d = edf,directed = FALSE)
  #V(g)[vdf[,1]]$x=vdf[,2];V(g)[vdf[,1]]$y=vdf[,3]
  # bug/memoryleak in igraph when trying to index this way ?
  #  -> do not need vertices coords in a first time if edges have length
  
  return(g)
}


##
# Function constructing gravity based graph with threshold filtering
#
gravityGraph<-function(g,beta,threshold=0,quant=NULL){
  # compute gravity potential
  gg=g
  E(gg)$t=1/((E(gg)$length)^beta)
  
  if(!is.null(quant)){
    threshold = quantile(E(gg)$t,quant)
  }
  
  # filter
  gg=subgraph.edges(gg,which(E(gg)$t>threshold),delete.vertices = FALSE)
  
  return(gg)
}


##
#  Indicator "link quantity" per cluster
#  Q_g = log(\sum t_ij)
linkQuantity<-function(gg){
  comps=components(gg)
  res=c()
  for(c in 1:length(comps$csize)){
    v = V(gg)[comps$membership==c]
    e = E(gg)[v%--%v]
    res=append(res,log(sum(e$t)))
  }
  return(res)
}




###
# Tests


# construct network

g=constructNetwork(buildings)

# test
gg=gravityGraph(g,beta=1.5,q=0.5)




# plot numcomp = f(quantile) at fixed betas
res=data.frame()
for(q in (1:100)/100){
show(q)
gg=gravityGraph(g,beta=1.5,quant=q)
res=rbind(res,c(q,length(components(gg)$csize)))
}
colnames(res)=c("q","comps")

gp=ggplot(res,aes(x=q,y=comps))
gp+geom_line()



# plot rank-size of link quantity for varying beta (paper plots)

q = 0.5
res=data.frame()
for(beta in seq(from=0.5,to=4.0,by=0.1)){
  for(q in c(0.2,0.4,0.6,0.8)){
  show(beta)
  gg=gravityGraph(g,beta=beta,quant=q)
  lq = linkQuantity(gg)
  lq[is.infinite(lq)]=0
  res=rbind(res,data.frame(rep(beta,length(lq)),rep(q,length(lq)),sort(lq,decreasing=TRUE),1:length(lq)))
  }
}
colnames(res)=c("beta","q","lq","g")

gp=ggplot(res,aes(x=g,y=lq,colour=beta,group=beta))
gp+geom_line()+facet_wrap(~q)






