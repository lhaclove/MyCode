# Library
library(ggplot2)

# mtcars data
head(mtcars)

# First type of color
ggplot(mtcars, aes(factor(cyl), mpg)) + 
  geom_violin(aes(fill = cyl))

# Second type
ggplot(mtcars, aes(factor(cyl), mpg)) +
  geom_violin(aes(fill = factor(cyl)))

#################Charge the vioplot library
library(vioplot)

# Create data
treatment=c(rep("A", 40) , rep("B", 40) , rep("C", 40) )
value=c( sample(2:5, 40 , replace=T) , sample(c(1:5,12:17), 40 , replace=T), sample(1:7, 40 , replace=T) )
data=data.frame(treatment,value)

# Draw the plot
with(data , vioplot( value[treatment=="A"] , value[treatment=="B"], value[treatment=="C"],  col=rgb(0.1,0.4,0.7,0.7) , names=c("A","B","C") ))


###################
# For the weatherAUS dataset.
library(rattle)

# To generate a density plot.
library(ggplot2)  
cities <- c("Canberra", "Darwin", "Melbourne", "Sydney")
ds <- subset(weatherAUS, Location %in% cities & ! is.na(Temp3pm))
p  <- ggplot(ds, aes(Temp3pm, colour=Location, fill=Location))
p  <- p + geom_density(alpha=0.55)
p


#############################
# ggplot2 library
library(ggplot2)

# Let's use the diamonds dataset
data(diamonds)
head(diamonds)

# plot 1: Density of price for each type of cut of the diamond:
ggplot(data=diamonds,aes(x=price, group=cut, fill=cut)) + 
  geom_density(adjust=1.5)

# plot 2: Density plot with transparency (using the alpha argument):
ggplot(data=diamonds,aes(x=price, group=cut, fill=cut)) + 
  geom_density(adjust=1.5 , alpha=0.2)

# plot 3: Stacked density plot:
ggplot(data=diamonds,aes(x=price, group=cut, fill=cut)) + 
  geom_density(adjust=1.5, position="fill")

# plot 4
ggplot(diamonds, aes(x=depth, y=..density..)) + 
  geom_density(aes(fill=cut), position="stack") +
  xlim(50,75) + 
  theme(legend.position="none")


#################

# Library
library(tidyverse)
library(plotly)

# A classic histogram for the iris data set (left)
ggplot(iris, aes(x=Sepal.Length)) +
  geom_histogram()

# Transform a litte bit the dataset to make dots
don = iris %>% 
  arrange(Sepal.Length) %>% # sort using the numeric variable that interest you
  mutate(var_rounded = (Sepal.Length+1) - ( (Sepal.Length+1) %% 0.2 ) ) %>% # This attributes a bin to each observation. Here 0.2 is the size of the bin.
  mutate(y=ave(var_rounded, var_rounded, FUN=seq_along)) # This calculates the position on the Y axis: 1, 2, 3, 4...

# Make the plot (middle)
ggplot(don, aes(x=var_rounded, y=y) ) +
  geom_point( size=6, color="skyblue" ) 

# Improve the plot, and make it interactive (right)
don=don %>% mutate(text=paste("ID: ", rownames(iris), "\n", "Sepal Length: ", Sepal.Length, "\n", "Species:: ", Species, sep="" )) 
p=ggplot(don, aes(x=var_rounded, y=y) ) +
  geom_point( aes(text=text), size=6, color="skyblue" ) +
  xlab('Sepal Length') +
  ylab('# of individual') +
  theme_classic() +
  theme(
    legend.position="none",
    axis.line.y = element_blank(),
    axis.text=element_text(size=15)
  )
p

# Use the magic of ggplotly to have an interactive version
ggplotly(p, tooltip="text")


##################################

# Library
library(igraph)

# Create data
set.seed(1)
data=matrix(sample(0:1, 100, replace=TRUE, prob=c(0.8,0.2)), nc=10)
network=graph_from_adjacency_matrix(data , mode='undirected', diag=F )

# Default network
par(mar=c(0,0,0,0))
plot(network)


#########################
# library
library(igraph)

# Create data
data=matrix(sample(0:1, 400, replace=TRUE, prob=c(0.8,0.2)), nrow=20)
network=graph_from_adjacency_matrix(data , mode='undirected', diag=F )

# When ploting, we can use different layouts:
par(mfrow=c(2,2), mar=c(1,1,1,1))
plot(network, layout=layout.sphere, main="sphere")
plot(network, layout=layout.circle, main="circle")
plot(network, layout=layout.random, main="random")
plot(network, layout=layout.fruchterman.reingold, main="fruchterman.reingold")

# See the complete list with
help(layout)
##########################


# library
library(igraph)

# create data:
links=data.frame(
  source=c("A","A", "A", "A", "A","J", "B", "B", "C", "C", "D","I"),
  target=c("B","B", "C", "D", "J","A","E", "F", "G", "H", "I","I")
)

# Turn it into igraph object
network=graph_from_data_frame(d=links, directed=F) 

# Count the number of degree for each node:
deg=degree(network, mode="all")

# Plot
plot(network, vertex.size=deg*6, vertex.color=rgb(0.1,0.7,0.8,0.5) )

############
# libraries
library(networkD3)

# Load data
data(MisLinks)
data(MisNodes)

# Plot
forceNetwork(
  Links = MisLinks, Nodes = MisNodes, 
  Source = "source", Target = "target",
  Value = "value", NodeID = "name",
  Group = "group", opacity = 0.8,
  linkDistance = JS('function(){d3.select("body").style("background-color", "#DAE3F9"); return 50;}')
)


##################
# Create data with the randomNames package :
library('randomNames') 
NUMOFLINKS = 100
relations = data.frame(source = randomNames(1000,which.names='both'), target = "")
relations = relations[rep(seq_len(nrow(relations)), sample(1:10,nrow(relations), replace=T)),]
relations = relations[sample(nrow(relations),NUMOFLINKS),] 
relations$target = sample(relations$source,nrow(relations), replace = T)
relations = relations[relations[,1]!=relations[,2], ] 

# Have a look to the input table !
head(relations)

## Plot the graph using the IGRAPH package
library("igraph")
vertices<-data.frame("name" = unique(unlist(relations))) # node names
g = graph.data.frame(relations, directed=F, vertices=vertices) # raw graph
vertices$group = edge.betweenness.community(g)$membership # betweeness centrality for each node for grouping


plot(g,
     #mark.groups=vertices$group, # group vertices by betweeness indicator (redish blob background)
     layout=layout.auto, 
     vertex.color = vertices$group, # color vertices by edge betweeness
     vertex.label=NA, # no vertex label (name)
     vertex.size=5,
     edge.arrow.size=0.8)

##############################################
### You need several libraries
library(circlize)
library(migest)
library(dplyr)

### Make data
m <- data.frame(order = 1:6,
                country = c("Ausralia", "India", "China", "Japan", "Thailand", "Malaysia"),
                V3 = c(1, 150000, 90000, 180000, 15000, 10000),
                V4 = c(35000, 1, 10000, 12000, 25000, 8000),
                V5 = c(10000, 7000, 1, 40000, 5000, 4000),
                V6 = c(7000, 8000, 175000, 1, 11000, 18000),
                V7 = c(70000, 30000, 22000, 120000, 1, 40000),
                V8 = c(60000, 90000, 110000, 14000, 30000, 1),
                r = c(255,255,255,153,51,51),
                g = c(51, 153, 255, 255, 255, 255),
                b = c(51, 51, 51, 51, 51, 153),
                stringsAsFactors = FALSE)
df1 <- m[, c(1,2, 9:11)]
m <- m[,-(1:2)]/1e04
m <- as.matrix(m[,c(1:6)])
dimnames(m) <- list(orig = df1$country, dest = df1$country)
#Sort order of data.frame and matrix for plotting in circos
df1 <- arrange(df1, order)
df1$country <- factor(df1$country, levels = df1$country)
m <- m[levels(df1$country),levels(df1$country)]


### Define ranges of circos sectors and their colors (both of the sectors and the links)
df1$xmin <- 0
df1$xmax <- rowSums(m) + colSums(m)
n <- nrow(df1)
df1$rcol<-rgb(df1$r, df1$g, df1$b, max = 255)
df1$lcol<-rgb(df1$r, df1$g, df1$b, alpha=200, max = 255)

### Plot sectors (outer part)
par(mar=rep(0,4))
circos.clear()

### Basic circos graphic parameters
circos.par(cell.padding=c(0,0,0,0), track.margin=c(0,0.15), start.degree = 90, gap.degree =4)

### Sector details
circos.initialize(factors = df1$country, xlim = cbind(df1$xmin, df1$xmax))

### Plot sectors
circos.trackPlotRegion(ylim = c(0, 1), factors = df1$country, track.height=0.1,
                       #panel.fun for each sector
                       panel.fun = function(x, y) {
                         #select details of current sector
                         name = get.cell.meta.data("sector.index")
                         i = get.cell.meta.data("sector.numeric.index")
                         xlim = get.cell.meta.data("xlim")
                         ylim = get.cell.meta.data("ylim")
                         
                         #text direction (dd) and adjusmtents (aa)
                         theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
                         dd <- ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
                         aa = c(1, 0.5)
                         if(theta < 90 || theta > 270)  aa = c(0, 0.5)
                         
                         #plot country labels
                         circos.text(x=mean(xlim), y=1.7, labels=name, facing = dd, cex=0.6,  adj = aa)
                         
                         #plot main sector
                         circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2], ytop=ylim[2], 
                                     col = df1$rcol[i], border=df1$rcol[i])
                         
                         #blank in part of main sector
                         circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2]-rowSums(m)[i], ytop=ylim[1]+0.3, 
                                     col = "white", border = "white")
                         
                         #white line all the way around
                         circos.rect(xleft=xlim[1], ybottom=0.3, xright=xlim[2], ytop=0.32, col = "white", border = "white")
                         
                         #plot axis
                         circos.axis(labels.cex=0.6, direction = "outside", major.at=seq(from=0,to=floor(df1$xmax)[i],by=5), 
                                     minor.ticks=1, labels.away.percentage = 0.15)
                       })

### Plot links (inner part)
### Add sum values to df1, marking the x-position of the first links
### out (sum1) and in (sum2). Updated for further links in loop below.
df1$sum1 <- colSums(m)
df1$sum2 <- numeric(n)

### Create a data.frame of the flow matrix sorted by flow size, to allow largest flow plotted first
df2 <- cbind(as.data.frame(m),orig=rownames(m),  stringsAsFactors=FALSE)
df2 <- reshape(df2, idvar="orig", varying=list(1:n), direction="long",
               timevar="dest", time=rownames(m),  v.names = "m")
df2 <- arrange(df2,desc(m))

### Keep only the largest flows to avoid clutter
df2 <- subset(df2, m > quantile(m,0.6))

### Plot links
for(k in 1:nrow(df2)){
  #i,j reference of flow matrix
  i<-match(df2$orig[k],df1$country)
  j<-match(df2$dest[k],df1$country)
  
  #plot link
  circos.link(sector.index1=df1$country[i], point1=c(df1$sum1[i], df1$sum1[i] + abs(m[i, j])),
              sector.index2=df1$country[j], point2=c(df1$sum2[j], df1$sum2[j] + abs(m[i, j])),
              col = df1$lcol[i])
  
  #update sum1 and sum2 for use when plotting the next link
  df1$sum1[i] = df1$sum1[i] + abs(m[i, j])
  df1$sum2[j] = df1$sum2[j] + abs(m[i, j])
}

#######################

#Create data
name=c(3,10,10,3,6,7,8,3,6,1,2,2,6,10,2,3,3,10,4,5,9,10)
feature=paste("feature ", c(1,1,2,2,2,2,2,3,3,3,3,3,3,3,4,4,4,4,5,5,5,5) , sep="")
dat <- data.frame(name,feature)
dat <- with(dat, table(name, feature))

# Charge the circlize library
library(circlize)

# Make the circular plot
chordDiagram(as.data.frame(dat), transparency = 0.5)







