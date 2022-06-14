##############################
##
## SECOND-DEGREE REWIRING MODEL
##
## AUTHOR: DAMIEN R. FARINE (dfarine@ab.mpg.de)
##
## REVISION DATE: 28 SEPTEMBER 2019
##
## PLEASE CITE: Farine, D.R. "Structural trade-offs can predict rewiring in shrinking social networks"
##
##############################



#### OVERVIEW OF THE CODE ####
#
# Section 1 provides general functions, including the 'shrinking network function' remove.network.node
#
# Section 2 provides code to simulate data as per the manuscript
#
# Section 3 provides plotting functions
#
# Section 4 provides calls to the plotting function for the different figures contained in the manuscript
##############################


##############################

#### Packages required

library(sna)
library(assortnet)
library(fields)

##############################


##############################

#### Section 1 — GENERAL FUNCTIONS ####

# inds is a dataframe containing N rows and traits is a vector of length N with trait values
# this function is an implementation of the model by Ilany & Ackay 2016 Nature Communications

update.network <- function(network, traits, pm, ps, pb=1, mu.rate=0) {

	# Calculate pop size
	N <- length(traits)

	# first delete a node at random
	del <- sample(1:N, 1)
	traits <- traits[-del]
	network <- network[-del,]
	network <- network[,-del]
	N <- N-1

	# now select a node to reproduce
	repro <- sample(1:N, 1)
	
	# add new trait value
	traits <- c(traits, traits[repro]+rnorm(1,mean=0,sd=mu.rate))
	N <- N + 1

	# calculate new relationships:
	rel <- rep(0, N)
	
	# with parent's associates
	edges <- which(network[repro,]==1)
	rel[edges] <- sample(c(0,1),length(edges),prob=c(1-pm,pm),replace=TRUE)
	# with non-associates
	no.edges <- which(network[repro,]==0)
	rel[no.edges] <- sample(c(0,1),length(no.edges),prob=c(1-ps,ps),replace=TRUE)

	# now over-write with parent
	rel[repro] <- sample(c(0,1),1,prob=c(1-pb,pb))

	# update network
	network <- rbind(network, rel[-N])
	network <- cbind(network, rel)	
	colnames(network) <- NULL

	return(list(network=network,traits=traits))

}


# THE function for doing the network rewiring. 
# requires: 
# network = an NxN network of edge connectivity
# traits = an N vector of individual traits (simply to delete the individual from the trait vector)
# pm = value of Pm
# ps = value of Ps
# id = the identity of the individual to remove (if NULL, then selects one at random)
# n.removed = the number of individuals to remove

remove.network.node <- function(network, traits, pm, ps, id=NULL, n.removed=1) {

	# Calculate pop size
	N <- length(traits)

	if (is.null(id)) {
		del <- sample(1:N,n.removed,replace=FALSE)
	} else {
		del <- id
	}

	# first capture edges present
	edges <- network[del,-del]

	# delete the node
	traits <- traits[-del]
	network <- network[-del,]
	network <- network[,-del]
	N <- N-1

	# now update the network
	# first get potential edges
	potentials <- which(network[which(edges==1),which(edges==1),drop=FALSE]==0,arr.ind=T)
	
	# then allocate a new edge vs not
	if (length(potentials) > 0) {
		potentials <- potentials[which(potentials[,1] < potentials[,2]),,drop=FALSE]  # to avoid duplicates
		new.edge <- sample(c(0,1),nrow(potentials),prob=c(1-pm,pm),replace=T)
		network[which(edges==1),which(edges==1)][cbind(potentials[,1],potentials[,2])] <- new.edge
		network[which(edges==1),which(edges==1)][cbind(potentials[,2],potentials[,1])] <- new.edge
	}

	# then randomly allocate edges between newly disconnected nodes and other nodes
	potentials <- which(network[which(edges==1),which(edges==0),drop=FALSE]==0,arr.ind=T)
	if (length(potentials) > 0) {
		new.edge <- sample(c(0,1),nrow(potentials),prob=c(1-ps,ps),replace=T)
		network[which(edges==1),which(edges==0)][cbind(potentials[,1],potentials[,2])] <- new.edge
		network[which(edges==0),which(edges==1)][cbind(potentials[,2],potentials[,1])] <- new.edge
	}

	# return
	return(list(network=network,traits=traits))
}

##############################

##############################

### Section 2 - CODE TO GENERATE NETWORKS AND SIMULATE SHRINKING PROCESSES

# Code below creates an R data file containing the results for each of the set of parameters

## Model parameter(s)
for (pm in c(0.2)) {
for (ps in c(0)) {

# Network parameters
N <- 50			# Nodes in the network
nodes.removed <- N/2 	# Nodes to remove
n.removed <- 1		# How many to remove at a time

# Simulation parameters
n.rep <- 2500
n.rep.rem <- 100

# Results storage
den.orig <- matrix(NA, n.rep, n.rep.rem)
den.shrunk <- matrix(NA, n.rep, n.rep.rem)
mean.deg.orig <- matrix(NA, n.rep, n.rep.rem)
mean.deg.shrunk <- matrix(NA, n.rep, n.rep.rem)
assort.orig <- matrix(NA, n.rep, n.rep.rem)
assort.shrunk <- matrix(NA, n.rep, n.rep.rem)
clust.orig <- matrix(NA, n.rep, n.rep.rem)
clust.shrunk <- matrix(NA, n.rep, n.rep.rem)
params.in <- data.frame(pm.in = rep(NA, n.rep), ps.in = rep(NA, n.rep), pb.in = rep(NA, n.rep), mu.rate = rep(NA, n.rep))

# Run replications
for (zz in 1:n.rep) {

	# Generate a random starting network
	tie.prob <- runif(1)
	network.orig <-  rgraph(N,tprob=tie.prob,mode="graph")

	# Run the Ilany & Ackay model
	traits.orig <- rnorm(N)
	params.in$pm.in[zz] <- runif(1,0.5,1)
	params.in$ps.in[zz] <- runif(1,0,0.05)
	params.in$pb.in[zz] <- runif(1,0.95,1)
	params.in$mu.rate <- runif(1,0,0.4)	

	burn.in <- 500
	for (i in 1:burn.in) {
		output <- update.network(network.orig, traits.orig, params.in$pm.in[zz], params.in$ps.in[zz], params.in$pb.in[zz], params.in$mu.rate)
		network.orig <- output$network
		traits.orig <- output$traits
	}

	assort.orig[zz,] <- assortment.continuous(network.orig, traits.orig, weighted=FALSE)$r
	den.orig[zz,] <- gden(network.orig, mode="graph")
	mean.deg.orig[zz,] <- mean(degree(network.orig, gmode="graph",ignore.eval=TRUE))
	clust.orig[zz,] <- gtrans(network.orig,mode="graph")

	# For computational efficiency - run multiple removals for each network
	for (zz2 in 1:n.rep.rem) {

		# temp storage
		traits1 <- traits.orig
		network1 <- network.orig

		for (i in 1:nodes.removed) {
			id <- sample(1:(length(traits1)),1)
			output1 <- remove.network.node(network1, traits1, pm, ps, id)
			network1 <- output1$network
			traits1 <- output1$traits
		}

		assort.shrunk[zz,zz2] <- assortment.continuous(network1, traits1, weighted=FALSE)$r
		den.shrunk[zz,zz2] <- gden(network1, mode="graph")
		mean.deg.shrunk[zz,zz2] <- mean(degree(network1, gmode="graph",ignore.eval=TRUE))
		clust.shrunk[zz,zz2] <- gtrans(network1,mode="graph")

	}
	#if (zz %% 100 == 0) {
	#	cat(paste("pm =",pm,", ps =",ps,", ",zz,"of",nrow(params.in),"\n",sep=" "))
	#}
}

if (n.removed == 1) {
	save(list=ls(),file=paste("pm_",pm,"_ps_",ps,"_N_",N,".RData",sep=""))
} else {
	save(list=ls(),file=paste("pm_",pm,"_ps_",ps,"_N_",N,"_removed_",n.removed,".RData",sep=""))
}

}
}

##############################

##############################

#### Section 3 - PLOTTING FUNCTIONS  ####


right.half.circle <- function(x,y,r,nsteps=100,...){  
  rs <- seq(-pi/2,pi/2,len=nsteps) 
  xc <- x+r*cos(rs) 
  yc <- y+r*sin(rs) 
  polygon(xc,yc,...) 
} 

left.half.circle <- function(x,y,r,nsteps=100,...){  
  rs <- seq(-pi/2,pi/2,len=nsteps) 
  xc <- x-r*cos(rs) 
  yc <- y-r*sin(rs) 
  polygon(xc,yc,...) 
} 

square <- function(x,y,r1,r2,...){  
  polygon(c(x-r1,x-r1,x+r1,x+r1),c(y-r2,y+r2,y+r2,y-r2),...) 
} 


# Requires pm and ps as vectors of length 3 each
plot_figure <- function(label, pms, pss, Ns, labels, output_path="", input_path="", loading_tail="") {

	
pdf(paste(output_path,"figure_",label,".pdf",sep=""), height=6, width=24)
par(mfrow=c(1,4), mar=c(5,6.5,1.5,1.5))

# Assort

load(paste(input_path,"pm_",pms[1],"_ps_",pss[1],"_N_",Ns,loading_tail,".RData",sep=""))

xs <- seq(-0.2,1,0.1)
xs.results <- matrix(NA,3,length(xs))
x.s <- max(xs)-min(xs)

for (i in 1:length(xs)) {
	xs.results[,i] <- quantile(assort.shrunk[assort.orig >= xs[i] - 0.05 & assort.orig <= xs[i] + 0.05], c(0.025,0.5,0.975), na.rm=TRUE)
}

plot(xs,xs, type='l', xlim=range(xs), ylim=range(xs.results,na.rm=T), xlab="Initial assortment", ylab="Resulting assortment", cex.lab=2.5, cex.axis=2)
for (i in 1:length(xs)) {
	lines(c(xs[i],xs[i])-0.01*x.s, c(xs.results[1,i], xs.results[3,i]),col="blue")
	lines(c(xs[i],xs[i]-0.02)-0.01*x.s, c(xs.results[1,i], xs.results[1,i]),col="blue")
	lines(c(xs[i],xs[i]-0.02)-0.01*x.s, c(xs.results[3,i], xs.results[3,i]),col="blue")
	left.half.circle(x=xs[i]-0.01*x.s,y=xs.results[2,i],r=0.03*x.s,col="blue",border=NA)
}

load(paste(input_path,"pm_",pms[2],"_ps_",pss[2],"_N_",Ns,loading_tail,".RData",sep=""))
xs <- seq(-0.2,1,0.1)
xs.results <- matrix(NA,3,length(xs))

for (i in 1:length(xs)) {
	xs.results[,i] <- quantile(assort.shrunk[assort.orig >= xs[i] - 0.05 & assort.orig <= xs[i] + 0.05], c(0.025,0.5,0.975), na.rm=TRUE)
	lines(c(xs[i],xs[i]), c(xs.results[1,i], xs.results[3,i]),col="orange")
	lines(c(xs[i],xs[i]), c(xs.results[1,i], xs.results[1,i]),col="orange")
	lines(c(xs[i],xs[i]), c(xs.results[3,i], xs.results[3,i]),col="orange")
	square(x=xs[i],y=xs.results[2,i],r1=0.009*x.s,r2=0.03*x.s,col="orange",border=NA)
}


load(paste(input_path,"pm_",pms[3],"_ps_",pss[3],"_N_",Ns,loading_tail,".RData",sep=""))
xs <- seq(-0.2,1,0.1)
xs.results <- matrix(NA,3,length(xs))

for (i in 1:length(xs)) {
	xs.results[,i] <- quantile(assort.shrunk[assort.orig >= xs[i] - 0.05 & assort.orig <= xs[i] + 0.05], c(0.025,0.5,0.975), na.rm=TRUE)
	lines(c(xs[i],xs[i])+0.01*x.s, c(xs.results[1,i], xs.results[3,i]),col="red")
	lines(c(xs[i],xs[i]+0.02)+0.01*x.s, c(xs.results[1,i], xs.results[1,i]),col="red")
	lines(c(xs[i],xs[i]+0.02)+0.01*x.s, c(xs.results[3,i], xs.results[3,i]),col="red")
	right.half.circle(x=xs[i]+0.01*x.s,y=xs.results[2,i],r=0.03*x.s,col="red",border=NA)
}

# Density

load(paste(input_path,"pm_",pms[1],"_ps_",pss[1],"_N_",Ns,loading_tail,".RData",sep=""))

xs <- seq(0,1,0.1)
xs.results <- matrix(NA,3,length(xs))
x.s <- max(xs)-min(xs)

for (i in 1:length(xs)) {
	xs.results[,i] <- quantile(den.shrunk[den.orig >= xs[i] - 0.05 & den.orig <= xs[i] + 0.05], c(0.025,0.5,0.975), na.rm=TRUE)
}

plot(xs,xs, type='l', xlim=range(xs), ylim=c(0,1), xlab="Initial density", ylab="Resulting density", cex.lab=2.5, cex.axis=2)
for (i in 1:length(xs)) {
	lines(c(xs[i],xs[i])-0.01*x.s, c(xs.results[1,i], xs.results[3,i]),col="blue")
	lines(c(xs[i],xs[i]-0.02)-0.01*x.s, c(xs.results[1,i], xs.results[1,i]),col="blue")
	lines(c(xs[i],xs[i]-0.02)-0.01*x.s, c(xs.results[3,i], xs.results[3,i]),col="blue")
	left.half.circle(x=xs[i]-0.01*x.s,y=xs.results[2,i],r=0.03,col="blue",border=NA)
}


load(paste(input_path,"pm_",pms[2],"_ps_",pss[2],"_N_",Ns,loading_tail,".RData",sep=""))
xs <- seq(0.1,1,0.1)
xs.results <- matrix(NA,3,length(xs))

for (i in 1:length(xs)) {
	xs.results[,i] <- quantile(den.shrunk[den.orig >= xs[i] - 0.05 & den.orig <= xs[i] + 0.05], c(0.025,0.5,0.975), na.rm=TRUE)
	lines(c(xs[i],xs[i]), c(xs.results[1,i], xs.results[3,i]),col="orange")
	lines(c(xs[i],xs[i]), c(xs.results[1,i], xs.results[1,i]),col="orange")
	lines(c(xs[i],xs[i]), c(xs.results[3,i], xs.results[3,i]),col="orange")
	square(x=xs[i],y=xs.results[2,i],r1=0.009*x.s,r2=0.03*x.s,col="orange",border=NA)
}


load(paste(input_path,"pm_",pms[3],"_ps_",pss[3],"_N_",Ns,loading_tail,".RData",sep=""))
xs <- seq(0,1,0.1)
xs.results <- matrix(NA,3,length(xs))

for (i in 1:length(xs)) {
	xs.results[,i] <- quantile(den.shrunk[den.orig >= xs[i] - 0.05 & den.orig <= xs[i] + 0.05], c(0.025,0.5,0.975), na.rm=TRUE)
	lines(c(xs[i],xs[i])+0.01*x.s, c(xs.results[1,i], xs.results[3,i]),col="red")
	lines(c(xs[i],xs[i]+0.02)+0.01*x.s, c(xs.results[1,i], xs.results[1,i]),col="red")
	lines(c(xs[i],xs[i]+0.02)+0.01*x.s, c(xs.results[3,i], xs.results[3,i]),col="red")
	right.half.circle(x=xs[i]+0.01*x.s,y=xs.results[2,i],r=0.03*x.s,col="red",border=NA)
}


# Clustering

load(paste(input_path,"pm_",pms[1],"_ps_",pss[1],"_N_",Ns,loading_tail,".RData",sep=""))

xs <- seq(0.2,1,0.1)
xs.results <- matrix(NA,3,length(xs))
x.s <- max(xs)-min(xs)

for (i in 1:length(xs)) {
	xs.results[,i] <- quantile(clust.shrunk[clust.orig >= xs[i] - 0.05 & clust.orig <= xs[i] + 0.05], c(0.025,0.5,0.975), na.rm=TRUE)
}

plot(xs,xs, type='l', xlim=range(xs), ylim=c(0,1), xlab="Initial clustering", ylab="Resulting clustering", cex.lab=2.5, cex.axis=2)
for (i in 1:length(xs)) {
	lines(c(xs[i],xs[i])-0.01*x.s, c(xs.results[1,i], xs.results[3,i]),col="blue")
	lines(c(xs[i],xs[i]-0.02)-0.01*x.s, c(xs.results[1,i], xs.results[1,i]),col="blue")
	lines(c(xs[i],xs[i]-0.02)-0.01*x.s, c(xs.results[3,i], xs.results[3,i]),col="blue")
	left.half.circle(x=xs[i]-0.01*x.s,y=xs.results[2,i],r=0.03*x.s,col="blue",border=NA)
}


load(paste(input_path,"pm_",pms[2],"_ps_",pss[2],"_N_",Ns,loading_tail,".RData",sep=""))
xs <- seq(0.2,1,0.1)
xs.results <- matrix(NA,3,length(xs))

for (i in 1:length(xs)) {
	xs.results[,i] <- quantile(clust.shrunk[clust.orig >= xs[i] - 0.05 & clust.orig <= xs[i] + 0.05], c(0.025,0.5,0.975), na.rm=TRUE)
	lines(c(xs[i],xs[i]), c(xs.results[1,i], xs.results[3,i]),col="orange")
	lines(c(xs[i],xs[i]), c(xs.results[1,i], xs.results[1,i]),col="orange")
	lines(c(xs[i],xs[i]), c(xs.results[3,i], xs.results[3,i]),col="orange")
	square(x=xs[i],y=xs.results[2,i],r1=0.009*x.s,r2=0.03*x.s,col="orange",border=NA)
}


load(paste(input_path,"pm_",pms[3],"_ps_",pss[3],"_N_",Ns,loading_tail,".RData",sep=""))
xs <- seq(0.2,1,0.1)
xs.results <- matrix(NA,3,length(xs))

for (i in 1:length(xs)) {
	xs.results[,i] <- quantile(clust.shrunk[clust.orig >= xs[i] - 0.05 & clust.orig <= xs[i] + 0.05], c(0.025,0.5,0.975), na.rm=TRUE)
	lines(c(xs[i],xs[i])+0.01*x.s, c(xs.results[1,i], xs.results[3,i]),col="red")
	lines(c(xs[i],xs[i]+0.02)+0.01*x.s, c(xs.results[1,i], xs.results[1,i]),col="red")
	lines(c(xs[i],xs[i]+0.02)+0.01*x.s, c(xs.results[3,i], xs.results[3,i]),col="red")
	right.half.circle(x=xs[i]+0.01*x.s,y=xs.results[2,i],r=0.03*x.s,col="red",border=NA)
}


# Degree


load(paste(input_path,"pm_",pms[1],"_ps_",pss[1],"_N_",Ns,loading_tail,".RData",sep=""))

xs <- seq(0,(Ns-2),4)
xs.results <- matrix(NA,3,length(xs))
x.s <- max(xs)-min(xs)

for (i in 1:length(xs)) {
	xs.results[,i] <- quantile(mean.deg.shrunk[mean.deg.orig >= xs[i] - 2 & mean.deg.orig <= xs[i] + 2], c(0.025,0.5,0.975), na.rm=TRUE)
}

plot(0:((Ns/2)-1),0:((Ns/2)-1), type='l', xlim=c(0,Ns), ylim=c(0,Ns-(Ns/5)), xlab="Initial mean degree", ylab="Resulting mean degree", cex.lab=2.5, cex.axis=2)
lines(((Ns/2)-1):(Ns-1),rep(((Ns/2)-1),((Ns/2)+1)))
for (i in 1:length(xs)) {
	lines(c(xs[i],xs[i])-0.01*x.s, c(xs.results[1,i], xs.results[3,i]),col="blue")
	lines(c(xs[i],xs[i]-0.02)-0.01*x.s, c(xs.results[1,i], xs.results[1,i]),col="blue")
	lines(c(xs[i],xs[i]-0.02)-0.01*x.s, c(xs.results[3,i], xs.results[3,i]),col="blue")
	left.half.circle(x=xs[i]-0.01*x.s,y=xs.results[2,i],r=0.03*x.s,col="blue",border=NA)
}


load(paste(input_path,"pm_",pms[2],"_ps_",pss[2],"_N_",Ns,loading_tail,".RData",sep=""))
xs <- seq(0,(Ns-2),4)
xs.results <- matrix(NA,3,length(xs))

for (i in 1:length(xs)) {
	xs.results[,i] <- quantile(mean.deg.shrunk[mean.deg.orig >= xs[i] - 2 & mean.deg.orig <= xs[i] + 2], c(0.025,0.5,0.975), na.rm=TRUE)
	lines(c(xs[i],xs[i]), c(xs.results[1,i], xs.results[3,i]),col="orange")
	lines(c(xs[i],xs[i]), c(xs.results[1,i], xs.results[1,i]),col="orange")
	lines(c(xs[i],xs[i]), c(xs.results[3,i], xs.results[3,i]),col="orange")
	square(x=xs[i],y=xs.results[2,i],r1=0.009*x.s,r2=0.03*x.s,col="orange",border=NA)
}


load(paste(input_path,"pm_",pms[3],"_ps_",pss[3],"_N_",Ns,loading_tail,".RData",sep=""))
xs <- seq(0,(Ns-2),4)
xs.results <- matrix(NA,3,length(xs))

for (i in 1:length(xs)) {
	xs.results[,i] <- quantile(mean.deg.shrunk[mean.deg.orig >= xs[i] - 2 & mean.deg.orig <= xs[i] + 2], c(0.025,0.5,0.975), na.rm=TRUE)
	lines(c(xs[i],xs[i])+0.01*x.s, c(xs.results[1,i], xs.results[3,i]),col="red")
	lines(c(xs[i],xs[i]+0.02)+0.01*x.s, c(xs.results[1,i], xs.results[1,i]),col="red")
	lines(c(xs[i],xs[i]+0.02)+0.01*x.s, c(xs.results[3,i], xs.results[3,i]),col="red")
	right.half.circle(x=xs[i]+0.01*x.s,y=xs.results[2,i],r=0.03*x.s,col="red",border=NA)
}

legend("topright", horiz=FALSE, legend=labels, pch=c(16,15,16,NA), lty=c(NA, NA, NA, 1), col=c("blue","orange","red","black"),cex=2,pt.cex=3.5)
mtext("(a)",side=3,outer=TRUE,at=0.008,line=-2.5,cex=2)
mtext("(b)",side=3,outer=TRUE,at=0.258,line=-2.5,cex=2)
mtext("(c)",side=3,outer=TRUE,at=0.508,line=-2.5,cex=2)
mtext("(d)",side=3,outer=TRUE,at=0.758,line=-2.5,cex=2)

dev.off()
# end Plot, end function

}

##############################

##############################

#### Section 4 - CODE TO PLOT THE FIGURES (requires RData files included in the supplemental materials) ####

# Figure 1
plot_figure("1",pms=c(0,0,0.2),pss=c(0,0.04,0),Ns=50, labels=c("No rewiring", "Random rewiring", "Second degree rewiring","No structural change"), input_path="~/Dropbox/Papers/2018_Network_shrinking/Revision2/Submitted/Data/", output_path="~/Dropbox/Papers/2018_Network_shrinking/Revision2/Submitted/")

# Figure S1
plot_figure("S1",pms=c(0,0,0.05),pss=c(0,0.01,0),Ns=50, labels=c("No rewiring", "Random rewiring", "Second degree rewiring","No structural change"), input_path="~/Dropbox/Papers/2018_Network_shrinking/Revision2/Submitted/Data/", output_path="~/Dropbox/Papers/2018_Network_shrinking/Revision2/Submitted/")

# Figure S2
plot_figure("S2",pms=c(0,0.05,0.2),pss=c(0,0.01,0.04),Ns=50, labels=c("No rewiring", "Both low rewiring", "Both high rewiring","No structural change"), input_path="~/Dropbox/Papers/2018_Network_shrinking/Revision2/Submitted/Data/", output_path="~/Dropbox/Papers/2018_Network_shrinking/Revision2/Submitted/")

# Figure S3
plot_figure("S3",pms=c(0,0.5,0.5),pss=c(0,0,0.01),Ns=50, labels=c("No rewiring", "Second degree rewiring", "Both rewiring","No structural change"), input_path="~/Dropbox/Papers/2018_Network_shrinking/Revision2/Submitted/Data/", output_path="~/Dropbox/Papers/2018_Network_shrinking/Revision2/Submitted/")

# Figure S4
plot_figure("S4",pms=c(0,0,0.2),pss=c(0,0.04,0),Ns=50, labels=c("No rewiring", "Random rewiring", "Second degree rewiring","No structural change"), input_path="~/Dropbox/Papers/2018_Network_shrinking/Revision2/Submitted/Data/", output_path="~/Dropbox/Papers/2018_Network_shrinking/Revision2/Submitted/", loading_tail=paste("_removed_",5,sep=""))

# Figure S5
plot_figure("S5",pms=c(0,0,0.2),pss=c(0,0.04,0),Ns=100, labels=c("No rewiring", "Random rewiring", "Second degree rewiring","No structural change"), input_path="~/Dropbox/Papers/2018_Network_shrinking/Revision2/Submitted/Data/", output_path="~/Dropbox/Papers/2018_Network_shrinking/Revision2/Submitted/")

# Figure S6
plot_figure("S6",pms=c(0,0,0.2),pss=c(0,0.04,0),Ns=24, labels=c("No rewiring", "Random rewiring", "Second degree rewiring","No structural change"), input_path="~/Dropbox/Papers/2018_Network_shrinking/Revision2/Submitted/Data/", output_path="~/Dropbox/Papers/2018_Network_shrinking/Revision2/Submitted/")

# Figure S7
plot_figure("S7",pms=c(0,0,0.1),pss=c(0,0.02,0),Ns=100, labels=c("No rewiring", "Random rewiring", "Second degree rewiring","No structural change"), input_path="~/Dropbox/Papers/2018_Network_shrinking/Revision2/Submitted/Data/", output_path="~/Dropbox/Papers/2018_Network_shrinking/Revision2/Submitted/")

# Figure S8
plot_figure("S8",pms=c(0,0,0.5),pss=c(0,0.1,0),Ns=24, labels=c("No rewiring", "Random rewiring", "Second degree rewiring","No structural change"), input_path="~/Dropbox/Papers/2018_Network_shrinking/Revision2/Submitted/Data/", output_path="~/Dropbox/Papers/2018_Network_shrinking/Revision2/Submitted/")


##############################

## END OF CODE

##############################
