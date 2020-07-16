#Code modified from GAMA:

# This function serves to compare results of GA clustering with a set of randomly generated solutions. 
# Uses the same method of random generations using pop.f.2
gama.random <- function(dataset, d2=NULL, popSize, k) {

  #Calculate Lower and upper bounds vector:
  lowers <- apply(dataset, 2, min)
  uppers <- apply(dataset, 2, max)

  lower <- unlist(lapply(lowers, function (x) { rep(x, k) } ))
  upper <- unlist(lapply(uppers, function (x) { rep(x, k) } ))

  #number of variables:
  nvars <- length(lower)
  cat("nvars: ",nvars, "\n")
  cat("ncol(dataset)* k = ",ncol(dataset)*k, "\n")
  #number of columns and rows in dataset (cols, rows):
  cols <- ncol(dataset)
  rows <- nrow(dataset)

  population <- matrix(as.double(NA), nrow = popSize, ncol = cols * k)
  cat("Generating initial population..")

  # generate the initial population based on uniform distribution
  for(j in 1:nvars) {
    #old method:
    #population[,j] <- runif(popSize, lower[j], upper[j])
    
    #Generate only pre-observed values for each variable (uniformly):
    sampled <- dataset[sample(1:rows,size=popSize, replace=TRUE),((j-1)%%cols)+1]
    population[,j] <- as.matrix(sampled)

  }
  
  cat("Done.\n")
  #return(population)
  #If not given, calculate the distance:
  if(is.null(d2)){
  	cat("Now calculating distance matrix..")
  	d <- dist(dataset, method = "euclidean", diag = FALSE, upper = FALSE)
  	d2 <- d^2
  	cat("Done.\n")
  }


  sol_evals <- data.frame()
  #Now Evaluate the solutions:
  cat("Now evaluating each of the solutions..\n")
  for(i in 1:popSize){
  	#cat("Solution",i)
  	solution <- matrix(population[i,], nrow = k, ncol = cols)
  	which.dists <- apply(Rfast::dista(dataset, solution, "euclidean", square = TRUE), 1, which.min)
	num.clusters <- length(unique(which.dists))
	#cat(". Number of clusters: ", num.clusters)
	if(num.clusters > 1){
		asw <- cluster::silhouette(which.dists, d2)
		#print(summary(asw))
		#print(summary(asw)$avg.width)
		sol_evals <- rbind(sol_evals, cbind(i, summary(asw)$avg.width, num.clusters))
		#cat(". ASW: ", summary(asw)$avg.width)
	}
	#cat("\n")

  }
  names(sol_evals) <- c("number", "asw", "clusters")
  cat("Completed evaluating ",i, " solutions.\n" )
  print(summary(sol_evals$asw))
  best <- sol_evals[which.max(sol_evals$asw),]
  cat("Best solutions: \n")
  print(best)
  s.d <- sd(sol_evals$asw)
  cat("\nStandard Deviation: ")
  print(s.d)

  #print(sol_evals)
  return (sol_evals)
}

naive.random <- function(dataset, d2=NULL, popSize, k){
  if(is.null(d2)){
  	cat("Now calculating distance matrix..")
  	d <- dist(dataset, method = "euclidean", diag = FALSE, upper = FALSE)
  	d2 <- d^2
  	cat("Done.\n")
  }
  sol_evals <- data.frame()
  #Now Evaluate the solutions:
  for(i in 1:popSize){
  	#cat("Solution",i)
  	which.dists <- sample(1:k,size=nrow(dataset), replace=TRUE)
	num.clusters <- length(unique(which.dists))
	#cat(". Number of clusters: ", num.clusters)
	if(num.clusters > 1){
		asw <- cluster::silhouette(which.dists, d2)
		#print(summary(asw))
		#print(summary(asw)$avg.width)
		sol_evals <- rbind(sol_evals, cbind(i, summary(asw)$avg.width, num.clusters))
		#cat(". ASW: ", summary(asw)$avg.width)
	}
	#cat("\n")

  }
  names(sol_evals) <- c("number", "asw", "clusters")
  cat("Completed evaluating ",i, " solutions.\n" )
  print(summary(sol_evals$asw))
  cat("Standard Deviation = ", sd(sol_evals$asw), "\n")
  max.sol <- sol_evals[which.max(sol_evals$asw),]
  cat("Best Solutions: " )
  print(max.sol)
  cat("\n")
  #print(sol_evals)
  return (sol_evals)
}

#A function that generates random populations and check the k that produces the 
#best solution, it returns a data frame of trend of k + generation producing the
#best k. 
find.best.k <- function(dataset, d2=NULL, popSize = 2000){
	if(is.null(d2)){
  	cat("Now calculating distance matrix..")
  	d <- dist(dataset, method = "euclidean", diag = FALSE, upper = FALSE)
  	d2 <- d^2
  	cat("Done.\n")
  }

  dataset.size <- nrow(dataset)
  k_evals <- data.frame()
  #Now Evaluate the solutions:
  max.asw = -1
  max.clusters = NULL
  best.pop = NULL
  #for(k in seq(2, dataset.size/2, by=250)){
  for(k in 2:6){
  	#cat("K: ", k, " | ")
  	population <- pop.f.2(dataset = dataset, popSize = popSize, k = k)

  	# Loop each of the solutions, find their best and mean ASW. 
  	sol_evals <- data.frame()
  	
  	for(i in 1:popSize){
  	#cat("Solution",i)
  	solution <- matrix(population[i,], nrow = k, ncol = ncol(dataset))
  	which.dists <- apply(Rfast::dista(dataset, solution, "euclidean", square = TRUE), 1, which.min)
	num.clusters <- length(unique(which.dists))
	if(num.clusters > 1){
		asw <- cluster::silhouette(which.dists, d2)
		sol_evals <- rbind(sol_evals, cbind(summary(asw)$avg.width, num.clusters))
	}
	#cat("\n")

  } # End solutions loop
  names(sol_evals) <- c("asw", "clusters")
  sd <- sd(sol_evals$asw)
  current.mean <- mean(sol_evals$asw)
  current.mean.clusters <- mean(sol_evals$clusters)
  current.max <- sol_evals[which.max(sol_evals$asw),]
  k_evals <- rbind(k_evals, cbind(k, current.mean, sd, current.mean.clusters, current.max[,1],current.max[,2] ))
  print(tail(k_evals, 1))
  cat("\n")

  if(current.max[,1] > max.asw){
  	max.asw = current.max[,1]
  	max.clusters = current.max[,2]
  	best.pop = population
  }

  rm(sol_evals)
  rm(population)
  rm(which.dists)
  } # End k loop
  names(k_evals) <- c("k", "mean.asw", "sd", "mean clusters", "max.asw", "max.asw.clusters")

  library(ggplot2)

  theme_set(theme_minimal())
  ggplot(k_evals, aes(x=k)) + 
  geom_line(aes(y = mean.asw, color = "Mean Silhouette")) + 
  geom_line(aes(y = max.asw, color="Max Silhouette")) +
  scale_colour_manual("", 
                      values = c("Mean Silhouette"="green", "Max Silhouette"="red")) +
  labs(x="K (Number of Clusters)", y = "Silhouette Width") +
  ggtitle(paste("Achieved Mean and Max Sihouette Width\nFor Population = ",popSize))
  ggsave(paste("find_best_k_pop_", popSize,".pdf"))

  return(list("asw" = max.asw, "clusters"= max.clusters, "population" = best.pop, "k_table" = k_evals))

} # End find.best.k function. 

pop.f.2 <- function(dataset = dataset, popSize = popSize, k = k){
  
  cols <- ncol(dataset)
  rows <- nrow(dataset)
  nvars <- cols * k
  population <- matrix(as.double(NA), nrow = popSize, ncol = cols * k)


  for(j in 1:nvars) {
    sampled <- dataset[sample(1:rows,size=popSize, replace=FALSE),((j-1)%%cols)+1]
    population[,j] <- as.matrix(sampled)

  }
  return(population)
}

