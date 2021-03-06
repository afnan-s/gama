# libraries are used along the code by explicit refer syntax package::function()
source("../gama_sources/gama/R/fitness.R")
source("../gama_sources/gama/R/bestk.R")

gama.env <- new.env(parent = emptyenv())

# generates an initial population of candidate centers for partitions
# using two possible methods: pop.f and pop.f.2

# generate randomly from the range lower:upper
pop.f <- function(object) {

  lower <- object@lower
  upper <- object@upper
  nvars <- length(lower)

  population <- matrix(as.double(NA), nrow = object@popSize, ncol = gama.env$dims * gama.env$k)

  # generate the initial population based on uniform distribution
  for(j in 1:nvars) {
    population[,j] <- runif(object@popSize, lower[j], upper[j])
  }

  return(population)
}

# generate randomly from only pre-observed values for each variable:
pop.f.2 <- function(object){
  lower <- object@lower
  upper <- object@upper
  nvars <- length(lower)

  population <- matrix(as.double(NA), nrow = object@popSize, ncol = gama.env$dims * gama.env$k)


  for(j in 1:nvars) {
    sampled <- gama.env$dataset[sample(1:gama.env$rows,size=object@popSize, replace=FALSE),((j-1)%%gama.env$dims)+1]
    population[,j] <- as.matrix(sampled)

  }
  return(population)
}

#penalty function, penalise based on shortage of number of clusters
enforce.k.penalty <- function(m.individual) {

  penalty <- 0
  # counts how many centroids successfully assigned members:
  which.dists <- apply(Rfast::dista(gama.env$dataset, m.individual, "euclidean", square = TRUE), 1, which.min)
  num.clusters <- length(unique(which.dists))
  difference <- abs(num.clusters - gama.env$k )
  #difference <- num.clusters
  #if there's a difference, the penalty is the proportion:
  if (difference > 0) {
    penalty <- difference / gama.env$k
    #penalty <- difference
    #penalty <- num.clusters
  }
  #cat("In penalty: ", penalty, "\n")
  return (penalty)
}

# entry point of the gama package
gama <- function(dataset = NULL, d2 = NULL, k = "broad", scale = FALSE, crossover.rate = 0.9,
                         mutation.rate = 0.01, elitism = 0.05, pop.size = 25,
                         generations = 100, seed.p = 42,
                         fitness.criterion = "ASW",
                         penalty.function = enforce.k.penalty,
                         plot.internals = TRUE, ...) {

  # --- arguments validation --- #
  Check <- ArgumentCheck::newArgCheck()

  if (is.null(dataset))
    ArgumentCheck::addError(
      msg = "'dataset' can not be NULL",
      argcheck = Check)

  if (class(dataset) != 'data.frame')
    ArgumentCheck::addError(
      msg = "'dataset' must be a data.frame object.",
      argcheck = Check)
  if(class(d2) != 'dist')
    ArgumentCheck::addError(
      msg = "'d2' must be a dist object (see 'dist' function in 'cluster package). Or send null and Distance Matrix will be calculated automatically. ",
      argcheck = Check)
  if (is.null(k))
    ArgumentCheck::addError(
      msg = "'k' can not be NULL",
      argcheck = Check)

  if (is.numeric(k)) {
    # forces k to be an integer value
    k <- as.integer(k)

    if (k < 2) {
      ArgumentCheck::addError(
        msg = "'k' must be a positive integer value greater than one (k > 1), or one of the methods to estimate it: 'minimal' or 'broad'.",
        argcheck = Check)
    }
  } else if (is.character(k)) {
    if (!(k %in% c('minimal', 'broad')))
      ArgumentCheck::addError(
        msg = "'k' must be a positive integer value greater than one (k > 1), or one of the methods to estimate it: 'minimal' or 'broad'.",
        argcheck = Check)
  }

  if (is.null(fitness.criterion)) {
    ArgumentCheck::addError(
      msg = "'fitness.criterion' can not be NULL.",
      argcheck = Check)

  } else if (!(fitness.criterion %in% c('ASW', 'CH', 'CI', 'DI'))) {
    ArgumentCheck::addError(
      msg = "'fitness.criterion' must be one of the criteria: 'ASW', 'CH', 'CI' or 'DI'.",
      argcheck = Check)
  }

  if (!is.null(penalty.function)) {
    if (class(penalty.function) != "function") {
      ArgumentCheck::addError(
        msg = "'penalty.function' must be a valid R function.",
        argcheck = Check)
    }
  }

  ArgumentCheck::finishArgCheck(Check)

  # --- final of arguments validation --- #

  call <- match.call()

  # choose the correct fitness function
  fitness.function <- switch (fitness.criterion,
                             "ASW" = fitness.asw,
                             "CH" = fitness.ch,
                             "CI" = fitness.ci,
                             "DI" = fitness.di,
                             fitness.asw)

  if (scale) {
    dataset <- scale(dataset, center = TRUE, scale = TRUE)
  }

  # if k exists
  if (!is.na(k)) {
    #if k is string (a method to choose the estimative)
    if (is.character(k)) {
      estimative.method <- k
      k <- gama.how.many.k(dataset = dataset, method = estimative.method)
    }
  }

  dims <- ncol(dataset)
  rows <- nrow(dataset)
  gama.env$dataset <- dataset
  gama.env$k = k
  gama.env$dims = dims
  gama.env$rows = rows

  #elit.rate = floor(pop.size * elitism)
  elit.rate = base::max(1, round(pop.size*elitism))

  if(is.null(d2)){
    d <- dist(dataset, method = "euclidean", diag = FALSE, upper = FALSE)
    d2 <- d^2
    rm(d); gc()
  }

  gama.env$d2 <- d2


  lowers <- apply(dataset, 2, min)
  uppers <- apply(dataset, 2, max)

  lower_bound <- unlist(lapply(lowers, function (x) { rep(x, k) } ))
  upper_bound <- unlist(lapply(uppers, function (x) { rep(x, k) } ))

  set.seed(seed.p)

  # as defined by experimental procedure
  s <- "gareal_lrSelection"
  c <- "gareal_blxCrossover"
  m <- "gareal_nraMutation"

  start.time <- Sys.time()
  gama.env$start.time <- start.time
  gama.env$checked.time <- start.time

  # call GA functions
  genetic <- GA::ga(type = "real-valued",
                    seed = seed.p,
                    population = pop.f.2,
                    selection = s,
                    mutation = m,
                    crossover = c,
                    popSize = pop.size,
                    elitism = elit.rate,
                    pmutation = mutation.rate,
                    pcrossover = crossover.rate,
                    maxiter = generations,
                    fitness = fitness.function,  
                    penalty.function = penalty.function,
                    lower = lower_bound,
                    upper = upper_bound,
                    parallel = TRUE,
                    monitor = TRUE)


  end.time <- Sys.time()
  rt <- end.time - start.time

  num_solutions = length(genetic@solution)/(k*dims)
  if (num_solutions == 1) {
    solution <- matrix(genetic@solution, nrow = k, ncol = dims)
  } else {
    # if there is more than a single solution (they are identical, in ASW,
    # and must be close for centroids values)
    solution <- matrix(genetic@solution[1,], nrow = k, ncol = dims)
  }

  # calculates the average silhouette width
  which.dists <- apply(Rfast::dista(dataset, solution, "euclidean", square = TRUE), 1, which.min)

  asw <- cluster::silhouette(which.dists, d2)
  ch <- clusterCrit::intCriteria(as.matrix(dataset), which.dists, "Calinski_Harabasz")
  ci <- clusterCrit::intCriteria(as.matrix(dataset), which.dists, "C_index")
  di <- clusterCrit::intCriteria(as.matrix(dataset), which.dists, "Dunn")
  print(class(asw))
  # builds the solution object
  solution.df <- as.data.frame(solution)
  colnames(solution.df) <- colnames(dataset)
  solution.df <- solution.df[with(solution.df, order(apply(solution.df, 1, sum))), ]

  object <- methods::new("gama",
                original.data = as.data.frame(dataset),
                centers = solution.df,
                cluster = as.vector(which.dists),
                silhouette = summary(asw),
                calinski_harabasz = ch$calinski_harabasz,
                c_index = ci$c_index,
                dunn_index = di$dunn,
                runtime = round(rt, 2),
                call = call)

  print(object)

  # plot the results
  if (plot.internals) {
    plot(genetic, main = "Evolution")

    lim <- 500

    if (nrow(gama.env$dataset) > lim) {
      cat("\nIMPORTANT!!!\nThe dataset contains ", nrow(gama.env$dataset), "rows. To improve the quality of the ASW graph, gama will generate a file 'asw_gama.pdf' in working directory.\n")
      grDevices::pdf(file = 'asw_gama.pdf')
    }

    plot(asw, main = "Average Silhouette Width")

    if (nrow(gama.env$dataset) > lim) {
      garbage <- grDevices::dev.off()
    }
  }

  # return an object of class 'gama'
  return (object)

}

setClass(Class = "gama",
         slots = c(original.data = "data.frame",
         centers = "data.frame",
         cluster = "vector",
         silhouette = "ANY",
         calinski_harabasz = "numeric",
         c_index = "numeric",
         dunn_index = "numeric",
         runtime = "ANY",
         call = "ANY"))

print.gama <- function(x, ...) {

  # cat("\nDetails for the object of class 'gama':\n")

  # cat("\nOriginal data (first rows):\n")
  # print(head(x@original.data))
  # cat("\nCluster partitions:\n")
  # print(x@cluster)
  # cat("\nCluster Centers:\n")
  # print(x@centers)
  # cat("\nAverage Silhouette Width index (ASW):\n")
  # print(x@silhouette$avg.width)
  # cat("\nCalinski Harabasz index (CH):\n")
  # print(x@calinski_harabasz)
  # cat("\nC-Index (CI):\n")
  # print(x@c_index)
  # cat("\nDunn index (DI):\n")
  # print(x@dunn_index)

  # cat("\nCall:\n")
  # print(x@call)
  # cat("\nRuntime:\n")
  # print(x@runtime)
  cat("\nSilhouette Scores:\n")
  print(x@silhouette)
  cat("\nNumber of clusters:\n")
  print(length(unique(x@cluster)))

}
setMethod("print", "gama", print.gama )
