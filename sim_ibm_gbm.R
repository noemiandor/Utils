sim_ibm_gbm <- function(ibm.env, to = 100, watch = T, outPrefix = NULL, gridRows = 50, gridCols = 38) {
  library(simecol)
  library(matlab)
  library(fields)
  
  # Neighborhood weight matrix for energy diffusion:
  assign("wdist", matrix(c(0.25,0.25,0.25,0.25,0.25,0.25,0.25,
                           0.25,0.50,0.50,0.50,0.50,0.50,0.25,
                           0.25,0.50,1.00,1.00,1.00,0.50,0.25,
                           0.25,0.50,1.00,1.00,1.00,0.50,0.25,
                           0.25,0.50,1.00,1.00,1.00,0.50,0.25,
                           0.25,0.50,0.50,0.50,0.50,0.50,0.25,
                           0.25,0.25,0.25,0.25,0.25,0.25,0.25), nrow=7), env=ibm.env); 
  ### Simecol model:
  ibm_gbm = getIBM_GBM(ibm.env, to = to, watch = watch, gridRows = gridRows, gridCols = gridCols)
  watch && !is.null(outPrefix) && pdf(paste0(outPrefix, ".pdf"))
  times(ibm_gbm) = c(to = to)
  ibm_gbm@parms$watch = watch
  
  ### Runs the simulation
  X_10 <- sim(ibm_gbm)
  watch && !is.null(outPrefix) && dev.off()
  
  return(X_10)
}

getIBM_GBM <- function(ibm.env, to = 100, watch = T, gridRows = 50, gridCols = 38) {
  ibm_gbm <- new("gridModel", main = function(time, init, parms) {
    ## Format of entries in init: p|d|G , where p = proliferation rate, d = death rate and G is a vector with 22 entries, each encoding the ploidy of an autosome Example: '0.11|0.05|3233333333333343333333'
    if (time == times[2]) {
      print("initializing")
      pr = local(growthRate, env = ibm.env)
      dr = local(deathRate, env = ibm.env)
      ## Start out with resection margins:
      # iR=nrow(init)/6
      # iC=ncol(init)/6
      # coord=cbind(3*iR + iR*sin(seq(0,2*pi,length.out=100)),3*iC + iC*cos(seq(0,2*pi,length.out=100)))
      # for(i in 1:nrow(coord)){
      #   init[round(coord[i,1]), round(coord[i,2])] = paste0(pr,"|",dr,"|",paste(parms$B_genome,collapse = ""))
      # }
      ## Start out with core:
      iR = round((nrow(init) * 4/9):(nrow(init) * 5/9))
      iC = round(ncol(init) * 4/9):(ncol(init) * 5/9)
      init[sample(iR, 5), sample(iC, 5)] = paste0(pr, "|", dr, "|", paste(parms$B_genome, collapse = ""))
    }
    init = survive(init)
    if (mod(time, 20) == 0) {
      print(paste("Go or grow @", time))
    }
    init = go_or_grow(init, parms)
  },
  parms = list(
    ibm.env = environment()$ibm.env, 
    B_genome = rep(round(local(ploidy, env=ibm.env)), 22), 
    maxNeighbors = eightneighbors(matrix(1, nrow = gridRows, ncol = gridCols)), 
    gir  = local(gir, env=ibm.env), 
    wdist = local(wdist, env=ibm.env),
    maxEnergy = local(Moore_E, env=ibm.env) * 
      neighbours(matrix(1, nrow = gridRows, ncol = gridCols), 
                 wdist = local(wdist, env=ibm.env), bounds = 0),
    watch = T
  ),
  init   = matrix(NA, gridRows, gridCols), 
  times  = c(from=0, to=to, by=1),
  solver = "iteration",
  equations   = list()
  )
  
  equations(ibm_gbm) <- list(
    genomeupdate = function(genome, parms) {
    ## If region does not contain CNV already, GI increase is certainty, otherwise, non-negligeable chance of GI decrease
    with(parms, {
      genomes <- genome
      mgenome <- apply(genomes, 1, function(genome) {
        x = rep(0, 1/parms$gir)
        x[1:length(genome)] = 1:length(genome)
        ## Cell can't gain deletion back:
        x[genome == 0] = 0
        ## Maximum of one CNV per division:
        ii = sample(x, 1)
        genome[ii] <- genome[ii] + sample(c(-1, 1), 1)
        ## Only one digit per chr:
        genome[genome > 9] = 9
        genome[genome < 0] = 0
        return(genome)
      })
      mgenome = t(mgenome)
    })
  }, isalive = function(inds) {
    state = apply(!is.na(inds), 2, as.numeric)
    return(state)
  }, survive = function(inds) {
    ## The survival function returns the subset of surviving individuals.
    state = isalive(inds)
    pdeath = getDeathRate(inds)
    zdie = ifelse(state >= 1 & runif(state) < pdeath, TRUE, FALSE)
    inds[which(zdie)] = NA
    return(inds)
  }, getFitness = function(inds) {
    return(getProliferationRate(inds) - getDeathRate(inds))
  }, getRemainingEnergy = function(inds, parms) {
    with(parms, {
      ploidymat <- getPloidy(inds)
      ploidymat[is.na(ploidymat)] = 0
      lockedEnergy <- neighbours(ploidymat, wdist = wdist, bounds = 0)
      return(maxEnergy - lockedEnergy)
    })
  }, go_or_grow = function(inds, parms) {
    ## Division is initiated in a space and fitness dependent manner and the resulting new genome and fitness are calculated.  Population growth occurs by setting the corresponding grid addresses.
    with(parms, {
      state <- isalive(inds)
      nb <- eightneighbors(state)
      prolifrate <- getProliferationRate(inds)
      ploidymat <- getPloidy(inds)
      ploidymat[is.na(ploidymat)] = 0
      
      lockedEnergy <- neighbours(ploidymat, wdist = wdist, bounds = 0)
      requiredEnergy <- ploidymat * 0.5 * nrow(wdist) * ncol(wdist)
      availableEnergy <- maxEnergy - lockedEnergy - requiredEnergy
      ## Decision btw. Go and Grow: Cells that are allowed to divide have dcan>0
      dcan <- prolifrate * sign(maxNeighbors - nb)
      ## Dividing/migrating cells #par(mfrow=c(2,1)); xgen=zgen; image(xgen);
      zgen <- ifelse(state == 1 & runif(state) < dcan, 1, 0)
      ## Migrating cells: if they can't grow, they go
      mgen <- zgen * as.matrix(availableEnergy <= 0)
      ## Dividing cells
      zgen <- zgen * as.matrix(availableEnergy > 0)
      
      ## Logistics: Dividing cells' index
      have.neo <- which(zgen == 1, arr.ind = T)
      ## Candidate locations for new cells
      are.neo <- which(sign(eightneighbors(zgen)) - state > 0, arr.ind = T)
      ## Migrating cells' index
      do.mov <- which(mgen == 1, arr.ind = T)
      ## Candidates for migrating cells' new locations
      are.mov <- which(sign(eightneighbors(mgen)) - state > 0, arr.ind = T)
      ## Make sure candidate locations between go vs. grow cells don't overlap:
      o = intersect_Matlab(are.mov, are.neo)
      if (!isempty(o$ia)) {
        are.neo = are.neo[-sample(o$ib, round(length(o)/2)), , drop = F]
        are.mov = are.mov[-intersect_Matlab(are.mov, are.neo)$ia, , drop = F]
      }
      ## Locations of daughter cells:
      matchNewLocations <- function(have.neo, are.neo) {
        if (isempty(are.neo)) {
          return(are.neo)
        }
        d <- flexclust::dist2(have.neo, are.neo)
        idx <- lapply(1:nrow(d), function(i) which(d[i, , drop = F] < 2))
        names(idx) <- rownames(have.neo) <- as.character(1:length(idx))
        ## reduce dividing/migrating cells count if space is insufficient
        idx <- idx[sapply(idx, length) > 0]
        have_are.neo <- lapply(names(idx), function(i) rbind(have.neo[i, , drop = F], are.neo[idx[[i]], , drop = F]))
        return(have_are.neo)
      }
      have_are.neo = matchNewLocations(have.neo, are.neo)
      do_are.mov = matchNewLocations(do.mov, are.mov)
      
      ###############################
      ###Deal with "grow"-scenario###
      if (!isempty(have_are.neo) && !is.na(have_are.neo)) {
        ## Locations of dividing and daughter cells:
        have.neo <- do.call(rbind, lapply(have_are.neo, function(x) x[1, ]))
        ## Randomly sample among location candidates (except 1st because that is the dividing cell)
        are.neo <- lapply(have_are.neo, function(x) x[1 + sample(nrow(x) - 1, 1), ])
        are.neo <- do.call(rbind, are.neo)
        watch && heatmap(getPloidy(inds), symm = F, scale = "none", Rowv = NA, Colv = NA)[[1]]
        ## Update genome
        oldgenomes <- getGenome(inds[have.neo], as.numeric = T)
        genomes1 <- genomeupdate(oldgenomes, parms)
        genomes2 <- 2 * oldgenomes - genomes1
        ploidy <- apply(genomes1, 1, .calcPloidy)
        ismut <- which(ploidymat[have.neo] - ploidy != 0)
        ## Fitness: update growth- AND death rate
        newdeathrates <- repmat(getDeathRate(inds[have.neo]), 1, 2)
        newprolifrates <- repmat(getProliferationRate(inds[have.neo]), 1, 2)
        for (no in 1:2) {
          newdeathrates[ismut, no] <- newdeathrates[ismut, no] * selectionEffect_Death(ploidymat[have.neo][ismut])
          newprolifrates[ismut, no] <- newprolifrates[ismut, no] * selectionEffect_Proliferation(length(ismut))
        }
        ## Create dividing cells (assymmetric divisions):
        inds[are.neo] <- paste0(newprolifrates[, 1], "|", newdeathrates[, 1], "|", apply(genomes1, 1, paste, collapse = ""))
        inds[have.neo] <- paste0(newprolifrates[, 2], "|", newdeathrates[, 2], "|", apply(genomes2, 1, paste, collapse = ""))
      }
      
      #############################
      ###Deal with "go"-scenario###
      if (!isempty(do_are.mov) && !is.na(do_are.mov)) {
        # Calculate smooth gaussian of cell density
        psm = mat2XY(image.smooth(ploidymat, theta = 2)$z)
        # Get positions of minimum density: X
        ii_min = psm[psm$v <= quantile(psm$v, 0.1), c("r", "c")]
        # Among all location candidates (including current location) choose the one that brings cell closest to any X Distance to low-density coordinates
        dx = lapply(do_are.mov, function(x) flexclust::dist2(x, ii_min))
        # Choose closest candidate
        dx = lapply(dx, function(x) which(x == min(x), arr.ind = T)[, "row"])
        # But make sure it's jsut one
        to = lapply(dx, function(x) x[sample(length(x), 1)])
        # Which is the one
        to = lapply(1:length(to), function(i) do_are.mov[[i]][to[[i]], ])
        to = do.call(rbind, to)
        from = do.call(rbind, lapply(do_are.mov, function(x) x[1, ]))
        # dx<-dx2<-image.smooth(state, theta=2)$z; image(dx); dx[t(to)]=0.5; dx2[t(from)]=0.6; image(dx2+dx) Do the actual moving:
        inds[to] <- inds[from]
        inds[from] <- NA
      }
      return(inds)
    })
  }, selectionEffect_Death = function(ploidies) {
    return(rnorm(length(ploidies), mean = 1 * local(ploidy, env = ibm.env)/ploidies, sd = 0.15))
  }, selectionEffect_Proliferation = function(n) {
    return(rnorm(n, mean = 0.95, sd = 0.15))
  })
  return(ibm_gbm)
}

######################
###Helper functions###
getProliferationRate = function(inds) {
  f = sapply(strsplit(inds, "|", fixed = T), function(x) x[1])
  f[is.na(f)] = 0
  return(as.numeric(f))
}

getDeathRate = function(inds) {
  f = sapply(strsplit(inds, "|", fixed = T), function(x) x[min(2, length(x))])
  f[is.na(f)] = 1
  return(as.numeric(f))
}

getGenome = function(inds, as.numeric = F) {
  genomes = strsplit(inds, "|", fixed = T)
  genomes = sapply(genomes, function(x) x[[length(x)]])
  if (as.numeric) {
    genomes = t(sapply(strsplit(genomes, ""), as.numeric))
  }
  return(genomes)
}

getPloidy = function(inds) {
  genomes = getGenome(inds, as.numeric = T)
  inds[T] = sapply(genomes, .calcPloidy)
  return(apply(inds, 2, as.numeric))
}

.calcPloidy <- function(x) {
  return(sum(x * (length(x):1))/sum(1:length(x)))
  # return(mean(x))
}

mat2XY <- function(state) {
  la = data.frame(r = as.vector(row(state)), c = as.vector(col(state)), v = as.vector(state))
  return(la)
}

intersect_Matlab <- function(A, B) {
  A = apply(A, 1, paste, collapse = "_")
  B = apply(B, 1, paste, collapse = "_")
  A = gsub(" ", "", A)
  B = gsub(" ", "", B)
  x = intersect(A, B)
  ia = match(x, A)
  ib = match(x, B)
  return(list(ia = ia, ib = ib))
}
