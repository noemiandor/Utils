runPAGODA<-function(d_cell, NCORES=10,NGOTERMS=100,label=""){
  # Pathway and Gene Set Overdispersion Analysislibrary("scde")
  
  d_cell$mat=apply(d_cell$mat,2, as.integer)
  rownames(d_cell$mat)=d_cell$barcodes; colnames(d_cell$mat)=d_cell$genes
  print("Running knn error model...")
  knn <- knn.error.models(t(d_cell$mat), k = round(nrow(d_cell$mat)/4), n.cores = NCORES, min.count.threshold = 1, min.nonfailed = 15, max.model.plots = 10)
  
  ##Normalizing variance
  #In order to accurately quantify excess variance or overdispersion, we must normalize out expected levels of technical and intrinsic biological noise. Briefly, variance of the NB/Poisson mixture processes derived from the error modeling step are modeled as a chi-squared distribution using adjusted degrees of freedom and observation weights based on the drop-out probability of a given gene. Here, we normalize variance, trimming 3 most extreme cells and limiting maximum adjusted variance to 5.
  print("Removing techical/biological noise...")
  varinfo <- pagoda.varnorm(knn, counts = t(d_cell$mat), trim = 3/nrow(t(d_cell$mat)), max.adj.var = 5, n.cores = NCORES, plot = TRUE)
  
  #The plot on the left shows coefficient of variance squared (on log10 scale) as a function of expression magnitude (log10 FPM). The red line shows local regression model for the genome-wide average dependency. The plot on the right shows adjusted variance (derived based on chi-squared probability of observed/genomewide expected ratio for each gene, with degrees of freedom adjusted for each gene). The adjusted variance of 1 means that a given gene exhibits as much variance as expected for a gene of such population average expression magnitude. Genes with high adjusted variance are overdispersed within the measured population and most likely show subpopulation-specific expression:
  # list top overdispersed genes
  # sort(varinfo$arv, decreasing = TRUE)[1:10]
  
  ##Controlling for sequencing depth
  #Even with all the corrections, sequencing depth or gene coverage is typically still a major aspects of variability. In most studies, we would want to control for that as a technical artifact (exceptions are cell mixtures where subtypes significantly differ in the amount of total mRNA).
  #Below we will control for the gene coverage (estimated as a number of genes with non-zero magnitude per cell) and normalize out that aspect of cell heterogeneity:
  print("Controlling for sequencing depth...")
  varinfo <- pagoda.subtract.aspect(varinfo, rowSums(d_cell$mat[rownames(knn),]>0))
  
  ###Evaluate overdispersion of pre-defined gene sets
  #In order to detect significant aspects of heterogeneity across the population of single cells, 'pagoda' identifies pathways and gene sets that exhibit statistically significant excess of coordinated variability. Specifically, for each gene set, we tested whether the amount of variance explained by the first principal component significantly exceed the background expectation. We can test both pre-defined gene sets as well as 'de novo' gene sets whose expression profiles are well-correlated within the given dataset.
  #For pre-defined gene sets, we'll use GO annotations. 
  library(org.Hs.eg.db)
  # translate gene names to ids
  ids <- unlist(lapply(mget(d_cell$genes, org.Hs.egALIAS2EG, ifnotfound = NA), function(x) x[1]))
  rids <- names(ids); names(rids) <- ids
  # convert GO lists from ids to gene names
  gos.interest <- unique(c(ls(org.Hs.egGO2ALLEGS)[1:NGOTERMS]))
  go.env <- lapply(mget(gos.interest, org.Hs.egGO2ALLEGS), function(x) as.character(na.omit(rids[x])))
  go.env <- clean.gos(go.env) # remove GOs with too few or too many genes
  go.env <- list2env(go.env) # convert to an environment
  
  #Now, we can calculate weighted first principal component magnitudes for each GO gene set in the provided environment.
  print("Calculating weighted first principal component magnitudes for each GO gene set...")
  pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components = 1, n.cores = NCORES)
  
  #We can now evaluate the statistical significance of the observed overdispersion for each GO gene set.
  df <- pagoda.top.aspects(pwpca, return.table = TRUE, plot = TRUE, z.score = 1.96)
  
  #Each point on the plot shows the PC1 variance (lambda1) magnitude (normalized by set size) as a function of set size. The red lines show expected (solid) and 95% upper bound (dashed) magnitudes based on the Tracey-Widom model.
  # head(df)
  
  # The z column gives the Z-score of pathway over-dispersion relative to the genome-wide model (Z-score of 1.96 corresponds to P-value of 5%, etc.).
  # "z.adj" column shows the Z-score adjusted for multiple hypothesis (using Benjamini-Hochberg correction).
  # "score" gives observed/expected variance ratio
  # "sh.z" and "adj.sh.z" columns give the raw and adjusted Z-scores of "pathway cohesion", which compares the observed PC1 magnitude to the magnitudes obtained when the observations for each gene are randomized with respect to cells. When such Z-score is high (e.g. for GO:0008009) then multiple genes within the pathway contribute to the coordinated pattern.
  
  ### Evaluate overdispersion of 'de novo' gene sets
  # We can also test 'de novo' gene sets whose expression profiles are well-correlated within the given dataset. The following procedure will determine 'de novo' gene clusters in the data, and build a background model for the expectation of the gene cluster weighted principal component magnitudes. Note the higher trim values for the clusters, as we want to avoid clusters that are formed by outlier cells.
  print("Clustering gene sets...")
  clpca <- pagoda.gene.clusters(varinfo, trim = 7.1/ncol(varinfo$mat), n.clusters = 50, n.cores = NCORES, plot = TRUE)
  # The plot above shows background distribution of the first principal component (PC1) variance (lambda1) magnitude. The blue scatterplot on the left shows lambda1 magnitude vs. cluster size for clusters determined based on randomly-generated matrices of the same size. The black circles show top cluster in each simulation. The red lines show expected magnitude and 95% confidence interval based on Tracy-Widom distribution. The right plot shows extreme value distribution fit of residual cluster PC1 variance magnitude relative to the Gumbel (extreme value) distribution.
  # Now the set of top aspects can be recalculated taking these de novo gene clusters into account:
  df <- pagoda.top.aspects(pwpca, clpca, return.table = TRUE, plot = TRUE, z.score = 1.96)
  # head(df)
  
  
  ### Visualize significant aspects of heterogeneity
  # To view top heterogeneity aspects, we will first obtain information on all the significant aspects of transcriptional heterogeneity. We will also determine the overall cell clustering based on this full information:
  # get full info on the top aspects
  tam <- pagoda.top.aspects(pwpca, clpca, n.cells = NULL, z.score = qnorm(0.01/2, lower.tail = FALSE))
  # determine overall cell clustering
  print("Clustering cells...")
  hc <- pagoda.cluster.cells(tam, varinfo)
  # Next, we will reduce redundant aspects in two steps. First we will combine pathways that are driven by the same sets of genes:
  tamr <- pagoda.reduce.loading.redundancy(tam, pwpca, clpca)
  # In the second step we will combine aspects that show similar patterns (i.e. separate the same sets of cells). Here we will plot the cells using the overall cell clustering determined above:
  tamr2 <- pagoda.reduce.redundancy(tamr, distance.threshold = 0.9, plot = TRUE, cell.clustering = hc, labRow = NA, labCol = NA, box = TRUE, margins = c(0.5, 0.5), trim = 0)
  # In the plot above, the columns are cells, rows are different significant aspects, clustered by their similarity pattern.The green-to-orange color scheme shows low-to-high weighted PCA scores (aspect patterns), where generally orange indicates higher expression. Blocks of color on the left margin show which aspects have been combined by the command above. Here the number of resulting aspects is relatively small. "top" argument (i.e. top = 10) can be used to limit further analysis to top N aspects.
  # We will view the top aspects, clustering them by pattern similarity (note, to view aspects in the order of increasing lambda1 magnitude, use row.clustering = NA).
  col.cols <- rbind(groups = cutree(hc, 3))
  
  outTiff=paste(label,"_PAGODEaspectsView.tiff",sep="")
  tiff(filename = outTiff, width=5.55, height=4.65, units="in", res=200)
  pagoda.view.aspects(tamr2, cell.clustering = hc, box = TRUE, labCol = NA, margins = c(0.5, 20), col.cols = substr(hc$labels,nchar(hc$labels[1]),nchar(hc$labels[1])))
  # While each row here represents a cluster of pathways, the row names are assigned to be the top overdispersed aspect in each cluster.
  dev.off()
  print(paste("Displayed clustering aspects under: ",outTiff))
  return(list(tamr=tamr,tamr2=tamr2,hc=hc,clpca=clpca,knn=knn))
}