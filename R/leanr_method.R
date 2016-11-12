## FUNCTIONS

# register the parallel backend: 
# Tries using library doMC (available on linux/MacOS)
# if this does not work parallel execution is disabled with a warning
register.backend<-function(cores=NULL){
  if (requireNamespace("doMC", quietly = TRUE)) {
    if (is.null(cores)) doMC::registerDoMC()else doMC::registerDoMC(cores=cores)
  } else {
      warning('Library doMC could not be found. Parallel execution disabled.')
  }
}

# load network from sif file
sif2graph<-function(siffile){
  tmpfile<-sprintf('%s/%s.ncol',tempdir(),basename(siffile))
	dat<-read.table(siffile,sep='\t',stringsAsFactors=F)
	scores_scaled<- rep(1,dim(dat)[1])
	write.table(cbind(dat[,c(1,3)],scores_scaled),tmpfile,sep='\t',quote=F,col.names=F,row.names=F)
	g <- read.graph(tmpfile,format='ncol')
	file.remove(tmpfile)
	g
}

subnet.simulation<-function(g,nmods=10,mod_lims=c(10,50),pval_scaling=0.1,mod_enrich_perc=0.5,spec='',prob_function=function(degs){degs/sum(degs)},create.files=T){
  # draw module sizes
  mod_sizes<-sample(mod_lims[1]:mod_lims[2],nmods,replace=T)
  # induce modules by preferential attachment (or other probability function) off of random seed nodes
  modules<-list()
  for (m in 1:length(mod_sizes)){
    mod_size<-mod_sizes[m]
    # get seed node and initialize module of size 1
    seed<-as.numeric(sample(V(g),1))
    currmod<-c(seed)
    sel<-c(NA)
    # iteratively add neighbors until either desired module size is reached or
    # no more neighbors are available
    while(length(currmod) < mod_size && length(sel)>0){
      # determine potential "attachment points", i.e. all neighbors of nodes in the current module
      attach_points<-neighborhood(g,1,currmod,'all')
      # remove attachment points that don't have any neighbors not already contained in the module
      tosel<-sapply(1:length(attach_points),function(i){length(setdiff(attach_points[[i]],currmod))>0})
      selmod<-currmod[tosel]
      # randomly pick on of the attachment points 
      # (with a given probability function depending on node degree)
      if (length(selmod)>0){
        attsel<-sample(1:length(selmod),1,prob=prob_function(as.numeric(igraph::degree(g,selmod,'all'))))
        sel<-attach_points[tosel][[attsel]]
        sel<-setdiff(sel,currmod)
        # add a randomly selected neighbor of the picked attachment point
        if (length(sel)>1){
          to_add<-sample(x=sel,size=1)
        }else {
          if (length(sel)==1){
            to_add<-sel	
          } else{
            sel<-c()
            to_add<-sel
          }
        }
        currmod<-union(currmod,to_add)
      }else{
        sel<-c()	
      }
    }
    modules[[m]]<-currmod
  }
  
  # draw gene p-values
  pval_file<-sprintf('./modsim%s.rtb',spec)
  mod_nodes<-unique(unlist(modules))
  # define scores with enrichment of module genes at high end
  all_scores<-runif(length(V(g)))
  for (mod in modules){
    enrich<-sample(mod,ceiling(length(mod)*mod_enrich_perc),replace=F)
    all_scores[enrich]<-runif(length(enrich))*pval_scaling	
  }
  names(all_scores)<-V(g)$name
  o<-order(all_scores)
  all_types<-rep('BG',length(V(g)))
  names(all_types)<-V(g)$name
  mod_affect<-sapply(mod_nodes,function(mn){
    bla<-as.character(which(sapply(1:length(modules),function(i){
      mn %in% modules[[i]]})));
    sprintf('M%s',paste(bla,collapse=','))
  })
  all_types[mod_nodes]<-mod_affect
  pvals<-data.frame(P.Value=all_scores[V(g)$name],NodeType=all_types[V(g)$name])
  rownames(pvals)<-V(g)$name
  pvals<-pvals[o,]
  if (create.files){
    write.table(pvals,pval_file,sep='\t',quote=F)
  }
  list(mods=modules,pvals=pvals,pvalfile=pval_file)
}


# compute ES score of a single star and return only the score
ptilde.score<- function(gene.list.scores, gs.idx){
	scores<-sort(gene.list.scores[gs.idx])
	m<-length(scores)
	ptildes<-sapply(1:m,function(k){log10(pbinom(k-1,m,scores[k],lower.tail=F,log.p=F))})
	min(ptildes)
}

# compute ES score of a single star and return:
# the score, the position k at which a minimum p~k occurs, the size of the star and the selected neighbors with their scores
ptilde.score.ext<- function(gene.list.scores, gs.idx){
	scores<-sort(gene.list.scores[gs.idx])
	m<-length(scores)
	ptildes<-sapply(1:m,function(k){log10(pbinom(k-1,m,scores[k],lower.tail=F,log.p=F))})
	ksel<-as.numeric(which.min(ptildes))
	c(min(ptildes),ksel,m,scores[ksel])
}

# wrapper function to compute (parallelized) ES scores for all stars mentioned in gsc.idx
compute.ptilde<-function(gsc.idx,gene.list.scores){
    gs=NULL # evade R CMD check notes for undefined "global" variables
	ptildes<-foreach(gs=1:length(gsc.idx),.combine='c') %dopar% {
		ptilde.score(gene.list.scores, gsc.idx[[gs]])
	}
	names(ptildes)<-names(gsc.idx)
	ptildes
}

# wrapper function to compute ES scores plus additional info for all stars mentioned in gsc.idx
compute.ptilde.ext<-function(gsc.idx,gene.list.scores){
	ptildes<-matrix(nrow=0,ncol=4)
	for (gs.idx in gsc.idx){
		ptildes<-rbind(ptildes,ptilde.score.ext(gene.list.scores, gs.idx))
	}
	rownames(ptildes)<-names(gsc.idx)
	colnames(ptildes)<-c('Ptilde','k','m','pk')
	ptildes
}


# returns a list of neighborhoods (as vectors of node names)
get.nhs<-function(g,gene.scores){
    nidx=NULL # evade R CMD check notes for undefined "global" variables
	node_idcs<-V(g)
	gsc.idx<-list();
	gsc.idx<-foreach (nidx=1:length(node_idcs),.packages = c('igraph')) %dopar% {
	  vc<-node_idcs[as.numeric(nidx)]
	  tmp<-unique(c(vc,igraph::neighbors(graph=g, v=V(g)[vc],mode='all')));
	  intersect(igraph::V(g)[tmp]$name,names(gene.scores))
	}
	names(gsc.idx)<-V(g)$name;
	gsc.idx<-gsc.idx[sapply(gsc.idx,length)>0]
	gsc.idx
}

# read in gene order determined e.g. by limma
# extract genes present in graph <g>
# return graph reduced to all genes mentioned in ranking file
reduce.graph<-function(g,infile,add.scored.genes=F,keep.nodes.without.scores=F,verbose=F){
	tmp<-read.table(infile,sep='\t',quote='',comment.char='',stringsAsFactors=F)
	if (is.null(rownames(tmp))|| rownames(tmp)[1]=='1'){
		gene.list.char<-tmp[,1]	
	}else{
		gene.list.char<-rownames(tmp)
	}
	if (!keep.nodes.without.scores){
	  matched_genes<-intersect(V(g)$name,gene.list.char)
	  g2<-induced.subgraph(g,matched_genes)
	}else{g2<-g}
	if (add.scored.genes){
	  add_genes<-setdiff(gene.list.char,V(g2)$name)
	  g2<-add_vertices(g2,length(add_genes),name=add_genes)
	 }
	if (verbose){
	  print(sprintf('Adapted graph with %i nodes %i edges to %i genes in ranking and %i edges',length(V(g)),length(E(g)),length(V(g2)),length(E(g2))))
	}
	g2
}

# determine genes with scores present in graph <g>
# return graph reduced to all genes contained in <gene.list.scores>
reduce.graph.fromdata<-function(g,gene.list.scores,add.scored.genes=F,keep.nodes.without.scores=F,verbose=F){
  gene.list.char<-names(gene.list.scores)
  if (!keep.nodes.without.scores){
    matched_genes<-intersect(V(g)$name,gene.list.char)
    g2<-induced.subgraph(g,matched_genes)
  }else {g2<-g}
  if (add.scored.genes){
    add_genes<-setdiff(gene.list.char,V(g2)$name);
    g2<-add_vertices(g2,length(add_genes),name=add_genes)
  }
  if (verbose){
    print(sprintf('Adapted graph with %i nodes %i edges to %i genes in ranking and %i edges',length(V(g)),length(E(g)),length(V(g2)),length(E(g2))))
  }
  g2
}

# retrieve gene pvalues for all genes present in graph g from tab-separated file
# possible input file formats are:
# 1) 2 columns, no column or row names: first column = gene names / second column = pvalues
# 2) output format of limma, including row and column names: rownames = gene names / column "P.Value" = pvalues
get.gene.scoring<-function(g,infile){
	tmp<-read.table(infile,sep='\t',quote='',comment.char='',stringsAsFactors=F)
	if (is.null(rownames(tmp)) || is.null(colnames(tmp)) || colnames(tmp)[1]=='V1' || rownames(tmp)[1]=='1'){
		gene.list.char<-tmp[,1]
		gene.list.scores<-as.numeric(tmp[,2])
	}else{
		gene.list.char<-rownames(tmp)
		gene.list.scores<-as.numeric(tmp[,'P.Value'])
	}
	names(gene.list.scores)<-gene.list.char
	gene.list.scores<-gene.list.scores[V(g)$name]
	gene.list.scores[gene.list.scores==1]<- (1 - 1e-5)
	gene.list.scores
}


# extract the local subnetwork around a given protein and write to sif file
# (can be loaded in Cytoscape)
write.ls.to.sif<-function(prot_id,LEANres,outfile){
    i=NULL# evade R CMD check notes for undefined "global" variables
    restab<-LEANres$restab;g<-LEANres$indGraph
    gsc.idx<-LEANres$nhs;gene_list_scores<-LEANres$gene.scores
    nh.scores<-sort(gene_list_scores[as.character(gsc.idx[[prot_id]])])
    k<-restab[prot_id,'k']
    sel_neighbors<-names(nh.scores)[1:k]
    left_neighbors<-setdiff(names(nh.scores),sel_neighbors)
    if (length(sel_neighbors)>0){
        sel_cons<-foreach(i=1:length(sel_neighbors),.combine=rbind) %do% c(prot_id,'selected',sel_neighbors[i])
        if (length(left_neighbors)>0){
            left_cons<-foreach(i=1:length(left_neighbors),.combine=rbind) %do% c(prot_id,'interacts',left_neighbors[i])
            siftab<-rbind(sel_cons,left_cons)
        }else{
            siftab<-sel_cons    
        }
    }else{
        if (length(left_neighbors)>0){
            siftab<-foreach(i=1:length(left_neighbors),.combine=rbind) %do% c(prot_id,'interacts',left_neighbors[i])
        }else{
            print(sprintf('Warning! No valid neighbors found for protein %s',prot_id))    
            return(1)
        }
    }
    write.table(siftab,outfile,sep='\t',row.names=F,col.names=F,quote=F)
}

# extract info about the local subnetwork around a given protein
get.ls.info<-function(prot_id,LEANres){
  i=NULL# evade R CMD check notes for undefined "global" variables
  restab<-LEANres$restab;g<-LEANres$indGraph
  gsc.idx<-LEANres$nhs;gene_list_scores<-LEANres$gene.scores
  nh.scores<-sort(gene_list_scores[as.character(gsc.idx[[prot_id]])])
  k<-restab[prot_id,'k']
  info.tab<-data.frame(ID=names(nh.scores),input.ps=nh.scores,selected=1:length(nh.scores)<=k)
  if (length(nh.scores)<1){
      print(sprintf('Warning! No valid neighbors found for protein %s',prot_id))    
      return(1)
  }
  info.tab
}

# MAIN function
run.lean<-function(ranking,network,ranked=F,add.scored.genes=F,keep.nodes.without.scores=F,verbose=F,n_reps=10000,bootstrap=F,ncores=NULL){
    mi=n=NULL # evade R CMD check notes for undefined "global" variables
  # try to initialize parallel backend
  register.backend(cores=ncores)
  if (is.character(network)){
	  if (verbose){
		  print('## Parsing network file...')
		  print(system.time(g <- sif2graph(network)))
	  }else{
		  g <- sif2graph(network)
	  }
  }else{
    if (is.igraph(network)) g<-network
    else{
      stop('Provided network is neither a valid file name nor an igraph object. Exiting')
    }
  }
	
	# adapt graph to nodes appearing in gene ranking
	if (verbose){
		print('## Adapting graph to scoring genes...')
	  if (is.character(ranking)){
		  print(system.time(g2<-reduce.graph(g,ranking,add.scored.genes,keep.nodes.without.scores,verbose)))
	    gene.list.scores<-get.gene.scoring(g2,ranking)
	  }else{
	    print(system.time(g2<-reduce.graph.fromdata(g,ranking,add.scored.genes,keep.nodes.without.scores,verbose)))  
	    gene.list.scores<-ranking
	  }
	} else{
	  if (is.character(ranking)){
	    g2<-reduce.graph(g,ranking,add.scored.genes,keep.nodes.without.scores,verbose)
	    gene.list.scores<-get.gene.scoring(g2,ranking)
	  }else{
	    g2<-reduce.graph.fromdata(g,ranking,add.scored.genes,keep.nodes.without.scores)
	    gene.list.scores<-ranking
	  }
	}
	
	# transform scores into uniform pval dists if <ranked>==T
	if (ranked){
		N<-length(gene.list.scores)
		gene.list.ranks<-rank(gene.list.scores,ties.method='random')
		gene.list.scores2<-sapply(gene.list.ranks,function(r){r/N})
		gene.list.scores<-gene.list.scores2-(min(gene.list.scores2)/2)
	}
	
	# define neighborhoods as gene sets
	if (verbose){
		print('## Defining neighborhoods as gene sets...')
		print(system.time(gsc.idx<-get.nhs(g2,gene.list.scores)))
		print(sprintf('Found %i neighborhoods with at least one measured input score',length(gsc.idx)))
	}else{
		gsc.idx<-get.nhs(g2,gene.list.scores)
	}
	g2m<-sapply(names(gsc.idx),function(x)length(gsc.idx[[x]]))
	all_ms<-sort(unique(as.numeric(g2m)),decreasing = TRUE)
	N<-length(gene.list.scores)
	if (verbose){
		print(sprintf('## Computing random p~ background dists with n=%i...',n_reps))
		print(sprintf('Number of different observed gene set sizes: %i',length(all_ms)))
		print(system.time(bgs<-foreach (mi=1:length(all_ms),.combine=rbind) %dopar% {
			sapply(1:n_reps,function(i){ptilde.score(gene.list.scores, sample(1:N,all_ms[mi],replace=bootstrap))})
		}))
	}else{
		bgs<-foreach (mi=1:length(all_ms),.combine=rbind) %dopar% {
			sapply(1:n_reps,function(i){ptilde.score(gene.list.scores, sample(1:N,all_ms[mi],replace=bootstrap))})}
	}
	rownames(bgs)<-as.character(all_ms)
	colnames(bgs)<-c()
	
	# calculate mean/sd of bg dists to be used in the calculation of z-scores later
	zpars<-list()
	zscore<-function(xs,pars){(xs-pars$mu)/pars$sigma}
	for (m in all_ms){
		mchar<-as.character(m)
		mu<-mean(bgs[mchar,])
		sigma<-sd(bgs[mchar,])
		zpars[[mchar]]<-list(mu=mu,sigma=sigma)
	}
	
	# compute enrichment for each gene neighborhood
	if (verbose){
		print('# Computing enrichment scores for each gene neighborhood...')
		print(system.time(ptilde.ext<-compute.ptilde.ext(gsc.idx,gene.list.scores)))
	}else{
		ptilde.ext<-compute.ptilde.ext(gsc.idx,gene.list.scores)
	}
	
	# create results table
	if (verbose){print('Compiling result table...')}
	tac<-system.time(comps_tab<-foreach(n=1:dim(ptilde.ext)[1],.combine=rbind) %do% {
		mchar<-as.character(g2m[n])
		true_ps<-ptilde.ext[n,'Ptilde']
		bg_es <- as.numeric(bgs[mchar,])
		mean_bg <- mean(bg_es,na.rm=T)
		z<-zscore(true_ps,zpars[[mchar]])
		c(true_ps,mean_bg,ptilde.ext[n,'k'],g2m[n],ptilde.ext[n,'pk'],sum(bg_es<=true_ps)/n_reps,z)
	})
	if (verbose){print(tac)}
	dimnames(comps_tab)<-list(names(gsc.idx),c('Ptilde','mean_bg_p','k','m','pk','pstar','z.score'))	
	comps_tab<-cbind(comps_tab,p.adjust(((comps_tab[,'pstar']*n_reps)+1)/(n_reps+1),'BH'))
	colnames(comps_tab)[dim(comps_tab)[2]]<-'PLEAN'
	
	# compute same scores based on random perm of gene scores
	if (verbose){
	  print('# Computing enrichment scores on permuted scores for each gene neighborhood...')
	  perm<-sample(1:length(gene.list.scores),length(gene.list.scores),replace=F)
	  gene.scores.permed<-gene.list.scores[perm]
	  names(gene.scores.permed)<-names(gene.list.scores)
	  print(system.time(ptilde.ext.rand<-compute.ptilde.ext(gsc.idx,gene.scores.permed)))
	}else{
	  perm<-sample(1:length(gene.list.scores),length(gene.list.scores),replace=F)
	  gene.scores.permed<-gene.list.scores[perm]
	  names(gene.scores.permed)<-names(gene.list.scores)
	  ptilde.ext.rand<-compute.ptilde.ext(gsc.idx,gene.scores.permed)
	}
	
	# create result table
	if (verbose){print('Compiling result table...')}
	tac<-system.time(comps_tab.rand<-foreach(n=1:dim(ptilde.ext.rand)[1],.combine=rbind) %do% {
		mchar<-as.character(g2m[n])
		true_ps<-ptilde.ext.rand[n,'Ptilde']
		bg_es <- as.numeric(bgs[mchar,])
		mean_bg <- mean(bg_es,na.rm=T)
		z<-zscore(true_ps,zpars[[mchar]])
		c(true_ps,mean_bg,ptilde.ext.rand[n,'k'],g2m[n],ptilde.ext.rand[n,'pk'],sum(bg_es<=true_ps)/n_reps,z)	
	})
	if (verbose){print(tac)}
	dimnames(comps_tab.rand)<-list(names(gsc.idx),c('Ptilde','mean_bg_p','k','m','pk','pstar','z.score'))	
	comps_tab.rand<-cbind(comps_tab.rand,p.adjust(((comps_tab.rand[,'pstar']*n_reps)+1)/(n_reps+1),'BH'))
	colnames(comps_tab.rand)[dim(comps_tab.rand)[2]]<-'PLEAN'
	
	list(restab=comps_tab,randtab=comps_tab.rand,indGraph=g2,nhs=gsc.idx,gene.scores=gene.list.scores)
}

# run.lean.fromdata<-function(gene.list.scores,g,ranked=F,add.scored.genes=F,keep.nodes.without.scores=F,verbose=F,n_reps=10000,bootstrap=F,ncores=NULL){
#   mi=n=NULL # evade R CMD check notes for undefined "global" variables
#   
#   # try to initialize parallel backend
#   register.backend(cores=ncores)
#   
# 	# transform scores into uniform pval dists if <ranked>==T
# 	if (ranked){
# 		N<-length(gene.list.scores)
# 		gene.list.ranks<-rank(gene.list.scores,ties.method='random')
# 		gene.list.scores2<-sapply(gene.list.ranks,function(r){r/N})
# 		gene.list.scores<-gene.list.scores2-(min(gene.list.scores2)/2)
# 	}
# 	
# 	# define neighborhoods as gene sets
# 	if (verbose){
# 		print('## Defining neighborhoods as gene sets...')
# 		print(system.time(gsc.idx<-get.nhs(g2,gene.list.scores)))
# 		print(sprintf('Found %i neighborhoods with at least one measured input score',length(gsc.idx)))
# 	}else{
# 		gsc.idx<-get.nhs(g2,gene.list.scores)
# 	}
# 	g2m<-sapply(names(gsc.idx),function(x)length(gsc.idx[[x]]))
# 	all_ms<-sort(unique(as.numeric(g2m)),decreasing = TRUE)
# 	N<-length(gene.list.scores)
# 	if (verbose){
# 		print(sprintf('## Computing random p~ background dists with n=%i...',n_reps))
# 		print(sprintf('Number of different observed gene set sizes: %i',length(all_ms)))
# 		print(system.time(bgs<-foreach (mi=1:length(all_ms),.combine=rbind,.packages = c("LEANR","igraph")) %dopar% {
# 			sapply(1:n_reps,function(i){ptilde.score(gene.list.scores, sample(1:N,all_ms[mi],replace=bootstrap))})
# 		}))
# 	}else{
# 		bgs<-foreach (mi=1:length(all_ms),.combine=rbind,.packages = c("LEANR","igraph")) %dopar% {
# 			sapply(1:n_reps,function(i){ptilde.score(gene.list.scores, sample(1:N,all_ms[mi],replace=bootstrap))})}
# 	}
# 	rownames(bgs)<-as.character(all_ms)
# 	colnames(bgs)<-c()
# 	
# 	# calculate mean/sd of bg dists to be used in the calculation of z-scores later
# 	zpars<-list()
# 	zscore<-function(xs,pars){(xs-pars$mu)/pars$sigma}
# 	for (m in all_ms){
# 		mchar<-as.character(m)
# 		mu<-mean(bgs[mchar,])
# 		sigma<-sd(bgs[mchar,])
# 		zpars[[mchar]]<-list(mu=mu,sigma=sigma)
# 	}
# 	
# 	# compute enrichment for each gene neighborhood
# 	if (verbose){
# 		print('# Computing enrichment scores for each gene neighborhood...')
# 		print(system.time(ptilde.ext<-compute.ptilde.ext(gsc.idx,gene.list.scores)))
# 	}else{
# 		ptilde.ext<-compute.ptilde.ext(gsc.idx,gene.list.scores)
# 	}
# 	
# 	# create results table
# 	if (verbose){print('Compiling result table...')}
# 	tac<-system.time(comps_tab<-foreach(n=1:dim(ptilde.ext)[1],.combine=rbind) %do% {
# 		mchar<-as.character(g2m[n])
# 		true_ps<-ptilde.ext[n,'Ptilde']
# 		bg_es <- as.numeric(bgs[mchar,])
# 		mean_bg <- mean(bg_es,na.rm=T)
# 		z<-zscore(true_ps,zpars[[mchar]])
# 		c(true_ps,mean_bg,ptilde.ext[n,'k'],g2m[n],ptilde.ext[n,'pk'],sum(bg_es<=true_ps)/n_reps,z)
# 	})
# 	if (verbose){print(tac)}
# 	dimnames(comps_tab)<-list(names(gsc.idx),c('Ptilde','mean_bg_p','k','m','pk','pstar','z.score'))	
# 	comps_tab<-cbind(comps_tab,p.adjust(((comps_tab[,'pstar']*n_reps)+1)/(n_reps+1),'BH'))
# 	colnames(comps_tab)[dim(comps_tab)[2]]<-'PLEAN'
# 	
# 	# compute same scores based on random perm of gene scores
# 	if (verbose){
# 		print('# Computing enrichment scores on permuted scores for each gene neighborhood...')
# 		perm<-sample(1:length(gene.list.scores),length(gene.list.scores),replace=F)
# 		gene.scores.permed<-gene.list.scores[perm]
# 		names(gene.scores.permed)<-names(gene.list.scores)
# 		print(system.time(ptilde.ext.rand<-compute.ptilde.ext(gsc.idx,gene.scores.permed)))
# 	}else{
# 		perm<-sample(1:length(gene.list.scores),length(gene.list.scores),replace=F)
# 		gene.scores.permed<-gene.list.scores[perm]
# 		names(gene.scores.permed)<-names(gene.list.scores)
# 		ptilde.ext.rand<-compute.ptilde.ext(gsc.idx,gene.scores.permed)
# 	}
# 	
# 	# create randomized result table
# 	if (verbose){print('Compiling result table...')}
# 	tac<-system.time(comps_tab.rand<-foreach(n=1:dim(ptilde.ext.rand)[1],.combine=rbind) %do% {
# 		mchar<-as.character(g2m[n])
# 		true_ps<-ptilde.ext.rand[n,'Ptilde']
# 		bg_es <- as.numeric(bgs[mchar,])
# 		mean_bg <- mean(bg_es,na.rm=T)
# 		z<-zscore(true_ps,zpars[[mchar]])
# 		c(true_ps,mean_bg,ptilde.ext.rand[n,'k'],g2m[n],ptilde.ext.rand[n,'pk'],sum(bg_es<=true_ps)/n_reps,z)	
# 	})
# 	if (verbose){print(tac)}
# 	dimnames(comps_tab.rand)<-list(names(gsc.idx),c('Ptilde','mean_bg_p','k','m','pk','pstar','z.score'))	
# 	comps_tab.rand<-cbind(comps_tab.rand,p.adjust(((comps_tab.rand[,'pstar']*n_reps)+1)/(n_reps+1),'BH'))
# 	colnames(comps_tab.rand)[dim(comps_tab.rand)[2]]<-'PLEAN'
# 	
# 	list(restab=comps_tab,randtab=comps_tab.rand,indGraph=g2,nhs=gsc.idx,gene.scores=gene.list.scores)
# }
