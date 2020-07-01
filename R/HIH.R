library("randomForestSRC")
library("stats")
library("spatstat")
library("RColorBrewer")
#############################################################################################
### PMD based weight method to update subRF, iteration always set to 1 in our simulations ###
#############################################################################################
wt.itr <- function(formula, data, pmd.initial,
                   obj.initial = NULL,
                   wt = function(btpmd,digpmd){log(1/btpmd/digpmd)}, iteration = 3){
  pmd <- pmd.initial
  bd.kink <- function(pmd){
    o.max <- nrow(pmd)
    sapply(seq(o.max),function(x)
      which.max(sort(pmd[x,-x],decreasing = FALSE)-
                  smooth.spline(seq(o.max-1),sort(pmd[x,-x],decreasing = FALSE))$y))
  }
  lapply(1:iteration,function(i){
    sim.intcopy <- pmd
    k <- bd.kink(sim.intcopy)
    btw <- mapply(function(x,k)
      mean(sort(sim.intcopy[x,-x], decreasing = FALSE)[1:k], na.rm=TRUE)
      , seq(nrow(sim.intcopy)), k)
    digw <- diag(pmd)
    
    if (is.null(names(diag))){id <- 1:ncol(pmd)} else{
      yvar_name=colnames(data)[which(sapply(1:ncol(data),function(i) {
        grepl(colnames(data)[i],as.character(formula)[[2]])})==TRUE)]
      id <- match(colnames(data[,which(colnames(data)!=yvar_name)]),names(digw))}
    btpmd <- btw[id]; digpmd <- digw[id]
    w <- wt(btpmd,digpmd)
    
    if (iteration!=1){
      sim.out <- rfsrc(formula = formula, data, predictorWt = w, statistics = F,
                       mtry = obj.initial$mtry, nodesize = obj.initial$nodesize)
      pmd <<- find.interaction(sim.out, method="maxsubtree", sorted = F, verbose = F) }
    out <- cbind(pmd, digw, btw, w)
    colnames(out) <- c(paste("pmd", rownames(pmd)),"digpmd", "btpmd", "wt")
    data.frame(out)
  })
}

##################
### Grow subRF ###
##################
subRF <- function(formula, data, w0, subvars = function(data){ceiling((ncol(data)-1)/5)}, n.RF = 5, wtRF = T){
  var.columns <- (1:ncol(data))[which(colnames(data) %in% names(w0))]
  w0 <- as.vector(w0)
  w0[which(w0 <= 0)] <- 0
  size <- subvars(data)
  lapply(1:n.RF, function(i){
    if (size == length(var.columns)) {var.pt <- var.columns} else {
      set.seed(i)
      var.pt <- sample(x = var.columns,size = size,replace = FALSE,prob = w0)}
    subdata <- data[,c(var.pt,setdiff(1:ncol(data),var.columns))]
    wts1 <- w0[var.pt]
    if (wtRF == T){
      sim.out <- rfsrc(formula, data = subdata, predictorWt = wts1)} else {
        turn.o <- tune(formula, data = subdata)
        sim.out <- rfsrc(formula,data = subdata,
                         mtry = turn.o$rf$mtry, nodesize = turn.o$rf$nodesize)}
    pmd <- find.interaction(sim.out, method = "maxsubtree", sorted = F, verbose = F)
    err.sub <- if (sim.out$family=="regr"){
      err.sub<-round(100*(1-sim.out$err.rate[1000]/var(sim.out$yvar,na.rm = TRUE)), 2)
    }
    list(subdata = subdata, pmd = pmd, subRF.o = sim.out, err.sub = err.sub)
  })
}

##########################
### Variable Selection ###
##########################
var.sel <- function(formula, data, subvars = function(data){ceiling((ncol(data)-1)/5)},
                    w.initial = "vimp", max.var=50,
                    wt = function(btpmd,digpmd){log(1/btpmd/digpmd)},
                    itrSub = 5, wtSub = T,
                    itrWt = 1, subQtl = function(btpmd,digpmd){which(btpmd < quantile(btpmd, probs = 0.1))}, verbose = TRUE,
                    obj = NULL){
  if (w.initial == "vimp"){
    w.initial <- function(formula, data, obj = NULL){
      if (is.null(obj)){
        turn.o <- tune(formula, data)
        sim.out <- rfsrc(formula, data, importance = "random",
                         mtry = turn.o$rf$mtry, nodesize = turn.o$rf$nodesize)}
      else {sim.out <- obj}
      if (sim.out$family == "class"){
        list(sim.out$err.rate, pmax(sim.out$importance[,1],0))
      } else {
        list(round(100*(1-sim.out$err.rate[1000]/var(sim.out$yvar,na.rm = TRUE)),2),
             pmax(sim.out$importance,0))
      }}
    w0 <- w.initial(formula, data, obj = obj)[[2]]
  } else if (w.initial == "md"){
    w.initial <- function(formula, data, obj = NULL){
      if (is.null(obj)){
        turn.o <- tune(formula, data)
        sim.out <- rfsrc(formula,data,importance = FALSE,
                         mtry = turn.o$rf$mtry, nodesize = turn.o$rf$nodesize)}
      else {sim.out <- obj}
      1/(max.subtree(sim.out, sub.order = FALSE, max.order = 1)$order+0.001)}
    w0 <- w.initial(formula, data, obj = obj)
  } else { w0 <- w.initial}
  
  subchain0 <- subRF(formula, data, w0 = w0, subvars = subvars, n.RF = itrSub, wtRF = wtSub)
  chainwts <- lapply(1:itrSub, function(i){
    if (verbose == TRUE) {cat("------", i,"/",itrSub, "------ \n")}
    wt.itr(formula = formula,data = subchain0[[i]]$subdata, pmd.initial = subchain0[[i]]$pmd,
           obj.initial = subchain0[[i]]$subRF.o,
           wt = wt, iteration = itrWt)
  })
  var.sl.list <- lapply(1:itrSub, function(i){
    names <- rownames(chainwts[[i]][[itrWt]])
    id <- subQtl(chainwts[[i]][[itrWt]]$btpmd,chainwts[[i]][[itrWt]]$digpmd)
    list(sl.w = chainwts[[i]][[itrWt]]$wt[id], sl.var = names[id])
  })
  #var.sl <- unique(unlist(lapply(1:itrSub, function(i){var.sl.list[[i]]$sl.var})))
  varl<-lapply(lapply(var.sl.list,`[`,2),`[[`,1)
  if (length(unique(unlist(varl)))<max.var){
    var.sl<-unique(unlist(varl))
  }
  else {
    freql<-table(unlist(varl))[order(table(unlist(varl)),decreasing=TRUE)]
    freql.stepsum<-sapply(unique(freql),function(i){
      sum(freql>=i)
    })
    freql.thresh<-unique(freql)[which(abs(freql.stepsum-max.var)==min(abs(freql.stepsum-max.var)))]
    var.sl<-names(freql)[which(freql>=freql.thresh)]
  }
  list(var.sl = var.sl,var.sl.list = var.sl.list, itrwt = chainwts,
       err.init = w.initial(formula,data,obj = obj)[[1]], 
       err.sub = unlist(lapply(1:itrSub, function(x){subchain0[[x]]$err.sub})))
}

##########################################################################
### Calculate PMD VIMP for all interactions, i.e, select 10 variables, ###
### minimum order is defined by function "nterms" as 3, then all inter-###
### actions will be (10 choose 3)+(10 choose 2)+(10 chooose 1)         ###
##########################################################################
pmdhvp <- function(obj, inter, verbose = T){
  pmd <- xx <- find.interaction(obj, sorted = F , method = "maxsubtree", verbose = F)
  name <- rownames(pmd)
  pmd[lower.tri(xx)] <- (xx[lower.tri(xx)]+xx[upper.tri(xx)])/2
  pmd[upper.tri(xx,diag=T)] <- 0
  pmd = pmd + t(pmd) + diag(diag(xx))
  
  k <- length(inter)
  pmdvp <- matrix(unlist(lapply(1:k,function(i){
    if (((i%%100 == 0)|(i == k))&(verbose == T)){cat("\r","----",i,"/",k,100*(i/k),"%----")}
    
    mean(apply(as.matrix(combn(inter[[i]],2)),2,function(x){
      pmd[x[1],x[2]] }))  }) ),k ,1)
  rownames(pmdvp) <- unlist(lapply(1:k,function(i){ paste(name[inter[[i]]], collapse = "_") }))
  list(pmdvp = pmdvp, pmd = pmd)
}

##########################################################################
### This function is particularly used after selecting top 150 pmdvimp ###
### ranked interactions. It is to get the HIVIMP and normalized HIVIMP ###
### for the 150 interactions.                                          ###
##########################################################################

hvp <- function(joint, obj, importance = "permute",
                block = 1, choice = c("hivimp", "acuvimp"))
{ allobj <- obj
x.var <- allobj$xvar.names
if ( any(choice == "hivimp")  ){
  obj <- vimp(allobj,allobj$xvar.names,importance = importance,joint = F, seed = -1, block = block) }
o <- do.call(rbind,lapply(1:nrow(joint),function(i){
  
  if ( any(choice == "hivimp")  ){
    if (allobj$family == "class"){
      joinvimp <- vimp(allobj,x.var[joint[i,]],importance = importance,joint = T, seed = -1,block=block)$importance[1]
      hivimp <- (joinvimp - sum(obj$importance[joint[i,],1]))}else{
        joinvimp <- vimp(allobj,x.var[joint[i,]],importance = importance,joint = T, seed = -1,block=block)$importance
        hivimp <- (joinvimp - sum(obj$importance[joint[i,]]))
      }
    joinvimp <- joinvimp
    HIvimp <- hivimp
    nmHIvimp <- hivimp/ncol(joint)
    if (length(choice) == 1) {acuvimp <- NA}
  }
  if ( any(choice == "acuvimp")  ){
    percond <- condY(allobj)
    fml <- as.formula(paste("y~",paste(allobj$xvar.names[joint[i,]],collapse = "+")))
    joiny <- rfsrc(fml, data = percond$daty,importance = F, seed = -1,
                   mtry = allobj$mtry, nodesize = allobj$nodesize)$predicted.oob
    lmo <- lm(fml,data = data.frame(y = log(abs(joiny)), log(abs(percond$cond))))
    acuvimp <- sum((lmo$residuals)^2)/lmo$df.residual
    if (length(choice) == 1) {joinvimp <- HIvimp <- nmHIvimp <- NA}
  }
  
  out <- data.frame(joinvimp = joinvimp,
                    HIvimp = HIvimp,
                    nmHIvimp = nmHIvimp,
                    acuvimp = acuvimp)
  out <- out[,!is.na(out)]
  out <- as.data.frame(out)
  rownames(out) <- paste(x.var[joint[i,]],collapse ="_")
  out
}))
o
}

##############################################################
### Randomly sample interactions to approximate normalized ###
### HIVIMP distribution density. Then use this density to  ###
### calculate quantile scores for selected interactions.   ###
##############################################################

hvp1 <- function(joint1, obj, importance = "permute",
                 block = 1, choice = c("hivimp", "acuvimp"))
{ allobj <- obj
x.var <- allobj$xvar.names
if ( any(choice == "hivimp")  ){
  obj <- vimp(allobj,allobj$xvar.names,importance = importance,joint = F, seed = -1, block = block) }
o1 <- do.call(rbind,lapply(1:nrow(joint1),function(i){
  
  if ( any(choice == "hivimp") ){
    if (allobj$family == "class"){
      joinvimp <- vimp(allobj,x.var[joint1[i,]],importance = importance,joint = T, seed = -1,block=block)$importance[1]
      hivimp <- (joinvimp - sum(obj$importance[joint1[i,],1]))}else{
        joinvimp <- vimp(allobj,x.var[joint1[i,]],importance = importance,joint = T, seed = -1,block=block)$importance
        hivimp <- (joinvimp - sum(obj$importance[joint1[i,]]))
      }
    nmHIvimp <- hivimp/ncol(joint1)
  }
  nmHIvimp
}))
o1
}

##############################################################
### Function to determine up to which order we will search ###
### for the interactions.                                  ###
##############################################################

nterms <- function(obj){
  if (nrow(obj$xvar)>=(100*ncol(obj$xvar))){
    stat.obj <- stat.split(obj)
    list <- rapply(stat.obj,classes="matrix",how="list",
                   f=function(x) x[,6,drop=FALSE])
    listidx <-seq(length(list))
    tree.depths <- sapply(listidx,function(x) max((unlist(list[[x]]))))
    Dt <- mean(tree.depths)
    All.mdp <- mean(max.subtree(obj, max.order = 0)$order)
    terms <- min(c(round(Dt-All.mdp),ncol(obj$xvar)))
    list(tree.depth = Dt, min.depth = All.mdp, terms = terms)
  } else {
    list(terms = min(c(floor(log2((nrow(obj$xvar))/obj$nodesize)),ncol(obj$xvar)))
    )}
}

################################################
### Main function to get everything together ###
################################################

HIvimp <- function(formula, data, method = "hivimp", terms = NULL, PMDorder = TRUE,
                   importance = "permute", block = 1,
                   subRF = TRUE, subvars = NULL,
                   w.initial = "vimp",
                   wt = function(btpmd,digpmd){log(1/btpmd/digpmd)},
                   itrSub = NULL, wtSub = TRUE,
                   itrWt = 1, subQtl = NULL, max.var = 50,
                   nHIvimp = function(all.cmbn){min(length(all.cmbn),150)},
                   qtlscr=TRUE,
                   verbose = TRUE){
  if (is.null(terms)){
    turn.o <- tune(formula, data)
    Allobj <-  rfsrc(formula, data, importance = "random",statistics = TRUE,
                     mtry = turn.o$rf$mtry, nodesize = turn.o$rf$nodesize)
    terms <- nterms(Allobj)$terms
  }
  
  if (is.null(subQtl)){
    subQtl=function(btpmd,digpmd){
      if (length(btpmd)<=140) {
        which(btpmd < sort(btpmd,decreasing=FALSE)
              [floor(0.5*which.max(sort(btpmd,decreasing = FALSE)-
                                     (((max(btpmd)-min(btpmd))/(length(btpmd)-1))*
                                        (seq(length(btpmd))-1)+min(btpmd)))+
                       0.5*which.max(sort(btpmd,decreasing=FALSE)-
                                       predict(lm(sort(btpmd,decreasing=FALSE)~seq(length(btpmd))))))])}
      else {
        which(btpmd < sort(btpmd,decreasing=FALSE)
              [which.max(sort(btpmd,decreasing = FALSE)-
                           (((max(btpmd)-min(btpmd))/(length(btpmd)-1))*
                              (seq(length(btpmd))-1)+min(btpmd)))])}
    }
  }
  
  if (is.null(subvars)){
    subvars = function(data){
      yvar.num<-length(which(sapply(1:ncol(data),function(i) {
        grepl(colnames(data)[i],as.character(formula)[[2]])})==TRUE)
      )
      if((ncol(data)-yvar.num)<=50) {ncol(data)-yvar.num}
      else {ceiling((ncol(data)-yvar.num)/log(ncol(data)-yvar.num))}}
  }
  
  if (subRF == T){
    if (is.null(itrSub)){itr=terms*ceiling(log10(ncol(data)-length(
      which(sapply(1:ncol(data),function(i) {
        grepl(colnames(data)[i],as.character(formula)[[2]])})==TRUE)
    )))} else{itr=itrSub}
    if (verbose == TRUE) {cat("Fitting", itr, "subRFs \n")}
    sub.x.obj <- var.sel(formula, data,
                         subvars = subvars,
                         w.initial = w.initial, 
                         max.var=max.var, wt = wt,
                         itrSub = itr, wtSub = wtSub,
                         itrWt = itrWt,
                         subQtl = subQtl, verbose = verbose,
                         obj = Allobj)
    yvar_name=colnames(data)[which(sapply(1:ncol(data),function(i) {
      grepl(colnames(data)[i],as.character(formula)[[2]])})==TRUE)]
    data <- data[,c(sub.x.obj$var.sl, yvar_name)]
    turn.o <- tune(formula, data, mtryStart = (ncol(data) - 1))
    allobj <- rfsrc(formula, data = data,importance = importance, seed = -1, mtry = (ncol(data)-1),
                    nodesize = turn.o$rf$nodesize)
  } else {
    sub.x.obj <- list()
    sub.x.obj$var.sl <- sub.x.obj$var.sl.list <- sub.x.obj$itrwt <- NULL
    if (verbose == TRUE) {cat("Fitting RF \n")}
    turn.o <- tune(formula, data, mtryStart = (ncol(data) - 1))
    allobj <- rfsrc(formula, data = data,importance = importance, mtry = (ncol(data)-1), seed = -1,
                    nodesize = turn.o$rf$nodesize)
  }
  cmbn <- lapply(2:terms,function(i){t(combn(1:length(allobj$xvar.names),i))})
  inter <- unlist(lapply(1:length(cmbn), function(i){
    lapply(1:nrow(cmbn[[i]]),function(j){ matrix(cmbn[[i]][j,],1,length(cmbn[[i]][j,])) })
  }), recursive = F )
  if (qtlscr==TRUE){
    inter1 <- inter[sample(1:length(inter),min(length(inter),300))]
  }
  
  if (PMDorder){
    n.vimp <- nHIvimp(inter)
    sub.itr.obj <- pmdhvp(allobj, inter, verbose = verbose)
    pmd <- sub.itr.obj$pmd
    if (length(inter) > n.vimp) {
      if (verbose == TRUE) {cat("Selecting", length(inter), "interactions \n")}
      it.sel <- which(rank(sub.itr.obj$pmdvp) <= n.vimp)
      inter <- inter[it.sel]
    } else {it.sel <- 1:length(inter)}
    pmdvimp.long <- sub.itr.obj$pmdvp
  }else{
    pmd <- pmdvimp.long <- NULL
    n.vimp <- length(inter)
  }
  
  if (verbose == TRUE) {cat("\n Calculating", n.vimp, "High Order Interaction VIMPs \n")}
  itrcts <- as.data.frame(do.call(rbind,lapply(1:length(inter),function(i){
    as.data.frame(hvp(joint = inter[[i]], allobj, importance = importance, block = block, choice = method
    ))
  })))
  if (qtlscr==TRUE){
    pval <- as.data.frame(do.call(rbind,lapply(1:length(inter1),function(i){
      as.data.frame(hvp1(joint1 = inter1[[i]],allobj, importance = importance, block = block, choice = method
      ))
    })))
  }
  if (PMDorder){
    itrcts$pmdvimp <- as.vector(sub.itr.obj$pmdvp[it.sel,1])}
  if (ncol(itrcts) == 1){ colnames(itrcts) = method }
  list(terms = terms,
       interaction = itrcts,
       rf.obj = allobj, pmd = pmd,
       var.sl = sub.x.obj$var.sl,
       var.sl.list = sub.x.obj$var.sl.list,
       itrwt = sub.x.obj$itrwt,
       pmdvimp.long = pmdvimp.long,
       err.init = sub.x.obj$err.init,
       err.sub = sub.x.obj$err.sub,
       err.all = round(100*(1-allobj$err.rate[1000]/var(allobj$yvar,na.rm = TRUE)),2),
       hivimplist = if(qtlscr==FALSE){NULL} else{pval})
}


#####################################################
### Function to output final interactions by 95th ###
### percentile cut. Plot will be generated if     ###
### option is turned on.                          ###
#####################################################

qtlscr.rank <- function(o,data.nm,qtlcut=0.95,plot=TRUE){
  f<-CDF(density(log(abs(o$hivimplist[,1]))))
  froot<-function(f,qtl,lower=-100,upper=100){
    while (f((lower+upper)/2)-qtl!=0){
      if ((f(lower)-qtl)*(f((lower+upper)/2)-qtl)<0){
        upper = (lower+upper)/2
      }
      else{
        lower = (lower+upper)/2
      }
    }
    (lower+upper)/2
  }
  thresh<- exp(froot(f,qtlcut))
  int<-o$interaction
  valnmh<-int[order(-abs(int$nmHIvimp)),]
  valnmhf<-valnmh[which(abs(valnmh[,3])>thresh),]
  seq<-sapply(1:nrow(valnmhf),function(i){
    1-f(log(abs(valnmhf[i,3])))
  })
  seq.copy<-seq
  names(seq.copy)<-rownames(valnmhf)
  cat("Interactions by 95th percentile of normalized high order interaction variable importance density estimation for data",data.nm,"and their quantile scores: \n")
  print(seq.copy)
  if (plot==TRUE){
    rmid<-c()
    for (i in 1:nrow(valnmhf)){
      temp<-c()
      if (i<nrow(valnmhf)){
        for (j in c(1:(i-1),(i+1):nrow(valnmhf))){
          indi<-unlist(strsplit(rownames(valnmhf)[i], split = "_"))
          indj<-unlist(strsplit(rownames(valnmhf)[j], split = "_"))
          if (sum(indi %in% indj)==length(indi)){
            temp<-c(temp,i)
          }
          else {
            temp<-temp
          }
        }
      }
      else if (i==nrow(valnmhf)){
        for (j in c(1:(i-1))){
          indi<-unlist(strsplit(rownames(valnmhf)[i], split = "_"))
          indj<-unlist(strsplit(rownames(valnmhf)[j], split = "_"))
          if (sum(indi %in% indj)==length(indi)){
            temp<-c(temp,i)
          }
          else {
            temp<-temp
          }
        }
      }
      if (length(temp)>0){
        rmid<-c(rmid,i)
      }
      else{
        rmid<-rmid
      }
    }
    kpid<-seq(nrow(valnmhf))
    kpid<-kpid[-rmid]
    grplst<-list()
    for (i in 1:length(kpid)){
      temp<-kpid[i]
      for (j in seq(nrow(valnmhf))[-kpid[i]]){
        indi<-unlist(strsplit(rownames(valnmhf)[kpid[i]], split = "_"))
        indj<-unlist(strsplit(rownames(valnmhf)[j], split = "_"))
        if (sum(indj %in% indi) == length(indj)){
          temp<-c(temp,j)
        }
        else {
          temp<-temp
        }
      }
      grplst[[i]]<-temp
    }
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    color<-sample(col_vector,length(grplst))
    for (i in seq(length(kpid))){
      cidx<-list()
      id<-grplst[[i]]
      for (k in seq(length(id))){
        temp<-c()
        for (j in seq(length(kpid))[-i]){
          if (id[k] %in% unlist(grplst[[j]])==TRUE){
            temp<-c(temp,color[j])
          }
          else {
            temp<-temp
          }
        }
        cidx[[k]]<-c(color[i],temp)
      }
      nbar<-sapply(1:length(cidx),function(i){length(cidx[[i]])})
      if (i==1){
        png(paste(data.nm,i,".png",sep=""),width=1800,height=100+30*length(grplst[[i]]),res=300)
        par(mar=c(0.6, 13.6, 1.1, 2.1))
        b1<-barplot((seq[grplst[[i]]]*(1/nbar))[length(grplst[[i]]):1],
                    main=paste("Quantile Scores of Top Ranked Interactions\n for ",data.nm, " Data",sep=""), 
                    horiz=TRUE,xaxt="n",cex.main=0.6,
                    ylab="",xlab=NULL,xlim=c(0,0.05),las=1,cex.names=0.5,
                    col=unlist(lapply(cidx,`[`,1)),border=NA)
        axis(2,at=b1[1:length(grplst[[i]])],labels=FALSE,tck=-0.5/length(grplst[[i]]))
        if (length(grplst[[i]])==1){
          mtext(text=bquote(bold(.(rownames(valnmhf)[grplst[[i]]][1]))), side=2,
                at=b1[length(grplst[[i]])],
                col="red",las=1,cex=0.6,line=0.6)
        }
        else {
          mtext(text=bquote(bold(.(rownames(valnmhf)[grplst[[i]]][1]))), side=2,
                at=b1[length(grplst[[i]])],
                col="red",las=1,cex=0.6,line=0.6)
          mtext(text=rownames(valnmhf)[grplst[[i]]][2:length(grplst[[i]])], side=2,
                at=b1[(length(grplst[[i]])-1):1],
                col="black",las=1,cex=0.5,line=0.6)
        }
        val<-(seq[grplst[[i]]]*(1/nbar))[length(grplst[[i]]):1]
        if (max(nbar)>1){
          for (l in seq(max(nbar)-1)){
            perc<-sapply(1:length(nbar),function(i){max(min(1-1/nbar[i]*l,1/nbar[i]),0)})
            b1<-barplot((seq[grplst[[i]]]*perc)[length(grplst[[i]]):1],
                        offset=val, add=T, axes=F, axisnames=F, horiz=T, 
                        col=rev(unlist(lapply(cidx,`[`,l+1))),border=NA)
            val<-val+(seq[grplst[[i]]]*perc)[length(grplst[[i]]):1]
          }
        }
      }
      else if (i<length(kpid)){
        png(paste(data.nm,i,".png",sep=""),width=1800,height=40+30*length(grplst[[i]]),res=300)
        par(mar=c(0.6, 13.6, 0.1, 2.1))
        b1<-barplot((seq[grplst[[i]]]*(1/nbar))[length(grplst[[i]]):1],
                    main=NULL,
                    horiz=TRUE,xaxt="n",cex.main=0.8,
                    ylab="",xlab=NULL,xlim=c(0,0.05),las=1,cex.names=0.5,
                    col=unlist(lapply(cidx,`[`,1)),border=NA)
        axis(2,at=b1[1:length(grplst[[i]])],labels=FALSE,tck=-0.5/length(grplst[[i]]))
        if (length(grplst[[i]])==1){
          mtext(text=bquote(bold(.(rownames(valnmhf)[grplst[[i]]][1]))), side=2,
                at=b1[length(grplst[[i]])],
                col="red",las=1,cex=0.6,line=0.6)
        }
        else {
          mtext(text=bquote(bold(.(rownames(valnmhf)[grplst[[i]]][1]))), side=2,
                at=b1[length(grplst[[i]])],
                col="red",las=1,cex=0.6,line=0.6)
          mtext(text=rownames(valnmhf)[grplst[[i]]][2:length(grplst[[i]])], side=2,
                at=b1[(length(grplst[[i]])-1):1],
                col="black",las=1,cex=0.5,line=0.6)
        }
        val<-(seq[grplst[[i]]]*(1/nbar))[length(grplst[[i]]):1]
        if (max(nbar)>1){
          for (l in seq(max(nbar)-1)){
            perc<-sapply(1:length(nbar),function(i){max(min(1-1/nbar[i]*l,1/nbar[i]),0)})
            b1<-barplot((seq[grplst[[i]]]*perc)[length(grplst[[i]]):1],
                        offset=val, add=T, axes=F, axisnames=F, horiz=T, 
                        col=rev(unlist(lapply(cidx,`[`,l+1))),border=NA)
            val<-val+(seq[grplst[[i]]]*perc)[length(grplst[[i]]):1]
          }
        }
      }
      else {
        png(paste(data.nm,i,".png",sep=""),width=1800,height=100+30*length(grplst[[i]]),res=300)
        par(mar=c(1.6, 13.6, 0.1, 2.1))
        b1<-barplot((seq[grplst[[i]]]*(1/nbar))[length(grplst[[i]]):1],
                    main=NULL,
                    horiz=TRUE,xaxt="n",cex.main=0.8,
                    ylab="",xlab="",xlim=c(0,0.05),las=1,cex.names=0.5,
                    col=unlist(lapply(cidx,`[`,1)),border=NA)
        axis(1,at=round(seq(0,0.05,0.05/5),2),labels=FALSE,tck=-0.5/length(grplst[[i]]))
        mtext(text=round(seq(0,0.05,0.05/5),2), side=1,
              at=round(seq(0,max(seq),max(seq)/5),2),las=1,cex=0.6,line=0.1)
        axis(2,at=b1[1:length(grplst[[i]])],labels=FALSE,tck=-0.5/length(grplst[[i]]))
        if (length(grplst[[i]])==1){
          mtext(text=bquote(bold(.(rownames(valnmhf)[grplst[[i]]][1]))), side=2,
                at=b1[length(grplst[[i]])],
                col="red",las=1,cex=0.6,line=0.6)
        }
        else {
          mtext(text=bquote(bold(.(rownames(valnmhf)[grplst[[i]]][1]))), side=2,
                at=b1[length(grplst[[i]])],
                col="red",las=1,cex=0.6,line=0.6)
          mtext(text=rownames(valnmhf)[grplst[[i]]][2:length(grplst[[i]])], side=2,
                at=b1[(length(grplst[[i]])-1):1],
                col="black",las=1,cex=0.5,line=0.6)
        }
        val<-(seq[grplst[[i]]]*(1/nbar))[length(grplst[[i]]):1]
        if (max(nbar)>1){
          for (l in seq(max(nbar)-1)){
            perc<-sapply(1:length(nbar),function(i){max(min(1-1/nbar[i]*l,1/nbar[i]),0)})
            b1<-barplot((seq[grplst[[i]]]*perc)[length(grplst[[i]]):1],
                        offset=val, add=T, axes=F, axisnames=F, horiz=T, 
                        col=rev(unlist(lapply(cidx,`[`,l+1))),border=NA)
            val<-val+(seq[grplst[[i]]]*perc)[length(grplst[[i]]):1]
          }
        }
      }
      dev.off()
    }
    png("legend.png",width=1800,height=50*length(grplst)/2,res=300)
    par(mar=c(0.1,0.1,0.1,0.1))
    frame()
    legend("top", rownames(valnmhf)[unlist(lapply(grplst,`[`,1))],
           fill=color[1:length(grplst)],border=color[1:length(grplst)],
           ncol=2,bty="n",cex=0.4)
    graphics.off()
  }
}






