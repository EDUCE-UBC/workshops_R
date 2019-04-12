# BDTT for non-ultrametric trees
#library(phytools)
#library(betapart)
#library(abind)
#library(Matrix)

#The first function is copied from phytools
getDescendants=function (tree, node, curr = NULL) 
{
  #if (!inherits(tree, "phylo")) 
   # stop("tree should be an object of class \\"phylo\\".")
  if (is.null(curr)) 
    curr <- vector()
  daughters <- tree$edge[which(tree$edge[, 1] == node), 2]
  curr <- c(curr, daughters)
  if (length(curr) == 0 && node <= length(tree$tip.label)) 
    curr <- node
  w <- which(daughters > length(tree$tip.label))
  if (length(w) > 0) 
    for (i in 1:length(w)) curr <- getDescendants(tree, daughters[w[i]], 
                                                  curr)
  return(curr)
}

BDTT=function(similarity_slices,tree,sampleOTUs,onlyBeta=T,metric="jac"){
  
  Betas=lapply(similarity_slices,getBDTT,tree=tree,sampleOTUs=sampleOTUs,onlyBeta=T,metric=metric)
  names(Betas)=similarity_slices
  res=do.call(function(...){abind(...,along=0)},Betas)
  return(res)
  }


getBDTT=function(similarity,tree,sampleOTUs,onlyBeta=T,metric="jac")
{
  New_OTUs_sample_matrix=getNew_OTUs_sample_matrix(similarity=similarity,sampleOTUs=sampleOTUs,tree=tree)
  BetasPA=getBeta(t(New_OTUs_sample_matrix[[1]]),ab=F)[c("jtu","jac"),,]
  #rownames(BetasPA)=paste(rownames(BetasPA),"_Simi_",similarity,sep="")
  BetasAb=getBeta(t(New_OTUs_sample_matrix[[1]]),ab=T)[c("bctu","bc"),,]  
  #rownames(BetasAb)=paste(rownames(BetasAb),"_Simi_",similarity,sep="")
  AllBetas=abind(BetasPA,BetasAb,along = 1)
  
  AllBetas=AllBetas[metric,,]
  
  if (onlyBeta==T){res=list(Beta_Div=AllBetas,NewOTU_table=New_OTUs_sample_matrix[[1]],New_to_old_OTUs=New_OTUs_sample_matrix[[2]])}
  if (onlyBeta==T){res=AllBetas}
  
  return(res)
}

getBeta=function(mat,ab=F)
{
  if (ab==T)
  {   
    h=bray.part(data.frame(mat))
    bctu=as.matrix(h[[1]])
    bcne=as.matrix(h[[2]])
    bc=as.matrix(h[[3]])
    res=abind(bctu,bcne,bc,along=0)
    dimnames(res)[[1]]=c("bctu","bcne","bc")
  }
  if (ab==F)
  {
    mat[mat>0]=1
    h=beta.pair(data.frame(mat), index.family="jaccard")
    hh=beta.pair(data.frame(mat), index.family="sorensen")
    jtu=as.matrix(h[[1]])
    jne=as.matrix(h[[2]])
    jac=as.matrix(h[[3]])
    stu=as.matrix(hh[[1]])
    sne=as.matrix(hh[[2]])
    sor=as.matrix(hh[[3]])  
    res=abind(jtu,jne,jac,stu,sne,sor,along=0)
    dimnames(res)[[1]]=c("jtu","jne","jac","stu","sne","sor")  
  }
  return(res)
}


getNew_OTUs_sample_matrix=function(similarity,sampleOTUs,tree)
{
  CorresMatrix=New_toOld_OTUs(similarity=similarity,tree=tree)
  sample=Matrix(t(sampleOTUs))
  sp=colnames(sample)
  Newsample=t(sample[,sp] %*% CorresMatrix[sp,])
  out=list(NewSample=as.matrix(Newsample),CorrepondanceOTUs=as.matrix(CorresMatrix))
  return(out)
}

getHnode=function(node,tree)
{ 
  tips=1:length(tree$tip.label)
  NH=node.depth.edgelength(tree)
  DescNodes=getDescendants(node=node,tree=tree)
  DescTips=DescNodes[DescNodes%in%tips]
  Hnode=max(NH[DescTips]-NH[node])# take the longest branches to define node height
  return(Hnode)
}

getHnodes=function(tree)
{
  allnodes=sort(unique(c(tree$edge)))
  names(allnodes)=as.character(allnodes)
  tips=1:length(tree$tip.label)
  nodes=allnodes[!allnodes%in%tips]
  
  Hnodes=sapply(nodes,getHnode,tree=tree)
  return(Hnodes)
}

multigetDescendants=function(node,tree){
  alltips=1:length(tree$tip.label)
  nodes=getDescendants(tree=tree,node=node)
  tips=nodes[nodes%in%alltips]
  return(tree$tip.label[tips])
}

New_toOld_OTUs=function(similarity,tree)
{
  Ntips=length(tree$tip.label)
  tips=1:Ntips
  NameTips=tree$tip.label
  Hnodes=getHnodes(tree) #get Nodes height 
  Hbranches=cbind(Hnodes[as.character(tree$edge[,1])],Hnodes[as.character(tree$edge[,2])]) #put them in a matrix of edges
  NodeToCluster=tree$edge[(Hbranches[,1]>similarity)&(Hbranches[,2]<similarity),2] #select the branches whom descandant node to be collapsed
  NodeToCluster=NodeToCluster[!is.na(NodeToCluster)] #remove tips
  
  DescendandTips=lapply(NodeToCluster,multigetDescendants,tree=tree) #get descendant tips
  names(DescendandTips)=as.character(NodeToCluster)
  
  #Lenght of the different categories of OTUs
  N_newOTUS=length(DescendandTips)
  N_collapsedOldOtus= sum(do.call(rbind, lapply(DescendandTips, function(x) length(x))))
  N_totalNewOTUS=Ntips-N_collapsedOldOtus+N_newOTUS
  
  print(paste(similarity," similarity provides ",N_totalNewOTUS," total new OTUs",sep=""))
  
  NewOTUs_OldOTUs_Matrix=NA
  if (N_totalNewOTUS>2) 
  {
    #Names of the different categories of OTUs
    collapsedOldOtus= unlist(DescendandTips)
    UncollapsedOTUs=NameTips[!NameTips%in%collapsedOldOtus]
    
    #Creating and filling New to Old matrix correpondance
    NewOTUSNames=c(UncollapsedOTUs,as.character(NodeToCluster))
    NewOTUs_OldOTUs_Matrix=Matrix(0,ncol=N_totalNewOTUS,nrow=Ntips,dimnames = list(tree$tip.label,NewOTUSNames))
    NewOTUs_OldOTUs_Matrix=Matrix(0,ncol=N_totalNewOTUS,nrow=Ntips,dimnames = list(tree$tip.label,NewOTUSNames))
    
    for (i in UncollapsedOTUs){NewOTUs_OldOTUs_Matrix[i,i]=1} # uncollapsed OTUS
    for (i in as.character(NodeToCluster)){NewOTUs_OldOTUs_Matrix[DescendandTips[[i]],i]=1} #collapsed OTUS
  } else{print("Only 2 OTUS (or less) present at this resolution -- are you sure this is meaningfull?")} 
  
  
  return(NewOTUs_OldOTUs_Matrix)
}