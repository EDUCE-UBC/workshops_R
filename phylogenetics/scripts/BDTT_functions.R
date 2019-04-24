New_toOld_OTUs_castor=function(similarity,tree)
{
  Ntips=length(tree$tip.label)
  NameTips=tree$tip.label
  
  #Collapse tips accorindg to the resolution 
  ColapsedNodes=collapse_tree_at_resolution(tree,similarity)
  DescendingTips=lapply(ColapsedNodes$collapsed_nodes,tree,FUN=function(i,tree){get_subtree_at_node(tree, i)$subtree$tip.label})
  names(DescendingTips)=as.character(ColapsedNodes$collapsed_nodes+Ntips)
  
  
  #Lenght of the different categories of OTUs
  N_newOTUS=length(DescendingTips)
  N_collapsedOldOtus= sum(do.call(rbind, lapply(DescendingTips, function(x) length(x))))
  N_totalNewOTUS=Ntips-N_collapsedOldOtus+N_newOTUS
  
  print(paste(similarity," similarity provides ",N_totalNewOTUS," total new OTUs",sep=""))
  
  NewOTUs_OldOTUs_Matrix=NA
  if (N_totalNewOTUS>2) 
  {
    #Names of the different categories of OTUs
    collapsedOldOtus= unlist(DescendingTips)
    UncollapsedOTUs=NameTips[!NameTips%in%collapsedOldOtus]
    
    #Creating New to Old matrix correpondance
    NewOTUSNames=c(UncollapsedOTUs,names(DescendingTips))
    NewOTUs_OldOTUs_Matrix=Matrix(0,ncol=N_totalNewOTUS,nrow=Ntips,dimnames = list(tree$tip.label,NewOTUSNames))
    
    #Filling  New to Old matrix correpondance / uncollapsed tips
    if (length(UncollapsedOTUs)>1){diag(NewOTUs_OldOTUs_Matrix[UncollapsedOTUs,UncollapsedOTUs])=1}     #just fill the diagonal as tips are unchanged
    if (length(UncollapsedOTUs)==1){NewOTUs_OldOTUs_Matrix[UncollapsedOTUs,UncollapsedOTUs]=1}
    
    #Filling  New to Old matrix correpondance / collapsed tips
    for (i in names(DescendingTips)){NewOTUs_OldOTUs_Matrix[DescendingTips[[i]],i]=1}
    
  } else{print("Only 2 OTUS (or less) present at this resolution -- are you sure this is meaningfull?")} 
  
  
  return(NewOTUs_OldOTUs_Matrix)
}

getNew_OTUs_sample_matrix=function(similarity,sampleOTUs,tree)
{
  CorresMatrix=New_toOld_OTUs_castor(similarity=similarity,tree=tree)
  sample=Matrix(t(sampleOTUs))
  sp=colnames(sample)
  Newsample=t(sample[,sp] %*% CorresMatrix[sp,])
  out=list(NewSample=as.matrix(Newsample),CorrepondanceOTUs=CorresMatrix)
  return(out)
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
    jtu=as.matrix(h[["beta.jtu"]])
    jne=as.matrix(h[["beta.jne"]])
    jac=as.matrix(h[["beta.jac"]])
    stu=as.matrix(hh[["beta.sim"]])
    sne=as.matrix(hh[["beta.sne"]])
    sor=as.matrix(hh[["beta.sor"]])  
    res=abind(jtu,jne,jac,stu,sne,sor,along=0)
    dimnames(res)[[1]]=c("jtu","jne","jac","stu","sne","sor")  
  }
  return(res)
}

#getBDTT=function(similarity,tree,sampleOTUs,onlyBeta=T,metric="jac")
#{
#  New_OTUs_sample_matrix=getNew_OTUs_sample_matrix(similarity=similarity,sampleOTUs=sampleOTUs,tree=tree)
  
  
#  BetasPA=getBeta(t(New_OTUs_sample_matrix[[1]]),ab=F)[c("jtu","jac"),,]
  #rownames(BetasPA)=paste(rownames(BetasPA),"_Simi_",similarity,sep="")
#  BetasAb=getBeta(t(New_OTUs_sample_matrix[[1]]),ab=T)[c("bctu","bc"),,]  
  #rownames(BetasAb)=paste(rownames(BetasAb),"_Simi_",similarity,sep="")
#  AllBetas=abind(BetasPA,BetasAb,along = 1)
  
#  AllBetas=AllBetas[metric,,]
  
#  if (onlyBeta==F){res=list(Beta_Div=AllBetas,NewOTU_table=New_OTUs_sample_matrix[[1]],New_to_old_OTUs=New_OTUs_sample_matrix[[2]])}
#  if (onlyBeta==T){res=AllBetas}
  
#  return(res)
#}


getBDTT=function(similarity,tree,sampleOTUs,onlyBeta=T)
{
  New_OTUs_sample_matrix=getNew_OTUs_sample_matrix(similarity=similarity,sampleOTUs=sampleOTUs,tree=tree)
  OTU_table=t(New_OTUs_sample_matrix[[1]])
  
  Bray=as.matrix(vegdist(OTU_table,method="bray"))
  Jac=as.matrix(vegdist(OTU_table,method="jac",binary = T))
  JacTT=as.matrix(designdist(OTU_table, "(2*pmin(b,c))/(a+2*pmin(b,c))",abcd=TRUE,terms = "binary"))

  AllBetas=abind(Bray,Jac,JacTT,along = 0)
  dimnames(AllBetas)[[1]]=c("Bray","Jac","Jac_TT")
  
  if (onlyBeta==F){res=list(Beta_Div=AllBetas,NewOTU_table=New_OTUs_sample_matrix[[1]],New_to_old_OTUs=New_OTUs_sample_matrix[[2]])}
  if (onlyBeta==T){res=AllBetas}
  
  return(res)
}

BDTT=function(similarity_slices,tree,sampleOTUs,onlyBeta=T){
  
  Betas=lapply(similarity_slices,getBDTT,tree=tree,sampleOTUs=sampleOTUs,onlyBeta=onlyBeta)
  names(Betas)=similarity_slices
  res=do.call(function(...){abind(...,along=0)},Betas)
  return(res)
}


