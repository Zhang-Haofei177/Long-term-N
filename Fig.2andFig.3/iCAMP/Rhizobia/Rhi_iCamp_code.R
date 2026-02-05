wd="~/iCamp/Rhi/Input/"
com.file="~/iCamp/Rhi/Input/Rhi_ASV(zhengshu).txt"
tree.file="~/iCamp/Rhi/Input/rhi_tree.unrooted.nwk"

clas.file <- "~/iCamp/Rhi/Input/Rhizibia_Tax.txt"

treat.file="~/iCamp/Rhi/Input/group.txt"
env.file="~/iCamp/Rhi/Input/env.txt"
save.wd="~/iCamp/Rhi/OutPut/"
if(!dir.exists(save.wd)){dir.create(save.wd)}

prefix="Rhi_1"
rand.time=100
nworker=8
memory.G=50
library(iCAMP)
library(ape)
setwd(wd)
comm=t(read.table(com.file, header = TRUE, sep = "\t", row.names = 1,
                  as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                  check.names = FALSE))
tree=read.tree(file = tree.file)

clas=read.table(clas.file, header = TRUE, sep = "\t", row.names = 1, 
                as.is = TRUE, stringsAsFactors = FALSE, comment.char = "", 
                check.names = FALSE)

treat=read.table(treat.file, header = TRUE, sep = "\t", row.names = 1,
                 as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                 check.names = FALSE)

env=read.table(env.file, header = TRUE, sep = "\t", row.names = 1,
               as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
               check.names = FALSE)

sampid.check=match.name(rn.list=list(comm=comm,treat=treat,env=env))
treat=sampid.check$treat
comm=sampid.check$comm
comm=comm[,colSums(comm)>0,drop=FALSE]
env=sampid.check$env


asv_id <- colnames(comm)
clas <- clas[asv_id,]
spid.check = match.name(cn.list = list(comm = comm), 
                        rn.list = list(clas = clas), 
                        tree.list = list(tree = tree))

comm = spid.check$comm
clas = spid.check$clas
tree = spid.check$tree

setwd(save.wd)
if(!file.exists("pd.desc")) {
  pd.big=iCAMP::pdist.big(tree = tree, wd=save.wd, nworker = nworker, memory.G = memory.G)
}else{
  pd.big=list()
  pd.big$tip.label=read.csv(paste0(save.wd,"/pd.taxon.name.csv"),row.names = 1,stringsAsFactors = FALSE)[,1]
  pd.big$pd.wd=save.wd
  pd.big$pd.file="pd.desc"
  pd.big$pd.name.file="pd.taxon.name.csv"}

setwd(save.wd)
niche.dif=iCAMP::dniche(env = env,comm = comm,method = "niche.value",
                        nworker = nworker,out.dist=FALSE,bigmemo=TRUE,
                        nd.wd=save.wd)

ds = 0.2
bin.size.limit = 24
if(!ape::is.rooted(tree))
{tree.rt = iCAMP::midpoint.root.big(tree = tree, 
                                    pd.desc = pd.big$pd.file, 
                                    pd.spname = pd.big$tip.label, 
                                    pd.wd = pd.big$pd.wd, 
                                    nworker = nworker)
tree = tree.rt$tree
}
phylobin = taxa.binphy.big(tree = tree, 
                           pd.desc = pd.big$pd.file, 
                           pd.spname = pd.big$tip.label, 
                           pd.wd = pd.big$pd.wd, 
                           ds = ds, 
                           bin.size.limit = bin.size.limit, 
                           nworker = nworker)

sp.bin=phylobin$sp.bin[,3,drop=FALSE]
sp.ra=colMeans(comm/rowSums(comm))
abcut=3
commc=comm[,colSums(comm)>=abcut,drop=FALSE]
dim(commc)
spname.use=colnames(commc)
binps=iCAMP::ps.bin(sp.bin = sp.bin,sp.ra = sp.ra,spname.use = spname.use,
                    pd.desc = pd.big$pd.file, pd.spname = pd.big$tip.label, pd.wd = pd.big$pd.wd,
                    nd.list = niche.dif$nd,nd.spname = niche.dif$names,ndbig.wd = niche.dif$nd.wd,
                    cor.method = "pearson",r.cut = 0.1, p.cut = 0.05, min.spn = 5)
if (file.exists(paste0(prefix, ".PhyloSignalSummary.csv"))) {
  appendy <- TRUE
  col.namesy <- FALSE
} else {
  appendy <- FALSE
  col.namesy <- TRUE
}

write.table(data.frame(ds=ds,n.min=bin.size.limit,binps$Index),file = paste0(prefix,".PhyloSignalSummary.csv"),
            append = appendy, quote=FALSE, sep=",", row.names = FALSE,col.names = col.namesy)
if(file.exists(paste0(prefix,".PhyloSignalDetail.csv"))){appendy2=TRUE;col.namesy2=FALSE}else{appendy2=FALSE;col.namesy2=TRUE}
write.table(data.frame(ds=ds,n.min=bin.size.limit,binID=rownames(binps$detail),binps$detail),file = paste0(prefix,".PhyloSignalDetail.csv"),
            append = appendy2, quote = FALSE, sep = ",", row.names = FALSE, col.names = col.namesy2)


ds = 0.2
bin.size.limits = c(12,24, 36, 48)

for (bin.size.limit in bin.size.limits) {
  

  if (!ape::is.rooted(tree)) {
    tree.rt = iCAMP::midpoint.root.big(tree = tree, 
                                       pd.desc = pd.big$pd.file, 
                                       pd.spname = pd.big$tip.label, 
                                       pd.wd = pd.big$pd.wd, 
                                       nworker = nworker)
    tree = tree.rt$tree
  }
  

  phylobin = taxa.binphy.big(tree = tree, 
                             pd.desc = pd.big$pd.file, 
                             pd.spname = pd.big$tip.label, 
                             pd.wd = pd.big$pd.wd, 
                             ds = ds, 
                             bin.size.limit = bin.size.limit, 
                             nworker = nworker)
  

  sp.bin = phylobin$sp.bin[, 3, drop = FALSE]
  sp.ra = colMeans(comm / rowSums(comm))
  abcut = 3  
  commc = comm[, colSums(comm) >= abcut, drop = FALSE]
  dim(commc)
  spname.use = colnames(commc)
  

  binps = iCAMP::ps.bin(sp.bin = sp.bin, sp.ra = sp.ra, spname.use = spname.use,
                        pd.desc = pd.big$pd.file, pd.spname = pd.big$tip.label, pd.wd = pd.big$pd.wd,
                        nd.list = niche.dif$nd, nd.spname = niche.dif$names, ndbig.wd = niche.dif$nd.wd,
                        cor.method = "pearson", r.cut = 0.1, p.cut = 0.05, min.spn = 5)
  

  if (file.exists(paste0(prefix, bin.size.limit, ".PhyloSignalSummary.csv"))) {
    appendy <- TRUE
    col.namesy <- FALSE
  } else {
    appendy <- FALSE
    col.namesy <- TRUE
  }
  

  write.table(data.frame(ds = ds, n.min = bin.size.limit, binps$Index), 
              file = paste0(prefix, bin.size.limit, ".PhyloSignalSummary.csv"),
              append = appendy, quote = FALSE, sep = ",", row.names = FALSE, col.names = col.namesy)
  

  if (file.exists(paste0(prefix, bin.size.limit, ".PhyloSignalDetail.csv"))) {
    appendy2 <- TRUE
    col.namesy2 <- FALSE
  } else {
    appendy2 <- FALSE
    col.namesy2 <- TRUE
  }
  

  write.table(data.frame(ds = ds, n.min = bin.size.limit, binID = rownames(binps$detail), binps$detail), 
              file = paste0(prefix, bin.size.limit, ".PhyloSignalDetail.csv"),
              append = appendy2, quote = FALSE, sep = ",", row.names = FALSE, col.names = col.namesy2)
}

bin.size.limit =24
sig.index = "Confidence" 
icres = iCAMP::icamp.big(comm = comm, 
                         pd.desc = pd.big$pd.file, 
                         pd.spname = pd.big$tip.label, 
                         pd.wd = pd.big$pd.wd, 
                         rand = rand.time, 
                         tree = tree, 
                         prefix = prefix, 
                         ds = 0.2, 
                         pd.cut = NA, 
                         sp.check = TRUE, 
                         phylo.rand.scale = "within.bin", 
                         taxa.rand.scale = "across.all", 
                         phylo.metric = "bMPD", 
                         sig.index = sig.index, 
                         bin.size.limit = bin.size.limit, 
                         nworker = nworker, 
                         rtree.save = FALSE, 
                         detail.save = TRUE, 
                         qp.save = FALSE, 
                         detail.null = FALSE, 
                         ignore.zero = TRUE, 
                         output.wd = save.wd, 
                         correct.special = TRUE, 
                         unit.sum = rowSums(comm), 
                         special.method = "depend", 
                         ses.cut = 1.96, 
                         rc.cut = 0.95, 
                         conf.cut = 0.975, 
                         omit.option = "no", 
                         meta.ab = NULL)


icbin=iCAMP::icamp.bins(icamp.detail = icres$detail,
                        treat = treat, 
                        clas=clas,
                        silent=FALSE, 
                        boot = TRUE, 
                        rand.time = rand.time,
                        between.group = TRUE)
save(icbin,file = paste0(prefix,".iCAMP.Summary.rda")) 
write.csv(icbin$Pt,
          file = paste0(prefix,".ProcessImportance_EachGroup.csv"),
          row.names = FALSE)
write.csv(icbin$Ptk,
          file = paste0(prefix,".ProcessImportance_EachBin_EachGroup.csv"),
          row.names = FALSE)
write.csv(icbin$Ptuv,
          file = paste0(prefix,".ProcessImportance_EachTurnover.csv"),
          row.names = FALSE)
write.csv(icbin$BPtk,
          file = paste0(prefix,".BinContributeToProcess_EachGroup.csv"),
          row.names = FALSE)
write.csv(data.frame(ID=rownames(icbin$Class.Bin),icbin$Class.Bin,stringsAsFactors = FALSE),
          file = paste0(prefix,".Taxon_Bin.csv"),
          row.names = FALSE)
write.csv(icbin$Bin.TopClass,
          file = paste0(prefix,".Bin_TopTaxon.csv"),
          row.names = FALSE)

i=1
treat.use=treat[,i,drop=FALSE]
icamp.result=icres$CbMPDiCBraya
icboot=iCAMP::icamp.boot(icamp.result = icamp.result,
                         treat = treat.use,
                         rand.time = rand.time, 
                         compare = TRUE,
                         silent = FALSE,
                         between.group = TRUE,
                         ST.estimation = TRUE)
save(icboot,file=paste0(prefix,".iCAMP.Boot.",colnames(treat)[i],".rda"))
write.csv(icboot$summary,file = paste0(prefix,".iCAMP.BootSummary.",colnames(treat)[i],".csv"),row.names = FALSE)
write.csv(icboot$compare,file = paste0(prefix,".iCAMP.Compare.",colnames(treat)[i],".csv"),row.names = FALSE)

df_tax<-icres$detail$taxabin$sp.bin
bin_number<-unique(df_tax$bin.id.new)
write.csv(icres$detail$taxabin$sp.bin,file = paste0(prefix,".iCamp.Bin.taxon",".csv"))
