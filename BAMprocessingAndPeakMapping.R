# ______________________________________________** * **______________________________________________ # 
#LIBRARIES AND CUSTOM FUNCTIONS

# ______________________________________________** * **______________________________________________ # 
library(seqinr)
library(RiboProfiling)
library(GenomicFeatures)
library(GenomicAlignments)
library(GenomicRanges)
library(foreach)
library(doParallel)
library(data.table)
library(scales)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(eif)
library(ggpubr)
library(ggseqlogo)
library(ape)
library(reshape)
library(stringdist)


loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

getCodonCoverage <- function(alignment, targetLengths, txdb, seqs, coresToUse)
{
  #parse exons from txdb
  cds <- cdsBy(txdb, by="tx", use.names=TRUE);
  exonGRanges <- exonsBy(txdb, by="tx", use.names=TRUE);
  cdsPosTransc <- orfRelativePos(cds, exonGRanges);
  
  #read aln to start
  grangeAln <- readsToStartOrEnd(aln = alignment, what = "start");#put as end just to process disomeseq from the 3' en monosome of the diesome,

  #filter by read length
  grangeAln <- grangeAln[which(!is.na(match(grangeAln$score, targetLengths)))];
  readLens <- t(as.matrix(table(qwidth(alignment)))); 
  barplot(readLens, xlab = "read length", ylab = "number of reads",
          main = "Read length distribution", names.arg = colnames(readLens),
          col = "white");
  readLens[which(!(as.numeric(colnames(readLens)) %in% targetLengths))] <- 0;
  barplot(readLens, xlab = "read length", ylab = "number of reads",
          main = "Read length distribution", names.arg = NULL,
          col = "firebrick4", border = NA, add = T);
  legend("topright", legend = "Selected lengths", fill = "firebrick4", cex = 1.5, bty = "n");
  
  ##to correct riboprofiling offset
  # listeRangesCDS <- lapply(seq_len(NROW(cdsPosTranscShifted)), 
  #                          function(ixTransc) {
  #                            if (shiftValue >= 0) {
  #                              max(1, cdsPosTranscShifted[ixTransc, 1]):min(newTranscWidth[ixTransc], 
  #                                                                           cdsPosTranscShifted[ixTransc, 2])
  #                            }
  #                            else {
  #                              c(min(1, cdsPosTranscShifted[ixTransc, 1] - 1):-1,1:min(newTranscWidth[ixTransc], 
  #                                                                                      cdsPosTranscShifted[ixTransc, 2]))
  #                            }
  #                          })
  
  #get individual offsets for each read length
  oneBinRanges <- aroundPromoter(txdb, grangeAln, percBestExpressed=0.05)
  listPromoterCov <-
    readStartCov(
      grangeAln,
      oneBinRanges,
      matchSize = targetLengths,
      fixedInterval = c(-20, 20),
      renameChr = "aroundTSS",
      charPerc = "perc"
    );
  print(plotSummarizedCov(listPromoterCov));
  
  # listPromoterCov <-
  #   readStartCov(
  #     grangeAln,
  #     oneBinRanges,
  #     matchSize=targetLengths,
  #     fixedInterval=c(-20, 20),
  #     renameChr="aroundTSS",
  #     charPerc="perc"
  #   )
  # head(grangeAln@elementMetadata@nrows)
  
  
  offsets <- rep(0, length(listPromoterCov));
  for(i in 1:length(listPromoterCov))
    offsets[i] = which.max(listPromoterCov[[i]]@elementMetadata@listData$values) - 21;
  
  #substract offsets and count reads (individual offsets, not a unique one)
  countedReads <- list();
  print("Estimated offsets: ");
  print(data.frame("length" = targetLengths, "offset" = offsets[-1]));
  print(paste("Overall offset: ", offsets[1], sep = ""));
  offsets <- offsets[-1]; #first offset comes from the summary of all used lengths
  lensToUse <- readline(prompt = "Keep all read lengths? (y/n) ");
  if(lensToUse == "n")
  {
    lensToUse <- readline(prompt = "Enter read lengths to keep: (e.g. 26,27,28) ");
    lensToUse <- as.numeric(strsplit(lensToUse, split = ",")[[1]]);
    offsets <- offsets[which(targetLengths %in% lensToUse)];
    targetLengths <- lensToUse;
  }
  offsetToUse <- readline(prompt = "Use estimated offsets? (y/n) ");
  if(offsetToUse == "n")
  {
    for(l in 1:length(targetLengths))
      offsets[l] <- as.integer(readline(prompt = paste("offset length: ", targetLengths[l], sep = "") ));
  }
  print(data.frame("length" = targetLengths, "offset" = offsets));
  uoffsets <- unique(offsets);
  for(i in 1:length(uoffsets))
  {
    tempgrangeAln <- grangeAln[which(!is.na(match(grangeAln$score, targetLengths[which(offsets == uoffsets[i])])))];
    
    ###filter by frame
    # temp <- as(aln,"GRanges")
    # dt <- mapToTranscripts(temp, transcripts = transcripts(txdbScer, use.names = T))
    # frames <- (dt@ranges@start -11)%%3
    # dt <- mapFromTranscripts(dt, transcripts = transcripts(txdbScer, use.names = T))
    # mcols(dt, level="within")$frame <- frames;
    # alignment <- as(dt[dt$frame == 0], "GAlignments")
    
    print(as.integer(uoffsets[i]));
    #count reads applying estimated offset (P-site centering)
    countedReads[[i]] <-
      countShiftReads(
        exonGRanges = exonGRanges[names(cdsPosTransc)],
        cdsPosTransc = cdsPosTransc,
        alnGRanges = tempgrangeAln,
        shiftValue = uoffsets[i],
        motifSize = as.integer(3)
      );
  }
  
  #parse output
  uoffsets <- length(uoffsets);
  ids <- as.character(seqs@elementMetadata$tx_name);
  uids <- unique(ids);
  cores <- makeCluster(coresToUse);
  registerDoParallel(cores);
  hora <- Sys.time();
  readsPerORF <- foreach( i = 1:length(uids), .packages = c("Biostrings","BSgenome", "seqinr") ) %dopar% {
    id <- uids[i];
    countOut <- data.frame("codonID" = NA, "nbrReads" = NA)
    for(j in 1:uoffsets)
    {
      if(! length(countedReads[[j]][[2]][[id]]))
        next;
      countOut <- merge(countOut, countedReads[[j]][[2]][[id]], by = "codonID", all = T);
    }
    if(length(countOut$codonID) > 1)
    {
      slen <- 0;
      j <- which(ids == id);
      jlen <- length(j);
      if(jlen == 1)
      {
        slen <- DNAString(seqs[[j]])@length;
        dnaSeq <-  as.character(DNAString(seqs[[j]]));
      } else if(jlen > 1)
      {
        exs <- exonGRanges[[ids[j[1]]]];
        exs <- exs[order(exs@elementMetadata$exon_rank, decreasing = F)]@elementMetadata$exon_id;
        tempseqs <- seqs[j];
        dnaSeq <- c();
        for(j in exs)
        {
          dnaSeq <- paste(dnaSeq, as.character(DNAString(tempseqs[[which(tempseqs@elementMetadata$exon_id == j)]])),
                          sep = "");
          slen <- slen + DNAString(tempseqs[[which(tempseqs@elementMetadata$exon_id == j)]])@length;
        }
      }else
        next;
      if( slen %% 3 )
      {
        print(c(i, id, slen))
        print(c(dnaSeq[1], dnaSeq[length(dnaSeq)]))
      }
      countOut <- countOut[-which(is.na(countOut$codonID)),];
      countOut$total <- rowSums(countOut[,2:length(countOut)], na.rm = T);
      temp <- data.frame("codon" = splitseq(s2c(dnaSeq), frame = 0, word = 3), stringsAsFactors = F);
      temp$codonID <- as.numeric(1:length(temp$codon));
      temp <- merge(countOut[,c(1,length(countOut))], temp, by = "codonID", all = T);
      temp$total[is.na(temp$total)] <- 0;
      colnames(temp) <- c("number", "reads", "codon");
      temp$normMean <- temp$reads/sum(temp$reads);
      to.readsPerORF <-  list("counts" = temp[c(1,3,2,4)], "orf" = id);
    }
  }
  stopCluster(cores);
  print(Sys.time() - hora);
  
  creads <- list();
  totalreads <- data.frame("orf" = uids, "reads" = rep(0, length(uids)), "length" = rep(0, length(uids)),
                           stringsAsFactors = F);
  for(i in 1:length(readsPerORF))
  {
    if(length(readsPerORF[[i]]))
    {
      j <- which(totalreads$orf == readsPerORF[[i]]$orf);
      creads[[readsPerORF[[i]]$orf]] <- readsPerORF[[i]]$counts;
      totalreads$reads[j] <- sum(readsPerORF[[i]]$counts$reads);
      totalreads$length[j] <- length(readsPerORF[[i]]$counts$number)/1000;
    }
  }
  totalreads$tpm <- totalreads$reads/totalreads$length;
  totalreads$tpm <- totalreads$tpm/(sum(totalreads$tpm, na.rm = T)/1000000);
  rownames(totalreads) <- totalreads$orf;
  ids <- names(creads);
  for(i in 1:length(creads))
  {
    creads[[i]]$normTPM <- creads[[i]]$reads/totalreads[ids[i], 4];
  }
  
  #return estimated offsets and counted reads
  return( list("offsets" = offsets, "counts" = creads, "summary" = totalreads) );
}

outlierByIsoForest <- function(processedReads, alpha, ntrees, samplesize)
{
  temp <- processedReads$summary$tpm;
  
  layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8))
  par(mar=c(0, 3.1, 1.1, 2.1))
  boxplot(temp , horizontal=TRUE , xaxt="n" , ylim = c(0,150),
          col = alpha("grey", 0.4), frame=F, pch = 16, border = "black")
  par(mar=c(4, 3.1, 1.1, 2.1))
  plot(density(temp[which(temp < 500)]), col = "firebrick4", lwd = 2,
       xlab = "TPM", ylim = c(0,0.03), xlim = c(0,150), bty = "n", main = NA)
  explims <- boxplot.stats(temp)$stats;
  abline(v = explims, lty = 2,col = "grey", lwd = 2);
  explims <- c(explims, 1e+06);
  
  for(i in 1:length(processedReads$counts))
    processedReads$counts[[i]]$normNumber <- processedReads$counts[[i]]$number/max(processedReads$counts[[i]]$number);
  orflengths <- vapply(processedReads$counts, nrow, 1L);
  orfnames <- names(processedReads$counts);
  dt2 <- rbindlist(processedReads$counts)[, orf := rep(orfnames, times = orflengths)][];
  dt2[,score := 0];
  dt2[,peak := 0];
  
  # for(orf in unique(merged_datasets$ORF))
  # {
  #   merged_datasets[ORF == orf, normNumber := position/(max(position)+3)];
  # }
  # dt2 <- merged_datasets[,.(ORF, normNumber, orfnorm.wtr1)];
  # colnames(dt2) <- c("orf", "normNumber", "normTPM")
  
  for(i in 1:(length(explims)- 1)){
    orfs <- processedReads$summary$orf[which(processedReads$summary$tpm > explims[i] & 
                                               processedReads$summary$tpm <= explims[i + 1] )];
    positions <- dt2[orf %in% orfs, which = T];
    iForest <- eif(X = as.matrix(dt2[positions, .(normNumber, normTPM)]), ntrees = as.integer(ntrees),
                   sample_size = as.integer(samplesize), ExtensionLevel = as.integer(1));
    dt2[positions, score := iForest$scores];
    dt2[positions, quantile := as.numeric(quantile(iForest$score, alpha)[[1]])];
  }
  dt2[, peak := 0];
  dt2[!is.na(score) & score > quantile, peak := 1];
  
  # dt2[,position := merged_datasets$position]
  # dt2[,id := paste(orf, position,sep="_")]
  # dt2[,peakj := peak]
  # wtdt <- merge(wtdt, dt2[,c("id", "peakj")], by = "id", all = F)
  
  return(dt2);
}

addScerDescriptors <- function(dt, sss, codonfeatures)
{
  aafeatures <- rbind(read.csv("~/Documents/tricNetwork/aaHydrophobicityAndVolume.csv",
                               header = T, stringsAsFactors = F), c("*", NA, NA, NA, NA, NA));
  dt[, first := strsplit(orf, split = "")];
  dt[,first := sapply(dt$first, function(x)x[[1]])];
  for(i in 1:length(codonfeatures$codon))
  {
    dt[codon == codonfeatures$codon[i] & first != "Q", aa := codonfeatures$aa[i]];
    dt[codon == codonfeatures$codon[i] & first == "Q", aa := codonfeatures$aam[i]];
    dt[codon == codonfeatures$codon[i], tai := codonfeatures$tai[i]];
    dt[codon == codonfeatures$codon[i], optimality := codonfeatures$optimality[i]];
    dt[codon == codonfeatures$codon[i], cu := codonfeatures$cu[i]];
    dt[codon == codonfeatures$codon[i], trnaCN := codonfeatures$trnaCN[i]];
  }
  for(i in 1:length(aafeatures$residue))
  {
    dt[aa == aafeatures$residue[i], hydrophobicity := aafeatures$hydrophobicity[i]];
    dt[aa == aafeatures$residue[i], volume := aafeatures$volume[i]];
  }
  orfs <- unique(sss[,orf]);
  orfs <- orfs[which(orfs %in% unique(dt[,orf]))];
  for(i in orfs)
  {
    positions <- dt[orf == i & aa != "*", which = T];
    positions2 <- sss[orf == i, which = T];
    if(length(positions) == length(positions2))
    {
      if(paste(dt[positions, aa], collapse = "") == paste(sss[positions2, aa], collapse = ""))
      {
        dt[positions, ss := sss[positions2, ss]];
        dt[positions, solubility := sss[positions2, solubility]];
      }else{
        print(paste("Different sequences for: ", i, sep = ""));
      }
    }else{
      print(paste("Different residue number for: ", i, sep = ""));
    }
  }
}

mergeReplicas <- function(wt1, wt2)
{
  wt1 <- wt1[,c(1,6,2,7,5,8,9)];
  colnames(wt1) <- c("number", "normNumber", "codon",
                     "orf", "ntpm1", "iscore1", "peak1");
  wt1[, id := as.character(paste(orf, number, sep = "_"))];
  wt2 <- wt2[,c(1,6,2,7,5,8,9)];
  colnames(wt2) <- c("number", "normNumber", "codon",
                     "orf", "ntpm2", "iscore2", "peak2");
  wt2[, id := as.character(paste(orf, number, sep = "_"))];
  wtdt <- merge(wt1, wt2[, 5:8], all = F, by = "id", sort = F);
  
  addScerDescriptors(wtdt, solandSS, codonfeatures);
  for(i in unique(wtdt[peak1 > 0 | peak2 > 0, orf]))
  {
    j <- wtdt[orf == i & peak1 > 0 & ntpm1 < mean(wtdt[orf == i, ntpm1]), which = T];
    if(length(j))
      wtdt[j, peak1 := -1];
    j <- wtdt[orf == i & peak2 > 0 & ntpm2 < mean(wtdt[orf == i, ntpm2]), which = T];
    if(length(j))
      wtdt[j, peak2 := -1];
    # print(i);
  }
  return(wtdt);
}

# ______________________________________________** * **______________________________________________ # 
#DATA PROCESSING -FROM BAM FILES-

setwd("~/Documents/RibosomeStalling/");

#for yeast, load the following:
txdbScer <- makeTxDbFromGFF("~/Downloads/ScerGTF.gtf");
exomeScer <- exons(txdbScer, columns=c("exon_id", "tx_name", "gene_id"));
seqsScer <- Views(BSgenome.Scerevisiae.UCSC.sacCer3, exomeScer);

fragLens <- list();
bamfiles <- list.files(path = "MappedReads/", pattern = ".bam$")
#use RiboProfiling to create lists of counts per sample/orf and save them individually
for(i in bamfiles)
{
  id <- as.character(strsplit(i, split = ".bam")[[1]]);
  alignment <- readGAlignments(paste("MappedReads/", i, sep = ""));
  
  readLens <- qwidth(alignment);
  fragLens[[id]] <- as.data.frame(table(readLens));
  barplot(t(as.matrix(fragLens[[id]])), xlab = "read length", ylab = "number of reads",
          main = id, names.arg = fragLens[[id]]$readLens);
  lowerlim <- as.integer(readline(prompt = "Enter large fragments lower limit: "));
  upperlim <- as.integer(readline(prompt = "Enter large fragments upper limit: "));
  
  countOut <- getCodonCoverage(alignment, lowerlim:upperlim, txdbScer, seqsScer, 10);
  do.call("<-", list(id, countOut));
  save(list = paste(id, sep=""), file = paste(id, "_counts2.RData", sep = ""));
}

# ______________________________________________** * **______________________________________________ # 
#OUTLIER IDENTIFICATION USING EXTENDED ISOLATION FORESTS/ADD DESCRIPTORS
load("~/Documents/tricNetwork/s2DandCamSol_Scerevisiae.RData");
load("~/Documents/tricNetwork/codonFeatures.RData");

#join replicates into a single data table (not actual data merging, but joining)
mergeReplicas <- function(wt1, wt2)
{
  wt1 <- wt1[,c(1,6,2,7,5,8,9)];
  colnames(wt1) <- c("number", "normNumber", "codon",
                     "orf", "ntpm1", "iscore1", "peak1");
  wt1[, id := as.character(paste(orf, number, sep = "_"))];
  wt2 <- wt2[,c(1,6,2,7,5,8,9)];
  colnames(wt2) <- c("number", "normNumber", "codon",
                     "orf", "ntpm2", "iscore2", "peak2");
  wt2[, id := as.character(paste(orf, number, sep = "_"))];
  wtdt <- merge(wt1, wt2[, 5:8], all = F, by = "id", sort = F);
  
  addScerDescriptors(wtdt, solandSS, codonfeatures);
  for(i in unique(wtdt[peak1 > 0 | peak2 > 0, orf]))
  {
    j <- wtdt[orf == i & peak1 > 0 & ntpm1 < mean(wtdt[orf == i, ntpm1]), which = T];
    if(length(j))
      wtdt[j, peak1 := -1];
    j <- wtdt[orf == i & peak2 > 0 & ntpm2 < mean(wtdt[orf == i, ntpm2]), which = T];
    if(length(j))
      wtdt[j, peak2 := -1];
    # print(i);
  }
  return(wtdt);
}

dt <- loadRData("processedCounts/WT1_counts.RData");
wt1 <- outlierByIsoForest(dt, 0.95, 200, 300);
dt <- loadRData("processedCounts/WT2_counts.RData");
wt2 <- outlierByIsoForest(dt, 0.95, 200, 300);
save(wt1, file = "WT1_iForestStalls_95q.RData");
save(wt2, file = "WT2_iForestStalls_95q.RData");
wt12 <- mergeReplicas(day21, day22);
save(wt12, file = "WT12_iForestStallsWDescriptors_WPauseScores_NoTrimming.RData");
# ______________________________________________** * **______________________________________________ # 

