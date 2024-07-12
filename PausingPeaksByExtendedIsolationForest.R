# ##GENERAL CONSIDERATIONS
# 1. There are two options to run this script. The first one starts from an annotated BAM file
# and in its current form works only for S. cerevisiae. You can edit and modify the needed
# txdb files and the getCodonCoverage function to adapt it to other organisms/translatomes.
# It should be fairly straightforward. The second option -recommended- requires a csv file
# as input, which contains the per-codon raw read counts. The format of the table is described
# below and the zenodo directory contains a full example for S. cerevisiae data.
# 
# 2. Parameters of the outlierByIsoForest function:
#   -> processedReads: a data.table containing the read counts, of the format:
#         number codon reads      orf
#              1   ATG   345  YAL069W
#              2   ATC     0  YAL069W
#              3   GTA  1944  YAL069W
#             .
#             .
#             .
#   -> alpha: quantile threshold for outlier classification. 0.95 recommended.
#   -> ntrees: number of trees for the ensemble. Recommended 200 or greater.
#   -> samplesize: sub-sampling size for training of each tree. Recommended 300 or greater.
#   -> limits: the mode of estimation of the TPM intervals. Each estimated interval
#             results in an independent EIF model (see Methods of manuscript). Values
#             are either "static" or "dynamic": static results in TPM intervals coming from the
#             0,0.2,0.4,0.6,0.8,1 quantiles, as in the manuscript; dynamic allocates the limits
#             based on the boxplot.stats function of the grDevices library.
# 
# 3. The outlierByIsoForest function returns a data.table, with the same columns as the
# input data.table plus the additional following ones:
#     ->normTPM: reads normalized by the ORF's TPM value
#     ->normNumber: per ORF, normalized codon position (0,1]
#     ->score: the EIF score for that position
#     ->peak: boolean, 0=no pausing peak; 1=pausing peak
#     ->quantile: estimated quantile threshold used (each TPM interval would have its own)
#     


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
library(grDevices)


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

outlierByIsoForest <- function(processedReads, alpha, ntrees, samplesize, limits)
{
  totalreads <- merge(processedReads[, sum(reads), by = "orf"],
                      processedReads[, max(number), by = "orf"], all = T, order = F, by = "orf");
  colnames(totalreads) <- c("orf", "reads", "length");
  totalreads[, tpm := reads/length];
  x <- sum(totalreads[,tpm], na.rm = T);
  totalreads[, tpm := (tpm/x)*1e+06];
  
  processedReads[, normTPM := 0];
  processedReads[, normNumber := number/max(number), by = "orf"]
  for(i in 1:totalreads[,.N])
    processedReads[orf == totalreads[i,orf], normTPM := reads/totalreads[i, tpm]];
  layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8))
  par(mar=c(0, 3.1, 1.1, 2.1))
  boxplot(totalreads[,tpm] , horizontal=TRUE , xaxt="n" , ylim = c(0,150),
          col = alpha("grey", 0.4), frame=F, pch = 16, border = "black")
  par(mar=c(4, 3.1, 1.1, 2.1))
  plot(density(totalreads[tpm < 500,tpm]), col = "firebrick4", lwd = 2,
       xlab = "TPM", ylim = c(0,0.03), xlim = c(0,150), bty = "n", main = NA)
  
  if(limits == "dynamic"){
    explims <- boxplot.stats(totalreads[,tpm])$stats;
    explims <- c(explims, 1e+06);
  }else{
    explims <- as.numeric(quantile(totalreads[,tpm],c(0.2,0.4,0.6,0.8)));
    explims <- c(0, explims, 1e+06)
  }
  abline(v = explims, lty = 2,col = "grey", lwd = 2);
  
  processedReads[,score := 0];
  processedReads[,peak := 0];
  
  # for(orf in unique(merged_datasets$ORF))
  # {
  #   merged_datasets[ORF == orf, normNumber := position/(max(position)+3)];
  # }
  # dt2 <- merged_datasets[,.(ORF, normNumber, orfnorm.wtr1)];
  # colnames(dt2) <- c("orf", "normNumber", "normTPM")
  
  for(i in 1:(length(explims)- 1)){
    orfs <- totalreads[tpm > explims[i] & tpm <= explims[i + 1], orf];
    
    positions <- processedReads[orf %in% orfs, which = T];
    iForest <- eif(X = as.matrix(processedReads[positions, .(normNumber, normTPM)]), ntrees = as.integer(ntrees),
                   sample_size = as.integer(samplesize), ExtensionLevel = as.integer(1));
    processedReads[positions, score := iForest$scores];
    processedReads[positions, quantile := as.numeric(quantile(iForest$score, alpha)[[1]])];
  }
  processedReads[, peak := 0];
  processedReads[!is.na(score) & score > quantile, peak := 1];
  
  # dt2[,position := merged_datasets$position]
  # dt2[,id := paste(orf, position,sep="_")]
  # dt2[,peakj := peak]
  # wtdt <- merge(wtdt, dt2[,c("id", "peakj")], by = "id", all = F)
  
  return(processedReads);
}

# ______________________________________________** * **______________________________________________ # 
#OUTLIER IDENTIFICATION USING EXTENDED ISOLATION FORESTS

# ______________________________________________** * **______________________________________________ # 

##RUN OPTION NUMBER 1 - TWO STEPS
##STARTING FROM BAM FILE - S.CEREVISIAE ONLY -

# 1) BAM file processing and annotation

txdbScer <- makeTxDbFromGFF("~/Downloads/ScerGTF.gtf"); #dowload from the zenodo repository of this paper
exomeScer <- exons(txdbScer, columns=c("exon_id", "tx_name", "gene_id"));
seqsScer <- Views(BSgenome.Scerevisiae.UCSC.sacCer3, exomeScer);

fragLens <- list();
bamfiles <- list.files(path = "mappedReads/", pattern = ".bam$"); #modify path accordignly
for(i in bamfiles)
{
  id <- as.character(strsplit(i, split = ".bam")[[1]]);
  alignment <- readGAlignments(paste("mappedReads/", i, sep = "")); #modify path accordignly
  
  readLens <- qwidth(alignment);
  fragLens[[id]] <- as.data.frame(table(readLens));
  barplot(t(as.matrix(fragLens[[id]])), xlab = "read length", ylab = "number of reads",
          main = id, names.arg = fragLens[[id]]$readLens);
  lowerlim <- as.integer(readline(prompt = "Enter fragments lower limit: "));
  upperlim <- as.integer(readline(prompt = "Enter fragments upper limit: "));
  
  countOut <- getCodonCoverage(alignment, lowerlim:upperlim, txdbScer, seqsScer, 10);
  do.call("<-", list(id, countOut));
  save(list = paste(id, sep=""), file = paste(id, "_countsTestScript.RData", sep = ""));
}

# 2) Outlier by Extended Isolation Forest

dt <- loadRData("WT1_countsTestScript.RData"); #A file from the previous step, modify path accordignly

orflengths <- vapply(dt$counts, nrow, 1L);
orfnames <- names(dt$counts);
dt <- rbindlist(dt$counts)[, orf := rep(orfnames, times = orflengths)][];

wt1 <- outlierByIsoForest(dt[,c(1:3,6)], 0.95, 200, 300, "static");#returns a data.table

# ______________________________________________** * **______________________________________________ # 

##RUN OPTION NUMBER 2 - ONE STEP, RECOMMENDED OPTION
##STARTING WITH A CSV TABLE - MAKE SURE IT FOLLOWS THE FORMAT:

# number codon reads      orf
#      1   ATG   345  YAL069W
#      2   ATC     0  YAL069W
#      3   GTA  1944  YAL069W
#     .
#     .
#     .
#    459   TGA   341  YAL069W
#      1   ATG    45  YAR014C
#      2   GGC   129  YAR014C
#     .
#     .
#     .

#WHERE:
    # number: codon position in orf - numeric
    # codon: codon string - character
    # reads: read counts -NOT NORMALIZED, RAW- for a given codon - numeric
    # orf: name of the ORF

#FOR EACH ORF, COLUMN NUMBER SHOULD GO FROM 1 (ATG) TO N-TH CODON (STOP CODON)
#AND THE ORF COLUMN SHOULD REPEAT THE SAME ORF NAME.
#ALL ORFS SHOULD BE INCLUDED, OTHERWISE TPM ESTIMATION WILL BE INACCURATE

#the example file -testScerevisiae.csv- can be found in the same zenodo directory as this file
dt <- as.data.table(read.csv("testScerevisiae.csv", header = T,
                              stringsAsFactors = F)); #csv file, modify if other format
wt1 <- outlierByIsoForest(dt, 0.95, 200, 300, "static");#returns a data.table

