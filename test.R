#cluster.dir <- "gs://gipseq-clusters/pools/36bp-hg38/GC039"
#cluster.name <- "GC039790"
#sample.name <- "GC039790-AR012"
examples.FromReadStarts <- function(cluster.dir,cluster.name, sample.name,partition) {
    partition.read.starts <- ReadReadStarts(cluster.dir,cluster.name,sample.name)
    print("Read starts read")
    partition.bin.counts <- GetRawPartitionBinCounts(partition,partition.read.starts)
    print("Bin counts calculated")
    partition.bin.counts <- AddGCcorrectedPartitionBinCounts(partition,partition.bin.counts)
    partition.bin.counts.raw <- AddPartitionBinCounts(partition,partition.bin.counts,GCC=FALSE)
    partition.bin.counts.gcc <- AddPartitionBinCounts(partition,partition.bin.counts,GCC=TRUE)
    WritePartitionBinCounts(partition.bin.counts.raw,cluster.dir,cluster.name,sample.name, GCC=FALSE)
    WritePartitionBinCounts(partition.bin.counts.gcc,cluster.dir,cluster.name,sample.name, GCC=TRUE)

    for (GCC in c(TRUE,FALSE)) {
        partition.bin.counts <- ReadPartitionBinCounts(cluster.dir,cluster.name,sample.name, GCC)

        for (window.size in c(500,1000)) {
	    print(paste("Processing GCC=",GCC,"window.size=",window.size))
            partition.gipseq.counts <- GetSamplePartitionGipSeqCounts(partition.bin.counts, c("mean","sd"), partition, window.size)
            WritePartitionGipSeqCounts(partition.gipseq.counts,cluster.dir,cluster.name,sample.name,window.size, GCC)
        }
    }
}
    
examples.GipSeqMetrics <- function() {
    partition <- LoadPartition(partition.file)
    roi <- LoadRegionsOfInterest(roi.file)
    chr.size <- LoadChrSize(chr.size.file)
    all.roi <- GetAllROI()
    cluster.name <- "work"
    cluster.dir <- "/uz/data/shortcuts/genomicscore/folders/raw/HiSeqComputed/apps/cme_svn/gc/apps/gipseq"
    sample.name <- gsub(".bam","",basename(my.bamfile))

    partition.read.starts <- GetPartitionReadStarts(my.samtools, my.bamfile, partition, sample.name, verbose=TRUE)
      #   CHR START NM    BIN         SAMPLE
      #1 chr1 10151  0 chr1.1 GC058001-AR001
      #2 chr1 16079  0 chr1.1 GC058001-AR001
    WriteReadStarts(partition.read.starts,cluster.dir,cluster.name, sample.name)

    partition.bin.counts <- GetRawPartitionBinCounts(partition,partition.read.starts)
      #      BINDEX RAW.COUNT         SAMPLE
      #1          1         0 GC058001-AR001
      #2          2         2 GC058001-AR001
    partition.bin.counts <- AddGCcorrectedPartitionBinCounts(partition,partition.bin.counts)
      #  BINDEX RAW.COUNT GCC.COUNT         SAMPLE
      #1      1         0 0.0000000 GC058001-AR001
      #2      2         2 3.3051143 GC058001-AR001
    partition.bin.counts.raw <- AddPartitionBinCounts(partition,partition.bin.counts,GCC=FALSE)
      #  BINDEX SAMPLE COUNT
      #1      1         GC058001-AR001     0.00000000
      #2      2         GC058001-AR001     0.05882353
    partition.bin.counts.gcc <- AddPartitionBinCounts(partition,partition.bin.counts,GCC=TRUE)
      #  BINDEX   SAMPLE COUNT
      #1      1         GC058001-AR001     0.00000000
      #2      2         GC058001-AR001     0.10021460
    WritePartitionBinCounts(partition.bin.counts.raw,cluster.dir,cluster.name,sample.name, GCC=FALSE)
    WritePartitionBinCounts(partition.bin.counts.gcc,cluster.dir,cluster.name,sample.name, GCC=TRUE)

    partition.bin.mutation.burden <- GetPartitionBinMutationBurden(partition, partition.read.starts, read.len=READ.LENGTH)
      #      BINDEX     BIN.MUT         SAMPLE
      #1          1 0.000000000 GC058001-AR001
      #2          2 0.000000000 GC058001-AR001
    WritePartitionBinMutationBurden(partition.bin.mutation.burden,cluster.dir,cluster.name,sample.name)

    normalised.roi.bin.counts <- GetRoiBinCounts(partition.read.starts, all.roi)
      #   CHR START      END                 DESCRIPTION COUNT      RPKM
      #1 chr1     1  2000000                       1pter  2059 0.1075401
      #2 chr1     1 27600000 1p36 microdeletion syndrome 72555 0.2746011
    WriteRoiBinCounts(normalised.roi.bin.counts,cluster.dir,cluster.name,sample.name)

    roi.bin.mutation.burden <- GetRoiBinMutationBurden(partition.read.starts, all.roi, read.len=READ.LENGTH)
      #   CHR  START      END                 DESCRIPTION     BIN.MUT
      # 1 chr1     1  2000000                       1pter 0.001268145
      # 2 chr1     1 27600000 1p36 microdeletion syndrome 0.001096486
    WriteRoiBinMutationBurden(roi.bin.mutation.burden,cluster.dir,cluster.name,sample.name)

    partition.gipseq.counts.raw <- GetSamplePartitionGipSeqCounts(partition.bin.counts.raw, c("mean","sd"), partition, NR.BINS.IN.WINDOW)
    partition.gipseq.counts.gcc <- GetSamplePartitionGipSeqCounts(partition.bin.counts.gcc, c("mean","sd"), partition, NR.BINS.IN.WINDOW)
      #  BINDEX  CHR       MEAN         SD
      #1      1 chr1 0.02486406 0.04976587
      #2      2 chr1 0.02438591 0.04939604
    WritePartitionGipSeqCounts(partition.gipseq.counts.raw,cluster.dir,cluster.name,sample.name,NR.BINS.IN.WINDOW,GCC=FALSE)
    WritePartitionGipSeqCounts(partition.gipseq.counts.gcc,cluster.dir,cluster.name,sample.name,NR.BINS.IN.WINDOW,GCC=TRUE)

    partition.chr.counts.raw <- GetSamplePartitionChrCounts(partition.bin.counts.raw, c("sum"), partition)
    partition.chr.counts.gcc <- GetSamplePartitionChrCounts(partition.bin.counts.gcc, c("sum"), partition)
    WritePartitionChrCounts(partition.chr.counts.raw,cluster.dir,cluster.name,sample.name, GCC=FALSE, blacklist.code=blacklist.regions.code)
    WritePartitionChrCounts(partition.chr.counts.gcc,cluster.dir,cluster.name,sample.name, GCC=TRUE, blacklist.code=blacklist.regions.code)

    partition.gipseq.burden <- GetSamplePartitionGipSeqMutationBurden(partition.bin.mutation.burden, c("mean","sd"), partition, NR.BINS.IN.WINDOW)
    WritePartitionGipSeqMutationBurden(partition.gipseq.burden,cluster.dir,cluster.name,sample.name,NR.BINS.IN.WINDOW)
    partition.bin.mutation.burden <- ReadPartitionBinMutationBurden(cluster.dir,cluster.name,sample.name)


    GCC <- FALSE
    partition.bin.counts.raw <- ReadPartitionBinCounts(cluster.dir,cluster.name,sample.name, GCC)
    for (window.size in c(500,1000)) {
	    print(paste(window.size,GCC))
            partition.gipseq.counts <- GetSamplePartitionGipSeqCounts(partition.bin.counts, c("mean","sd"), partition, window.size)
            WritePartitionGipSeqCounts(partition.gipseq.counts,cluster.dir,cluster.name,sample.name,window.size,GCC)
    }

    GCC <- FALSE
    for (window.size in c(100,500,1000)) {
	 cpgc <- GetClusterPartitionGipSeqCounts(cluster.dir, cluster.name, partition, window.size, GCC)
	 WritePartitionGipSeqCounts(cpgc,cluster.dir,cluster.name,cluster.name,window.size,GCC)
    }

}

##################
# PARTITION BINS #
##################

GetSamplePartitionGipSeqMutationBurden <- function(partition.bin.mutation.burden, agg.functions, partition, window.size) {
  #' @title Get GipSeq burden of a sample
  #' @description Calculate aggregated partition bin burden over a fixed-sized sliding window
  #' @param partition.bin.burden data.frame with burden per bin; columns: "BINDEX (SAMPLE)" and column "BIN.MUT" 
  #' @param agg.functions The vector of aggregation functions to be applied, should be subset of: c("median","mad","mean","sd")
  #' @param partition Partition bin name and GC content data frame as obtained from LoadPartition
  #' @param  window.size Number of subsequent bins over which to aggregate
  #' @return  GipSeq-burden as data.frame with columns "BINDEX CHR" and requested columns from "MEDIAN MAD MEAN SD", sorted on BINDEX 

    partition.bin.mutation.burden <- EnsureBindexHasChr(partition.bin.mutation.burden,partition)
    gipseq.burden <- partition[c("BINDEX","CHR")]

    for (chr in unique(partition.bin.mutation.burden$CHR)) { # process 1 chr at at time to avoid multichr windows
        tmp.chr <- subset(partition.bin.mutation.burden,partition.bin.mutation.burden$CHR == chr);

        for (fun in agg.functions) {
            if (fun == "median") { gipseq.burden$MEDIAN[gipseq.burden$CHR == chr] <- ChrRollapply(tmp.chr$BIN.MUT, fun, window.size) }
            else if (fun == "mad") { gipseq.burden$MAD[gipseq.burden$CHR == chr] <- ChrRollapply(tmp.chr$BIN.MUT, fun, window.size) }
            else if (fun == "mean") { gipseq.burden$MEAN[gipseq.burden$CHR == chr] <- ChrRollapply(tmp.chr$BIN.MUT, fun, window.size) }
            else if (fun == "sd") { gipseq.burden$SD[gipseq.burden$CHR == chr] <- ChrRollapply(tmp.chr$BIN.MUT, fun, window.size) }
        }
    }

    gipseq.burden <- gipseq.burden[order(gipseq.burden$BINDEX), ]
    return(gipseq.burden)
}

example.GetClusterPartitionGipSeqCounts <- function(partition) {
    cluster.dir <- "gs://gipseq-clusters/manual/nano36.290"
    for (cluster.name in c("all","female","male")) {
        samples <- ReadClusterMembers(cluster.dir,cluster.name)

        for (GCC in c(TRUE,FALSE)) {
            for (window.size in c(100,500,1000)) {
                print(paste("Processing:",cluster.name,window.size, GCC))
    	        partition.gipseq.counts <- GetClusterPGSCIterative(samples,AddWinSize(partition,window.size),window.size,GCC)
                WritePartitionGipSeqCounts(partition.gipseq.counts,cluster.dir,cluster.name,cluster.name,window.size, GCC)
	    }
        }
    }
}

UpdateClusterPartitionGipSeqCounts <- function(cluster.dir, cluster.name,partition, window.size, GCC) {
    print(paste("Processing:",cluster.name,window.size, GCC))
    samples <- ReadClusterMembers(cluster.dir,cluster.name)
    partition.gipseq.counts <- GetClusterPGSCIterative(samples,AddWinSize(partition,window.size),window.size,GCC)
    WritePartitionGipSeqCounts(partition.gipseq.counts,cluster.dir,cluster.name,cluster.name,window.size, GCC)
}

UpdateClusterWithSamplePartitionGipSeqCounts <- function(cluster.dir, cluster.name, sample.cluster.dir, sample.cluster.name, sample.name, partition, window.size, GCC, max.in.cluster) {
    print(paste("Updating cluster",cluster.name,"with sample",sample.name))
    members <- ReadClusterMembers(cluster.dir, cluster.name)
    result <- FALSE
    old.nr.members <- nrow(members)
    is.new <- nrow(subset(members,members$CLUSTER.DIR==sample.cluster.dir & members$CLUSTER.NAME==sample.cluster.name & members$SAMPLE.NAME==sample.name)) == 0

    if (is.new) {
        if (old.nr.members < max.in.cluster) {
            members <- rbind(members,data.frame(CLUSTER.DIR=sample.cluster.dir, CLUSTER.NAME=sample.cluster.name, SAMPLE.NAME=sample.name))
            WriteClusterMembers(members,cluster.dir,cluster.name)
            print(paste("Updated nr members in",cluster.name,":", nrow(members)))
            clus.gipseq.counts <- ReadPartitionGipSeqCounts(cluster.dir, cluster.name, cluster.name, window.size,GCC)
            sam.gipseq.counts <- ReadPartitionGipSeqCounts(sample.cluster.dir, sample.cluster.name, sample.name, window.size,GCC)
            clus.gipseq.counts <- GetJointPartitionGipSeqCounts(clus.gipseq.counts,old.nr.members,sam.gipseq.counts,1, AddWinSize(partition,window.size))
            WritePartitionGipSeqCounts(clus.gipseq.counts,cluster.dir,cluster.name,cluster.name,window.size, GCC)
            result <- TRUE
        } else {
            members <- ReadClusterMembers(cluster.dir, cluster.name, extra=TRUE)
            members <- unique(rbind(members,data.frame(CLUSTER.DIR=sample.cluster.dir, CLUSTER.NAME=sample.cluster.name, SAMPLE.NAME=sample.name)))
            WriteClusterMembers(members,cluster.dir,cluster.name, extra=TRUE)
            print(paste("Members",cluster.name,"reached max of",max.in.cluster,"so no changes; updated nr extra members:", nrow(members)))
        }
    } else {
        print(paste("No changes to cluster",cluster.name,"as",old.nr.members,"already include",sample.name))
    }

    return(result)
}


GetClusterPartitionGipSeqCounts <- function(cluster.dir, cluster.name, partition, window.size, GCC) {
  #' @title Get GipSeq counts of a cluster
  #' @description Calculate aggregated partition bin counts over a fixed-sized sliding window
  #' @param cluster.dir Path to the root dir of all clusters
  #' @param cluster.name Name of the cluster
  #' @param partition Partition bin name and GC content data frame as obtained from LoadPartition
  #' @param  window.size Number of subsequent bins over which to aggregate
  #' @param GCC If true normalise the GCcorrected counts, otherwise normalise the raw counts
  #' @return  GipSeq-counts as data.frame with columns "BINDEX CHR MEAN SD",  sorted on BINDEX 
    samples <- ReadClusterMembers(cluster.dir,cluster.name)
    return(GetClusterPGSCIterative(samples,AddWinSize(partition,window.size),window.size,GCC))
}


#rclus.dir <- "gs://gipseq-clusters/manual/nano36-290"
#rclus.name <- "female"
GetClusterPGSCIterative <- function(samples,partition,window.size,GCC,start.N=0,start.gipseq.counts=NULL) {
    if (!("WIN.SIZE" %in% colnames(partition))) {partition <- AddWinSize(partition,window.size)}
    clus.gipseq.counts <- start.gipseq.counts
    first <- 1

    if (is.null(start.gipseq.counts)) {
        sam <- samples[1,]
        clus.gipseq.counts <- ReadPartitionGipSeqCounts(sam$CLUSTER.DIR, sam$CLUSTER.NAME, sam$SAMPLE.NAME,window.size,GCC)
        first <- 2
    }

    if (nrow(samples) >= first) {
        for (j in first:nrow(samples)) {
	    sam <- samples[j,]
	    print(paste("Adding sample",j,":",sam$SAMPLE.NAME))
            sam.gipseq.counts <- ReadPartitionGipSeqCounts(sam$CLUSTER.DIR, sam$CLUSTER.NAME, sam$SAMPLE.NAME,window.size,GCC)

            if (nrow(sam.gipseq.counts[!is.finite(sam.gipseq.counts$MEAN),]) > 0) {
	        print(paste(" PROBLEM WITH SAMPLE",j,":",sam$SAMPLE.NAME))
                examples.FromReadStarts(sam$CLUSTER.DIR, sam$CLUSTER.NAME, sam$SAMPLE.NAME,partition)
                sam.gipseq.counts <- ReadPartitionGipSeqCounts(sam$CLUSTER.DIR, sam$CLUSTER.NAME, sam$SAMPLE.NAME,window.size,GCC)
            } 

            if (nrow(sam.gipseq.counts[!is.finite(sam.gipseq.counts$MEAN),]) > 0) {
	        print(paste(" STILL PROBLEM WITH SAMPLE",j,":",sam$SAMPLE.NAME))
	        break
            } else {
                clus.gipseq.counts <- GetJointPartitionGipSeqCounts(clus.gipseq.counts,start.N+j-1,sam.gipseq.counts,1, partition)
            }

	    if (nrow(clus.gipseq.counts[!is.finite(clus.gipseq.counts$MEAN),]) > 0) {
	        print(" Stopping: problem with last sample")
                break
            }
        }
    }

    return(clus.gipseq.counts)
}
    

GetClusterPGSCRecursive <- function(samples,partition,window.size,GCC) {
    nr.samples <- nrow(samples)
    print(paste("Nr samples:", nr.samples))

    if (nr.samples == 1) {
        sam <- samples[1,]
        result <- ReadPartitionGipSeqCounts(sam$CLUSTER.DIR, sam$CLUSTER.NAME, sam$SAMPLE.NAME,window.size,GCC)
    } else {
        first <- nr.samples %/% 2
        result <- GetJointPartitionGipSeqCounts(GetClusGSCRecursive(samples[1:first,],partition,window.size,GCC),first,
                                                GetClusGSCRecursive(samples[(first+1):nr.samples,],partition,window.size,GCC),nr.samples-first,partition)
    }

    return(result)
}

GetJointPartitionGipSeqCounts <- function(gipseq.counts.1, nr.samples.1, gipseq.counts.2, nr.samples.2, partition) {
  #' @title Combine two partition GipSeq counts from either GipSeq clusters or GipSeq samples
  #' @description Combine the GipSeq counts of GipSeq cluster/sample 1 with the GipSeq counts of GipSeq-cluster/sample 2 and produce the GipSeq counts of GipSeq cluster 1+2; restrict stats to "MEAN" and "SD".
  #' @param gipseq.counts.1 the GipSeq-counts of the first GipSeq-cluster or GipSeq-sample; should contain at least columns "BINDEX","CHR","MEAN","SD"
  #' @param nr.samples.1 number of samples in the first GipSeq-cluster or GipSeq-sample (> 0)
  #' @param gipseq.counts.2 the GipSeq-counts of the second GipSeq-cluster or GipSeq-sample; should contain at least columns "BINDEX","CHR","MEAN","SD"
  #' @param nr.samples.2 number of samples in the second GipSeq-cluster or GipSeq-sample (> 0)
  #' @param partition Dataframe with columns "BINDEX","CHR","START","END", "WIN.SIZE" as produced by GetBindexPositions and extended by AddWinSize
  #' @return Partition GipSeq counts as data frame with columns "BINDEX","CHR","MEAN","SD"
    # calculate the nr of bins from both collections
    n1 <- nr.samples.1 * partition$WIN.SIZE
    n2 <- nr.samples.2 * partition$WIN.SIZE
    # join the 2 GipSeq-clusters, assuming both are sorted by "BINDEX"
    m <- CalculateCombinedMeanSd(n1,gipseq.counts.1$MEAN,gipseq.counts.1$SD,n2,gipseq.counts.2$MEAN,gipseq.counts.2$SD)
    # construct the new data frame
    gipseq.counts.new <- gipseq.counts.1[c("BINDEX","CHR")]
    gipseq.counts.new$MEAN <- m[,1]
    gipseq.counts.new$SD <- m[,2]
    return(gipseq.counts.new)
}

CalculateCombinedMeanSd <- function(n1,mean1,sd1,n2,mean2,sd2) {
  #' @title Calculated combined mean and standard deviation
  #' @description Combine mean and standard dev of two vectors of means and two vectors of standard deviations; the four vectors should have equal length
  #' @param n1 the number of counts aggregated to mean1 and sd1
  #' @param mean1 a first vector of means, e.g. containing one number per bin
  #' @param sd1 a first vector of standard deviations, e.g. containing one number per bin
  #' @param n2 the number of counts aggregated to mean2 and sd2
  #' @param mean2 a second vector of means, e.g. containing one number per bin
  #' @param sd2 a second vector of standard deviations, e.g. containing one number per bin
  #' @return A matrix with as many rows as the length of the input mean and sd vectors and two columns with combined stats: col1 contains the combined means, col2 the combined sds
    newMean <- (mean1*n1 + mean2*n2) / (n1+n2)
        # http://stats.stackexchange.com/questions/55999/is-it-possible-to-find-the-combined-standard-deviation
    newSD <- sqrt( ( (n1-1)*sd1**2  + (n2-1)*sd2**2 + n1*(mean1-newMean)**2 + n2*(mean2-newMean)**2 ) / (n1+n2-1) )
    return(matrix(c(newMean,newSD), ncol=2, byrow=FALSE))
}

GetSamplePartitionChrCounts <- function(partition.bin.counts, agg.functions, partition) {
  #' @title Get per chromosome counts of a sample
  #' @description Calculate aggregated partition bin counts per chromosome
  #' @param partition.bin.counts data.frame with counts per bin; columns: "BINDEX (SAMPLE)" and "COUNT"
  #' @param agg.functions The vector of aggregation functions to be applied, should be subset of: c("sum","median","mad","mean","sd")
  #' @param partition Partition bin name and GC content data frame as obtained from LoadPartition
  #' @param GCC If true normalise the GCcorrected counts, otherwise normalise the raw counts
  #' @return  Chromosome-counts as data.frame with columns "CHR" and requested columns from "MEDIAN MAD MEAN SD",  sorted on CHR 

    partition.bin.counts <- EnsureBindexHasChr(partition.bin.counts,partition)
    chromosomes <- unique(partition.bin.counts$CHR)
    chr.counts <- data.frame(CHR=chromosomes)

    for (fun in agg.functions) {
        chr.fun.counts <- aggregate(partition.bin.counts$COUNT, by=list(CHR=partition.bin.counts$CHR), FUN=fun)
	names(chr.fun.counts) <- c("CHR",fun)
	chr.counts <- merge(chr.counts, chr.fun.counts)
    }

    return(chr.counts)
}

GetSamplePartitionGipSeqCounts <- function(partition.bin.counts, agg.functions, partition, window.size) {
  #' @title Get GipSeq counts of a sample
  #' @description Calculate aggregated partition bin counts over a fixed-sized sliding window, handling chromosomes separately
  #' @param partition.bin.counts data.frame with counts per bin; columns: "BINDEX (SAMPLE)" and "COUNT"
  #' @param agg.functions The vector of aggregation functions to be applied, should be subset of: c("median","mad","mean","sd")
  #' @param partition Partition bin name and GC content data frame as obtained from LoadPartition
  #' @param  window.size Number of subsequent bins over which to aggregate
  #' @param GCC If true normalise the GCcorrected counts, otherwise normalise the raw counts
  #' @return  GipSeq-counts as data.frame with columns "BINDEX CHR" and requested columns from "MEDIAN MAD MEAN SD",  sorted on BINDEX 

    partition.bin.counts <- EnsureBindexHasChr(partition.bin.counts,partition)
    gipseq.counts <- partition[c("BINDEX","CHR")]

    for (chr in unique(partition.bin.counts$CHR)) { # process 1 chr at at time to avoid multichr windows
	#print(chr)
        tmp.chr <- subset(partition.bin.counts,partition.bin.counts$CHR == chr);

        for (fun in agg.functions) {
            if (fun == "median") { gipseq.counts$MEDIAN[gipseq.counts$CHR == chr] <- ChrRollapply(tmp.chr$COUNT, fun, window.size) }
            else if (fun == "mad") { gipseq.counts$MAD[gipseq.counts$CHR == chr] <- ChrRollapply(tmp.chr$COUNT, fun, window.size) }
            else if (fun == "mean") { gipseq.counts$MEAN[gipseq.counts$CHR == chr] <- ChrRollapply(tmp.chr$COUNT, fun, window.size) }
            else if (fun == "sd") { gipseq.counts$SD[gipseq.counts$CHR == chr] <- ChrRollapply(tmp.chr$COUNT, fun, window.size) }
        }
    }

    gipseq.counts <- gipseq.counts[order(gipseq.counts$BINDEX), ]
    return(gipseq.counts)
}

GetPartitionBinMutationBurden <- function(partition, partition.read.starts, read.len=READ.LENGTH) {
  #' @title Get partition bin mutation burden
  #' @description Calculate the average edit distance of reads per bin
  #' @param partition Partition bin name and GC content data frame obtained from LoadPartition
  #' @param partition.read.starts data.frame with 1 record per selected read and columns "CHR START NM BIN (SAMPLE)", sorted alphabetically on CHR
  #' @param read.len length of read being sequenced
  #' @return data.frame with mutation burden per bin; columns: "BINDEX BIN.MUT (SAMPLE)"  
    bin.mutation.burden <- as.data.frame(aggregate(x=partition.read.starts$NM, by=list(partition.read.starts$BIN), FUN=sum)$x*1000/(table(partition.read.starts$BIN)*read.len))
    colnames(bin.mutation.burden) <- c("BIN", "BIN.MUT")
    bin.mutation.burden <- merge(partition, bin.mutation.burden, by="BIN", all=TRUE)
    bin.mutation.burden$BIN.MUT[!is.finite(bin.mutation.burden$BIN.MUT)] <- 0
    bin.mutation.burden <- bin.mutation.burden[order(bin.mutation.burden$BINDEX), ]
    
    bin.mutation.burden <- bin.mutation.burden[c("BINDEX","BIN.MUT")]
    if ("SAMPLE" %in% colnames(partition.read.starts)) { bin.mutation.burden$SAMPLE <- head(partition.read.starts$SAMPLE,1) }

    return(bin.mutation.burden)
}

GetRawPartitionBinCounts <- function(partition,partition.read.starts) {
  #' @title Get partition bin counts
  #' @description Count the number of reads in each partition bin
  #' @param partition Partition bin name and GC content data frame as obtained from LoadPartition
  #' @param partition.read.starts data.frame with 1 record per selected read and columns "CHR  START  NM BIN (SAMPLE)", sorted alphabetically on CHR
  #' @return data.frame with counts per bin added to the partition table; columns: "BINDEX RAW.COUNT (SAMPLE)"
    partition.bin.counts <- as.data.frame(table(partition.read.starts$BIN)) #does the actual counting per bin
    names(partition.bin.counts) <- c("BIN","RAW.COUNT")
    partition.bin.counts <- merge(partition, partition.bin.counts, by="BIN",all=TRUE) #all=TRUE to make sure all bins are represented
    partition.bin.counts$RAW.COUNT[!is.finite(partition.bin.counts$RAW.COUNT)] <- 0
    partition.bin.counts <- partition.bin.counts[order(partition.bin.counts$BINDEX), ]

    partition.bin.counts <- partition.bin.counts[c("BINDEX","RAW.COUNT")]
    if ("SAMPLE" %in% colnames(partition.read.starts)) { partition.bin.counts$SAMPLE <- head(partition.read.starts$SAMPLE,1) }

    return(partition.bin.counts)
}

AddGCcorrectedPartitionBinCounts <- function(partition,partition.bin.counts) {
  #' @title Add the GC corrected counts to the partition bin counts
  #' @param partition Partition bin name and GC content data frame as obtained from LoadPartition
  #' @param partition.bin.counts data.frame with counts per bin; columns: "BINDEX RAW.COUNT (SAMPLE)"
  #' @return data.frame partition.bin.counts with the GC corrected count as an extra value
    partition.bin.counts <- merge(partition, partition.bin.counts, by="BINDEX",all=TRUE) #all=TRUE to make sure all bins are represented
    for.LOESS.calc <- subset(partition.bin.counts, partition.bin.counts$CHR != "chrX" & partition.bin.counts$CHR != "chrY" & partition.bin.counts$CHR != "chrM" & partition.bin.counts$RAW.COUNT != 0)
    median.FREQ.per.GC <- tapply(for.LOESS.calc$RAW.COUNT,for.LOESS.calc$GC.CONTENT, median)
    try(correction.factor <- predict( loess(median.FREQ.per.GC ~ as.numeric(names(median.FREQ.per.GC))), partition.bin.counts$GC.CONTENT))

    if (exists("correction.factor")) {
        partition.bin.counts$GCC.COUNT <- partition.bin.counts$RAW.COUNT * median(partition.bin.counts$RAW.COUNT)/correction.factor
        partition.bin.counts$GCC.COUNT[!is.finite(partition.bin.counts$GCC.COUNT)] <- 0
    } else {
        partition.bin.counts$GCC.COUNT <- 0
    }

    if ("SAMPLE" %in% colnames(partition.bin.counts)) {
	partition.bin.counts <- partition.bin.counts[c("BINDEX","RAW.COUNT","GCC.COUNT","SAMPLE")]
    } else {
        partition.bin.counts <- partition.bin.counts[c("BINDEX","RAW.COUNT","GCC.COUNT")]
    }

    return(partition.bin.counts)
}

AddPartitionBinCounts <- function(partition,partition.bin.counts,GCC,normalization='mean') {
  #' @title Add normalised count to partition bin counts
  #' @description The autosome-based normalised count for one of the absolute counts (raw or GC corrected) is added
  #' @param partition Partition bin name and GC content data frame as obtained from LoadPartition
  #' @param partition.bin.counts data.frame with counts per bin; columns: "BINDEX RAW.COUNT GCC.COUNT (SAMPLE)"
  #' @param GCC, if true normalise the GCcorrected counts, otherwise normalise the raw counts
  #' @param normalization The method used for normalization, can be 'mean' or 'median' of autosomal counts
  #' @return data.frame partition.bin.counts with the normalised counts as an extra column
    partition.bin.counts <- EnsureBindexHasChr(partition.bin.counts,partition)
    autosomes <- paste0("chr", seq(1,22))
    auto.partition.bin.counts <- subset(partition.bin.counts, partition.bin.counts$CHR %in% autosomes)
    
    if (GCC) {
        normalization.factor <- ifelse(normalization=='median', 
					median(auto.partition.bin.counts$GCC.COUNT, na.rm=TRUE),
					sum(auto.partition.bin.counts$GCC.COUNT)/1000000)
        partition.bin.counts$COUNT <- partition.bin.counts$GCC.COUNT / normalization.factor
    } else {
        normalization.factor <- ifelse(normalization=='median', 
					median(auto.partition.bin.counts$RAW.COUNT, na.rm=TRUE),
					sum(auto.partition.bin.counts$RAW.COUNT)/1000000)
        partition.bin.counts$COUNT <- partition.bin.counts$RAW.COUNT / normalization.factor
    } 

    if ("SAMPLE" %in% colnames(partition.bin.counts)) {
	partition.bin.counts <- partition.bin.counts[c("BINDEX","COUNT","SAMPLE")]
    } else {
	partition.bin.counts <- partition.bin.counts[c("BINDEX","COUNT")]
    }

    partition.bin.counts$COUNT[!is.finite(partition.bin.counts$COUNT)] <- 0
    return(partition.bin.counts)
}


GetChrDensity <- function(partition.read.starts,bandwidth,regionsinfo.file,itp.file=NA,autosomeonly=FALSE,exclude.centromere=TRUE) {
  #' @title Get density from read starts distribution
  #' @description Generate density of a sample from read starts using Gaussian kernel
  #' @param partition.read.starts data.frame with 1 record per selected read and columns "CHR START NM BIN (SAMPLE)", sorted alphabetically on CHR
  #' @param bandwidth Bandwidth of smoothing kernel
  #' @param regionsinfo.file The file containing chromosome(region) information
  #' @param itp.file Bed file that contains a list of genomic positions to be interpolated from the density
  #' @param autosomeonly Whether generate density only for autosomes
  #' @param exclude.centromere Exclude sampling positions that fall in centromere regions

  #set bandwidth for density calculation; bw: smoothing kernel bandwidth; final bandwidth to be used in density calculation: adjust*bandwidth
    bandwidth <- as.numeric(bandwidth)
    #print(paste("bandwidth set to: ", bandwidth))

    regionsinfo <- LoadRegionInfo(regionsinfo.file)

    if (autosomeonly) {
        chroms <- paste0("chr", seq(1,22))
        nn <- nrow(regionsinfo)
    } else {
        chroms <- paste0("chr", seq(1,22))
        chroms <- c(chroms, "chrX", "chrY")
        nn <- nrow(regionsinfo)
    }

    if (is.na(itp.file)) {
        denlist <- lapply(1:nn, function(x) {
            denchr <- regionsinfo$CHR[x]
            denstart <- regionsinfo$START[x]
            denend <- regionsinfo$END[x]
            centrostart <- regionsinfo$CENSTART[x]
            centroend <- regionsinfo$CENEND[x]
            readstartschr <- subset(partition.read.starts, CHR==denchr)

            #testwidth can be changed; testwidth infers the roughly distance between two sampling points
            #calculate grid size based on size of chromosome and sampling points
            testwidth <- 50000
            npower <- round(log((denend-denstart)/testwidth)/log(2))
            npower <- max(npower, 1)
            ninterpolation <- 2^npower
            regionden <- density(readstartschr$START, adjust=1, bw=bandwidth, from=denstart, to=denend, n=ninterpolation)
            denest <- as.data.frame(regionden$y)

            #nearest integer position; might not be integer position with fixed width interpolation
            pos <- as.data.frame(round(regionden$x, digits=0))
            estdf <- cbind(pos, denest)
            estdf$CHR <- denchr
            colnames(estdf) <- c("POS","DENSITY","CHR")
            if (exclude.centromere) {
                estdf <- estdf[estdf$POS < centrostart | estdf$POS > centroend,]
            }
            return(estdf)
        })
    } else {
        #specify interpolation position file
        itpos <- LoadInterpolationPos(itp.file)

        denlist <- lapply(1:nn, function(x) {
            denchr <- regionsinfo$CHR[x]
            denstart <- regionsinfo$START[x]
            denend <- regionsinfo$END[x]
            readstartschr <- subset(partition.read.starts, CHR==denchr)
            itposchr <- subset(itpos, CHR==denchr & START>=denstart & END<=denend)
            testwidth <- 50000
            npower <- round(log((denend-denstart)/testwidth)/log(2))
            npower <- max(npower, 1)
            ninterpolation <- 2^npower
            regionden <- approxfun(density(readstartschr$START, adjust=1, bw=bandwidth, from=denstart, to=denend, n=ninterpolation))
            denest <- as.data.frame(regionden(itposchr$END))
            pos <- as.data.frame(itposchr$END)
            estdf <- cbind(pos, denest)
            estdf$CHR <- denchr
            colnames(estdf) <- c("POS","DENSITY","CHR")
            return(estdf)
        })
    }

    denest <- Reduce(function(x, y) {rbind(x, y)}, denlist)
    denest <- denest[,c("CHR", "POS", "DENSITY")]
    denest$DENSITY[is.na(denest$DENSITY)] <- 0
    return(denest)
}

GetConcateDensity <- function(partition.read.starts,bandwidth,concateinfo,itp.file=NA,autosomeonly=TRUE,exclude.centromere=TRUE) {
  #' @title Get density from concatenated chromosome read starts distribution
  #' @description Generate density of a sample from read starts using Gaussian kernel
  #' @param partition.read.starts data.frame with 1 record per selected read and columns "CHR START NM BIN (SAMPLE)", sorted alphabetically on CHR
  #' @param cluster.dir Path to the root dir of all clusters
  #' @param cluster.name Name of the cluster
  #' @param sample.name Name of the sample
  #' @param bandwidth Bandwidth of smoothing kernel
  #' @param concateinfo Data frame containing concatenated chromosome(region) information
  #' @param itp.file Bed file that contains a list of genomic positions to be interpolated from the density
  #' @param autosomeonly Whether generate density only for autosomes
  #set bandwidth for density calculation; bw: smoothing kernel bandwidth; final bandwidth to be used in density calculation: adjust*bandwidth
    bandwidth <- as.numeric(bandwidth)
    partition.read.starts <- partition.read.starts[,c("CHR","START")]
    ##concate only autosomes
    if (autosomeonly) {
        chroms <- paste0("chr", seq(1,22))
        nn <- length(chroms)
    } else {
        chroms <- paste0("chr", seq(1,22))
        chroms <- c(chroms, "chrX", "chrY")
        nn <- length(chroms)
    }

    concateinfo <- concateinfo[concateinfo$CHR %in% chroms,]
    concateinfo$CONSTART <- as.numeric(concateinfo$CONSTART)
    concateinfo$CONEND <- as.numeric(concateinfo$CONEND)
    concateinfo$START <- as.numeric(concateinfo$START)
    denstart <- min(concateinfo$CONSTART)
    denend <- max(concateinfo$CONEND)

    concaterslist <- lapply(1:nn, function(x) {
        denchr <- concateinfo$CHR[x]
        denconstart <- concateinfo$CONSTART[x]
        chrstart <- concateinfo$START[x]
        readstartschr <- subset(partition.read.starts, CHR==denchr)

        #minus chrstart in case when the chromosome is acrocentric (13,14,15,21,22)
        rsconchr <- readstartschr$START + denconstart - chrstart
        rsconchr <- as.data.frame(rsconchr)
        conchrdf <- cbind(readstartschr, rsconchr)
        colnames(conchrdf) <- c("CHR", "START", "CONSTART")
        return(conchrdf)
    })
    conReadStarts <- Reduce(function(x, y) {rbind(x, y)}, concaterslist)
    conReadStarts$CONSTART <- as.numeric(conReadStarts$CONSTART)

    testwidth <- 50000
    npower <- round(log((denend-denstart)/testwidth)/log(2))
    npower <- max(npower, 1)
    ninterpolation <- 2^npower

    if (is.na(itp.file)) {
        catden <- density(conReadStarts$CONSTART, adjust=1, bw=bandwidth, from=denstart, to=denend,n=ninterpolation)
        pos1 <- round(catden$x, digits=0)
        pos2 <- pos1 + 1
        estdf <- data.frame(pos1, pos2, catden$y)
        colnames(estdf) <- c("CONSTART", "CONEND", "DENSITY")
        estdf <- data.table(estdf)

        #make sure start and end in 64 bit
        estdf$CONSTART <- as.integer64(estdf$CONSTART)
        estdf$CONEND <- as.integer64(estdf$CONEND)

        #find overlaps between concatenated density and concatenated chromosomes to translate positions back to chromosomal genomic positions

        concateinfo <- data.table(concateinfo)
        concateinfo$CONSTART <- as.integer64(concateinfo$CONSTART)
        concateinfo$CONEND <- as.integer64(concateinfo$CONEND)
        setkey(concateinfo, CONSTART, CONEND)
        dfop <- foverlaps(estdf, concateinfo, type="within", nomatch=NULL)
        dfop$POS <- dfop$i.CONSTART-dfop$CONSTART+dfop$START
        fdf <- dfop[,c("CHR","POS","DENSITY")]
        if (exclude.centromere) {
            excenlist <- lapply(1:nrow(concateinfo), function(x) {
                curchr <- concateinfo$CHR[x]
                centrostart <- concateinfo$CENSTART[x]
                centroend <- concateinfo$CENEND[x]
                curdf <- subset(fdf, CHR==curchr & (POS < centrostart | POS > centroend))
                return(curdf)
            })
            fdf <- Reduce(function(x, y) {rbind(x, y)}, excenlist)
        }
        fdf$DENSITY[is.na(fdf$DENSITY)] <- 0
    } else {
        itpos <- LoadInterpolationPos(itp.file)
        itpos <- itpos[itpos$CHR %in% chroms,]
        itpos$CONSTART <- as.numeric(itpos$CONSTART)

        catden <- approxfun(density(conReadStarts$CONSTART, adjust=1, bw=bandwidth, from=denstart, to=denend,n=ninterpolation))
        denest <- as.data.frame(catden(itpos$CONSTART))

        estdf <- cbind(itpos, denest)
        colnames(estdf) <- c("CHR","POS","END","CONSTART","CONEND","DENSITY")
        fdf <- estdf[,c("CHR", "POS", "DENSITY")]
        fdf$DENSITY[is.na(fdf$DENSITY)] <- 0
    }
    return(fdf)
}


GetPartitionReadStarts <- function(samtools,crambam.file,partition,sample.name="NA",skip.secondary.reads=TRUE, skip.duplicates=TRUE, min.mapping.quality=20, use.fragment.starts=TRUE, whitelist.bed=NA, verbose=FALSE) {
  #' @title Get start positions and a partition bin id of selected reads from bam (or cram) file
  #' @description Iterate over all bam file entries that meet the filter criteria and record the start position (and bin name) of each read
  #' @param samtools Full path to the samtools command
  #' @param crambam.file Full path to the cram or bam file of the sequenced sample
  #' @partition Partition bin name and GC content data frame as obtained from LoadPartition
  #' @param sample.name Sample name to be added to the result table
  #' @param skip.secondary.reads Flag -F256; those reads are mapped twice, usually a very small nr; they are skipped to make sure the totals match 
  #' @param skip.duplicates Flag -F1024; skip marked duplicates; assumes duplicates have been marked in the bam file
  #' @param min.mapping.quality Option -q used to exclude reads with bad mapping quality
  #' @param use.fragment.starts Correct read starting position of reverse reads by adding length of read, so that all starting positions indicate fragment starts (=cutting point cfDNA)
  #' @param whitelist.bed Bed file with regions to include
  #' @param verbose If TRUE print intermediate statements on progress
  #' @return Dataframe with 1 line per selected read and columns "CHR  START  NM BIN SAMPLE", sorted alphabetically on CHR
  #     CHR: chromosome name
  #     START: 1-based coordinate of read start
  #     NM: edit distance from NM tag in bam file
  #     BIN: id of the bin in format 'CHR.ID', with ID being the zero-based number of the bin within CHR
    view.command <- GetSamtoolsViewCommand(samtools,crambam.file, skip.secondary.reads, skip.duplicates, min.mapping.quality, use.fragment.starts,whitelist.bed,whitelist.bed=whitelist.bed)

    #NM tag column number might change if alignment changes; find NM column number
    NMcol <- as.numeric(fread(cmd=paste(view.command, "| head -n 1 | awk -v RS=\"\\t\" '/^NM:/{print NR}'"))[1,1])
    colstokeep <- paste(3, 4, NMcol, sep = ",")

    if(verbose) {print(paste("Extract all CHR,START,NM records from the bamfile, using command:", view.command))}
    bin.reads <- try(fread(cmd=paste(view.command, "| cut -f", colstokeep), header=FALSE, colClasses= c("factor", "integer", "character"),data.table=FALSE), silent=TRUE)

    if (!is.data.frame(bin.reads)) { #in case no reads, create empty dataframe
        bin.reads <- data.frame(CHR=character(), START=integer(), NM=character(), stringsAsFactors=FALSE)
    }

    names(bin.reads)=c("CHR","START","NM")
    bin.reads <- subset(bin.reads,bin.reads$CHR %in% unique(partition$CHR)) #typically: select chr 1-22, X,Y   and remove MT,contigs, etc...
    bin.reads <- bin.reads[order(as.character(bin.reads$CHR)), ]
    
    bin.reads$NM <- gsub("NM:i:", "", bin.reads$NM)
    bin.reads$NM <- as.numeric(as.character(bin.reads$NM))

    # here the actual binning happens
    bin.reads <- AddBinToPartitionReadStarts(bin.reads)
    #bin.reads$BIN <- paste(bin.reads$CHR, ifelse((bin.reads$START %% BIN.SIZE)==0, bin.reads$START %/% BIN.SIZE - 1, bin.reads$START %/% BIN.SIZE), sep = ".")

    if (!is.na(sample.name)) { bin.reads$SAMPLE <- sample.name }

    return(bin.reads)
}

AddBinToPartitionReadStarts <- function(partition.read.starts) {
    partition.read.starts$BIN <- paste(partition.read.starts$CHR, ifelse((partition.read.starts$START %% BIN.SIZE)==0, partition.read.starts$START %/% BIN.SIZE - 1, partition.read.starts$START %/% BIN.SIZE), sep = ".")
    return(partition.read.starts)
}

############
# ROI BINS #
############

GetRoiBinMutationBurden <- function(partition.read.starts, roi, read.len=36) {
  #' @title Get ROI bin mutation burden
  #' @description Calculate the average edit distance of reads per ROI
  #' @param partition.read.starts data.frame with 1 record per selected read and columns "CHR  START  NM BIN (SAMPLE)", sorted alphabetically on CHR
  #' @param roi data.frame with regions of interest; columns: CHR, START, END, DESCRIPTION
  #' @return data.frame with columns: CHR, START, END, DESCRIPTION, BIN.MUT
    BIN.MUT <- apply(roi, 1, 
		function(x){
		  prs <- subset(partition.read.starts,partition.read.starts$CHR==x[1] &
					   partition.read.starts$START>=as.numeric(x[2]) &
					   partition.read.starts$START<as.numeric(x[3]))
		  result <- 0 
		  if (nrow(prs) > 0) {
		    result <- sum(prs$NM) * 1000/(nrow(prs)*read.len)
		  }

		  return(result)
		}
	       )
    
    roi.mut.df <- cbind(roi, BIN.MUT)
    roi.mut.df <- roi.mut.df[order(roi.mut.df$CHR, roi.mut.df$START, roi.mut.df$END), ]
    
    return(roi.mut.df)
}

GetRoiBinCounts <- function(partition.read.starts, roi) {
  #' @title Get the normalised ROI bin counts for a sample
  #' @description Calculate the number of reads and RPKM for a list of regions of interest
  #' @param samtools Full path to the samtools command
  #' @param crambam.file Full path to the cram or bam file of the sequenced sample
  #' @param roi data.frame with regions of interest; columns: CHR, START, END, DESCRIPTION
  #' @param chr.size data.frame with chromosome sizes; columns: CHR, SIZE
  #' @return data.frame with columns: CHR,START,END,DESCRIPTION,COUNT,RPKM,
    totmapped <- 0

    totmapped <- nrow(partition.read.starts)
    #print(paste("Found", totmapped, "mapped reads"))
    permillion <- totmapped/1000000
    #roi.stats.df <- data.frame(CHR=character(),START=integer(),END=integer(),DESCRIPTION=character(),COUNT=integer(),RPKM=numeric())

    COUNT <- apply(roi, 1, function(x){nrow(partition.read.starts[(partition.read.starts$CHR==x[1] & partition.read.starts$START>=as.numeric(x[2]) & partition.read.starts$START<as.numeric(x[3])),])})
    RPKM <- COUNT/(permillion*((roi$END-roi$START)/1000))
    roi.stats.df <- cbind(roi, COUNT, RPKM)

    roi.stats.df <- roi.stats.df[order(roi.stats.df$CHR,roi.stats.df$START,roi.stats.df$END), ]

    return(roi.stats.df)
}


#############
# UTILITIES #
#############

GetSamtoolsViewCommand <- function(samtools,crambam.file, skip.secondary.reads, skip.duplicates, min.mapping.quality, use.fragment.starts,whitelist.bed=NA,region.of.interest=NA) {
  #' @title Construct samtools view command
  #' @description Use parameter setting to construct samtools command for selection/counting of reads
  #' @param samtools Full path to the samtools command
  #' @param crambam.file Full path to the cram or bam file of the sequenced sample
  #' @param skip.secondary.reads Flag -F256; those reads are mapped twice, usually a very small nr; they are skipped to make sure the totals match 
  #' @param skip.duplicates Flag -F1024; skip marked duplicates; assumes duplicates have been marked in the bam file
  #' @param min.mapping.quality Option -q used to exclude reads with bad mapping quality
  #' @param use.fragment.starts Correct read starting position of reverse reads by adding length of read, so that all starting positions indicate fragment starts (=cutting point cfDNA)
  #' @param whitelist.bed Bed file with regions to include
  #' @region.of.interest Region of interest in format CHR:START-END
  #' @return Samtools command

    view.command <- paste(samtools,"view")
    if (!is.na(region.of.interest)) {view.command <- paste(view.command,"-c")}
    if (skip.secondary.reads) {view.command <-  paste(view.command,"-F 256")}
    if (skip.duplicates) {view.command <-  paste(view.command,"-F 1024")}
    if (min.mapping.quality > 0) {view.command <- paste(view.command,"-q",min.mapping.quality)}
    if (!is.na(whitelist.bed)) {view.command <- paste(view.command,"-L",whitelist.bed)}

    if (startsWith(crambam.file,"gs://")) {
        view.command <- paste("gsutil cat", crambam.file, "|", view.command)
    } else {
        view.command <- paste(view.command, crambam.file)
    }

    if (is.na(region.of.interest)) {
        if (use.fragment.starts) { #add read length to reverse reads' start position, to get to the start of the fragment
	    view.command <- paste(view.command, "| awk -F'\t' 'BEGIN{OFS=FS}{if(and($2,0x10)==16){$4=$4+length($10)-1};print($0)}'")
        }
    } else {
	view.command <- paste(view.command, region.of.interest) # format 1:1-10 to get first 11 positions from chr1
    }

    return(view.command)
}

GetAllROI <- function() {
  #' @title Get a lost of all regions of interest
  #' @description Combine listed ROI with chromosome arms and full chromosomes
  #' @return data.frame with info on all ROIs in columns: "CHR START END DESCRIPTION"
    roi <- LoadRegionsOfInterest(roi.file)
    chr.arm <- LoadRegionsOfInterest(chr.arm.file)
    chr.size <- LoadChrSize(chr.size.file)
    result.df <- rbind(roi,chr.arm)

    for (i in 1:nrow(chr.size)) {
	chr.info <- chr.size[i,]
	result.df <- rbind(result.df, data.frame(CHR=chr.info$CHR, START=1, END=chr.info$SIZE, DESCRIPTION=chr.info$CHR, stringsAsFactors=FALSE))
    }

    result.df <- result.df[order(result.df$CHR, result.df$START, result.df$END), ]
    return(result.df)
}
