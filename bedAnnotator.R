#!/usr/local/bin/Rscript
####################################################################################################################
#
# This script annotates regions in bed files (e.g. from ChiP-Seq) with supplied annotation files.
# @author Jens Preußner<jens.preussner@mpi-bn.mpg.de>
# @version 1.0
# @license MIT
# @requirements getopt, ChIPpeakAnno
#
# Copyright (c) 2014, Jens Preussner
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
#####################################################################################################################

require('getopt');
require('ChIPpeakAnno');

options = matrix(c(
  'peaks','p',1,'character','A bed file with peaks to be annotated (required).',
  'help','h',0,'logical','Provides command line help.',
  'gff','g',2,'character','A gff file with the annotation data (null).',
  'bed','b',2,'character','A bed file with the annotation data (null).',
  'gtf','f',2,'character','A gtf file with the annotation data (null).',
  'annotation','a',2,'character','A short description for the annotation (annotation).',
  'strategy','s',2,'character','Annotation strategy. Can be overlap or nearest (overlap).',
  'maxgap','m',2,'integer','Intervals are allowed to be separated by at most <maxgap> bases with strategy overlap (5000).',
  'id_attribute','i',2,'character','Which gff3/gtf attribute (column 9) to use as feature ID (ID).',
  'select_strategy','t',2,'character','Define how to deal with multiple overlaps [all,first,last,arbitrary] (all).',
  'peaklocation','u',2,'character','Define location in peak for distance calculation in nearest strategy [middle,start] (middle).',
  'featurelocation','v',2,'character','Define location in feature for distance calculation in nearest strategy [middle, TSS, geneEnd] (TSS).',
  'bidirectional','d',0,'logical','Indicates wether the annotation can be treated as promoters with strategy overlap (false).',
  'peakheader','k',0,'logical','Set if peak file has a header line (false).',
  'annotationheader','l',0,'logical','Set if annotation file (gff or bed) has a header line (false).',
  'omit_full_annotation','n',0,'logical','Omit the full annotation as extra column from GFF/GTF files (false).',
  'writethrough','w',0,'logical','Set if write-through from stdin (false).'),ncol=5,byrow=T)

opt = getopt(options);

####################################################################################################################
#
# HELPER FUNCTIONS
#
####################################################################################################################
parseIDFromGffAttributes = function(x,i) {
  unlist(lapply(strsplit(x,";",fixed=T),function(y) {
    d = unlist(strsplit(y,"=",fixed=T))
    t = which(d==i)
    if(length(t) == 1) { return(d[t+1]); }
    else { return(NA); }
  }))
}

parseIDFromGftAttributes = function(x,i) {
  unlist(lapply(strsplit(x,";",fixed=T),function(y) {
    d = unlist(strsplit(y," ",fixed=T))
    t = which(d==i)
    if(length(t) == 1) { return(d[t+1]); }
    else { return(NA); }
  }))
}

printOverlap = function(overlap, format, strategy, omit_additional_infos=FALSE) {
  strand = overlap$strand;
  if(strategy == "overlap-nonbi") {
    chr = overlap$chr;
    start = overlap$peaks_start;
    end = overlap$peaks_end;
    peaks_name = overlap$peaks;
    feature_name = as.character(overlap$annotation);
    feature_strand = overlap$strand1;
    feature_position = overlap$overlapFeature;
    feature_start = overlap$annotation_start;
    feature_end = overlap$annotation_end;
  }
  else {
    chr = as.character(overlap$space);
    start = start(overlap);
    end = end(overlap);
    peaks_name = overlap$peak;
    feature_name = overlap$feature;
    feature_strand = overlap$strand;
    feature_position = overlap$insideFeature;
    feature_start = overlap$start_position;
    feature_end = overlap$end_position;
  }

  if(format == "gff") {
    AttributesData = read.table(file.path(opt$gff),header=opt$annotationheader, sep="\t", colClasses = c(rep("NULL",8),"character"))
    attributes = parseIDFromGffAttributes(AttributesData$V9,opt$id_attribute)
  }
  if(format == "gtf") {
    AttributesData = read.table(file.path(opt$gtf),header=opt$annotationheader, sep="\t", colClasses = c(rep("NULL",8),"character"))
    attributes = parseIDFromGftAttributes(AttributesData$V9,opt$id_attribute)
  }
  if(format %in% c("gff","gtf")) {
    fullAnnotation = AttributesData[as.numeric(feature_name),"V9"];
    feature_name = attributes[as.numeric(feature_name)];
  }

  if(format %in% c("gff","gtf") && !omit_additional_infos) {
    cat(paste(chr,start,end,peaks_name,strand,opt$annotation,feature_start,feature_end,feature_name,feature_strand,feature_position,fullAnnotation,sep="\t"),sep="\n")
  }
  else {
    cat(paste(chr,start,end,peaks_name,strand,opt$annotation,feature_start,feature_end,feature_name,feature_strand,feature_position,sep="\t"),sep="\n")
  }
}

####################################################################################################################
#
# SCRIPT STARTS HERE
#
####################################################################################################################
# help was asked for.
if ( !is.null(opt$help) ) {
  cat(getopt(options, usage=TRUE),file = stderr());
  q(status=1);
}

# set some reasonable defaults for the options that are needed,
# but were not specified.
if ( is.null(opt$omit_full_annotation) ) { opt$omit_full_annotation = FALSE; }
if ( is.null(opt$strategy) ) { opt$strategy = 'overlap'; }
if ( is.null(opt$select_strategy) ) { opt$select_strategy = 'all'; }
if ( is.null(opt$peaklocation) ) { opt$peaklocation = 'middle'; }
if ( is.null(opt$featurelocation) ) { opt$featurelocation = 'TSS'; }
if ( is.null(opt$id_attribute) ) { opt$id_attribute = 'ID'; }
if ( is.null(opt$annotation) ) { opt$annotation = 'annotation'; }
if ( is.null(opt$bidirectional) ) { opt$bidirectional = FALSE; }
if ( is.null(opt$peakheader) ) { opt$peakheader = FALSE; }
if ( is.null(opt$annotationheader) ) { opt$annotationheader = FALSE; }
if ( is.null(opt$writethrough) ) { opt$writethrough = FALSE; }
if ( is.null(opt$maxgap) ) { opt$maxgap = 5000; }
if ( is.null(opt$peaks) ) {
  cat(getopt(options, usage=TRUE),file = stderr());
  q(status=2);
}

#Sanity check
if( !opt$select_strategy %in% c("all","first","last","arbitrary") ) {opt$select_strategy = 'first'; }
if( !opt$peaklocation %in% c("middle","start") ) { opt$peaklocation = 'first'; }
if( !opt$featurelocation %in% c("middle","TSS","geneEnd")) { opt$featurelocation = 'TSS'; }
if( !opt$strategy %in% c("overlap","nearest")) { opt$strategy = 'overlap'; }

# Load peaks in BED format
if (!is.null(opt$peaks)) {
  peaks = BED2RangedData(file.path(opt$peaks),header=opt$peakheader)
}

annotationFormat = ifelse((!is.null(opt$gff) || !is.null(opt$gtf) && is.null(opt$bed)),ifelse(!is.null(opt$gff),"gff","gtf"),"bed")

# Load annotation data
if (annotationFormat %in% c("gff","gtf")) {
  if (annotationFormat == "gff") { gffgtf = opt$gff; }
  if (annotationFormat == "gtf") { gffgtf = opt$gtf; }
  annotation = GFF2RangedData(file.path(gffgtf),header=opt$annotationheader,comment.char="#",sep="\t")
} else {
  if (annotationFormat == "bed") {
    annotation = BED2RangedData(file.path(opt$bed),header=opt$annotationheader,sep="\t")
  } else {
    cat("None or more than one annotation files supplied, cannot decide which to take...\n",file = stderr())
    if (opt$writethrough) { q(status=2); }
  }
}

# Perform overlap or nearest-search depending on strategy
if( opt$strategy == 'overlap' ) {
  if ( !opt$bidirectional ) {
    ovl_data = findOverlappingPeaks(peaks, annotation, maxgap=opt$maxgap, select=opt$select_strategy, annotate=1, NameOfPeaks1="peaks", NameOfPeaks2="annotation")
    printOverlap(ovl_data$OverlappingPeaks, annotationFormat, 'overlap-nonbi', opt$omit_full_annotation)
  } else {
    ovl_data = peaksNearBDP(peaks, AnnotationData=annotation, MaxDistance=opt$maxgap, PeakLocForDistance = "middle", FeatureLocForDistance = "TSS")
    printOverlap(ovl_data$peaksWithBDP, annotationFormat, 'overlap-bi', opt$omit_full_annotation)
  }
} else {
  if ( opt$strategy == 'nearest' ) {
    ovl = annotatePeakInBatch(peaks, AnnotationData = annotation, maxgap = opt$maxgap, select=opt$select_strategy, PeakLocForDistance = opt$peaklocation, FeatureLocForDistance = opt$featurelocation)
    printOverlap(ovl, annotationFormat, 'nearest', opt$omit_full_annotation)
  } else {
    cat(paste("Unrecognized strategy",opt$strategy,"..."),file = stderr())
    if (opt$writethrough) { q(status=2); }
  }
}

####################################################################################################################
#
# std-in write through
#
####################################################################################################################
# Enable write-through. All input that was piped to the tool will be forwarded.
if (opt$writethrough) {
  f <- file("stdin")
  open(f)
  while(length(line <- readLines(f,n=1)) > 0) {
    write(line, stderr())
  }
}

