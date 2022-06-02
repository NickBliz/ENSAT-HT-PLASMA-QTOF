##################################################
## Title: xcmsAlignment.R
## Description: The script performs alignment using xcms and feature grouping and annotation using CAMERA. 
## Input variables:
## inputDirPath: Directory path containing the data files
## outputDirPath: Directory path where the feature table file will be exported
## fileType: A string containing the data file type 
## acquisitionMode: "positive" or "negative"
## Example run:
## To process mzData file:
## xcmsAlignment("D:/TML/Metabolomics/NGMS_v_0_5/XCMS_using_R/170421_mzData/", "D:/TML/Metabolomics/NGMS_v_0_5/XCMS_using_R/", "*.xml", "positive")
## To process mzML files
## xcmsAlignment("D:/TML/Metabolomics/NGMS_v_0_5/XCMS_using_R/170421_mzML/", "D:/TML/Metabolomics/NGMS_v_0_5/XCMS_using_R/", "*.mzML", "positive")
## Date: Tue Jan 15 13:41:06 2019
## Author: Purva Kulkarni
##################################################

library(xcms)
library(RColorBrewer)
library(multtest)
library(CAMERA)



xcmsAlignment <- function(inputDirPath, outputDirPath, fileType, acquisitionMode)
{
  # List absoulte paths of raw data files to be processed
  dataFiles <-
    list.files(
      inputDirPath,
      pattern = fileType,
      recursive = FALSE,
      full.names = TRUE
    )
  
  # Set outputDirPath as the current directory
  setwd(outputDirPath)
  
  # For EIC construction and peak picking in batch mode
  # xcmsSet: This function handles the construction of xcmsSet objects. It finds peaks in batch mode and pre-sorts files from subdirectories into different classes suitable for grouping
  # Feature detection method: centwave - The centWave algorithm perform peak density and wavelet based chromatographic peak detection for high resolution LC/MS data in centroid mode
  # ppm: maximal tolerated m/z deviation in consecutive scans, in ppm (parts per million) for initial ROI definition
  # peakwidth: range with minimum and maximum chromatographic peak width in seconds. c(min,max)
  # snthresh: Signal/Noise ratio threshold
  # mzdiff: minimum difference in m/z for peaks with overlapping retention times, can be negative to allow overlap. During peak post-processing, peaks defined to be overlapping are reduced to the one peak with the largest signal.
  # integrate: Integration method. If =1 peak limits are found through descent on the mexican hat filtered data, if =2 the descent is done on the real data. Method 2 is very accurate but prone to noise, while method 1 is more robust to noise but less exact.
  # prefilter: Prefilter step for the first analysis step. Mass features are only retained if they contain at least k peaks with intensity >= I. c(k, I)
  # noise: optional argument which is useful for data that was centroided without any intensity threshold. Allows to set a minimum intensity required for centroids to be considered in the first analysis step. Centroids with intensity < noise are omitted from ROI detection
  
  xset <- xcmsSet(
    dataFiles,
    method = "centWave",
    ppm = 15,
    peakwidth = c(5, 10),
    snthresh = 10,
    mzdiff = 0.01,
    integrate = 2,
    prefilter = c(3, 500),
    noise = 0
  )
  
  # Retention time correction
  # retcor: Used to correct differences between retention times between different samples. The obiwrap retention time correction method calculates retention time deviations for each sample.
  # It is based on the code at http://obi-warp.sourceforge.net/. However, this function is able to align multiple samples, by a center-star strategy.
  # profStep: step size (in m/z) to use for profile generation from the raw data files
  
  xset_retcor <- retcor(xset,
                        method = "obiwarp",
                        profStep = 1)
  
  # Alignment
  # Peaks representing the same analyte are grouped together across samples using overlapping m/z bins and calculation of smoothed peak distributions in chromatographic time
  # bw: Allowable retention time deviations, in seconds. In more detail: bandwidth (standard deviation or half width at half maximum) of gaussian smoothing kernel to apply to the peak density chromatogram
  # minfrac: minimum fraction of samples necessary in at least one of the sample groups for it to be a valid group
  # mzwid: width of overlapping m/z slices to use for creating peak density chromatograms and grouping peaks across samples
  # minsamp: minimum number of samples necessary in at least one of the sample groups for it to be a valid group
  # max: maximum number of groups to identify in a single m/z slice
  
  xset_retcor_group <- group(
    xset_retcor,
    bw = 3,
    minfrac = 0.1,
    mzwid = 0.1,
    minsamp = 1,
    max = 1000
  )
  
  # Fill peaks
  # For each sample, identify peak groups where that sample is not represented. For each of those peak groups, integrate the signal in the region of that peak group and create a new peak.
  # After peak grouping, there will always be peak groups that do not include peaks from every sample. This method produces intensity values for those missing samples by integrating raw data in peak group region.
  # According to the type of raw-data there are 2 different methods available. For filling gcms/lcms data the method "chrom" integrates raw-data in the chromatographic domain.
  # Whereas "MSW" is used for peaklists without retention-time information like those from direct-infusion spectra.
  
  xset_retcor_group_fillPeaks <- fillPeaks(xset_retcor_group)
  
  # Construct xsAnnotate object
  # This function transforms a xcmsSet object with peaks from multiple LC/MS samples into a set of annotation results.
  
  xset_retcor_group_fillPeaks_annotate <-
    xsAnnotate(xset_retcor_group_fillPeaks)
  
  # Group peaks
  # This function groups peaks of an xsAnnotate object according to their retention time into pseudospectra-groups. It returns xsAnnotate object with pseudospectra information.
  
  xset_retcor_group_fillPeaks_annotate_group <-
    groupFWHM(xset_retcor_group_fillPeaks_annotate)
  
  # Find isotopes
  # Annotates isotopes peaks for xsAnnotate object according to a set of rules. This generates a list of all possible isotopes, which is stored in object@isotopes. This isotope information will be used in the groupCorr funtion in the next step.
  
  xset_retcor_group_fillPeaks_annotate_group_isotopes <-
    findIsotopes(xset_retcor_group_fillPeaks_annotate_group,
                 ppm = 5,
                 mzabs = 0.015)
  
  # High correlation peak grouping
  # This function groups peaks that have high correlation amoungst each other based on correlation of intensities across samples (needs more than 3 samples), EIC correlation between peaks inside a sample and isotope cluster information.
  
  xset_retcor_group_fillPeaks_annotate_group_isotopes_corr <-
    groupCorr(xset_retcor_group_fillPeaks_annotate_group_isotopes)
  
  # Find adducts
  # This function annotates adducts peaks and fragments for a xsAnnotate object.
  
  xset_retcor_group_fillPeaks_annotate_group_isotopes_corr_adducts <-
    findAdducts(
      xset_retcor_group_fillPeaks_annotate_group_isotopes_corr,
      ppm = 5,
      mzabs = 0.015,
      polarity = acquisitionMode
    )
  
  # Extract annotated feature table
  
  annotatedPeakList <- getPeaklist(
    xset_retcor_group_fillPeaks_annotate_group_isotopes_corr_adducts
  )
  
  # Change required column names
  colnames(annotatedPeakList)[1] <- "mzmed"
  colnames(annotatedPeakList)[4] <- "rtmed"
  
  # Convert retention time in minute format
  
  annotatedPeakList[4] <- annotatedPeakList[4]/60
  annotatedPeakList[5] <- annotatedPeakList[5]/60
  annotatedPeakList[6] <- annotatedPeakList[6]/60
  
  # Export annotated feature table to a .tsv file
   
  write.table(
    annotatedPeakList,
    file = "XCMSR.annotated.Report.tsv", col.names = NA, sep = "\t", quote = FALSE
  )
}