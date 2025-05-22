########################################################################################
### required installations, libraries and workplaces
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("MSnbase")
library("MSnbase")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("RforProteomics")
library("RforProteomics")
########################################################################################
### Values to define
### MSMS tolerance in Da - might be adapted based on the instrument used
### Minimum mass shift denotes filter for mass shifts (in Da) to be eliminated if they are too low or negative
### Minimum fragement number indicate filter for numbe of reporter ions, can be adapted accordingly
FragTol<-0.01
FragTol<-FragTol^2
Min_mass_shift<-40
Min_frag_number<-8
########################################################################################
########################################################################################
### reading in mgf files - replace filename with the 
MR1_data<-readMgfData(filename = "HFX221_IPP345_TTS_Enas_Rep0_Msmeg_Infection_negative_control.mgf", 
                      pdata = NULL, centroided = TRUE, smoothed = FALSE,
                      verbose = TRUE, cache = 1)
########################################################################################
### creating checkfiles:
### listing m/z values to be expected for the peptide sequence
expected_fragements<-c(523.942589, 785.410246,
                       518.267073, 776.896971,
                       755.383499, 626.340906,
                       529.288142, 373.187031,
                       302.149918, 205.097154,
                       378.195388, 313.674091,
                       265.147709, 187.097154,
                       151.578597, 103.052215,
                       116.034219, 203.066248,
                       302.134661, 403.18234,
                       559.283451, 687.342029)
MGF_to_IONcheck <- function(Wt) {
  IonCheck<-matrix(nrow = length(Wt), ncol = 31)
  IonCheck_Absolute<-matrix(nrow = length(Wt), ncol = 31)
  colnames(IonCheck)<-c("p3+", "p2+","p3+ -NH2", "p2+ -NH2",
                        "y6", "y5", "y4", "y3", "y2", "y1", 
                        "y6+2", "y5+2", "y4+2", "y3+2", "y2+2", "y1+2",
                        "b1", "b2", "b3", "b4", "b5", "b6", 
                        "SpectrumNr", "#Fragments", "RetentionTime", "Precursor_mz", 
                        "Charge", "Delta_Mass", "Delta_Mass_Oxygen", "Fragment_Sum", 
                        "Fragment_Contribution")
  colnames(IonCheck_Absolute)<-c("p3+", "p2+","p3+ -NH2", "p2+ -NH2",
                                 "y6", "y5", "y4", "y3", "y2", "y1", 
                                 "y6+2", "y5+2", "y4+2", "y3+2", "y2+2", "y1+2",
                                 "b1", "b2", "b3", "b4", "b5", "b6", 
                                 "SpectrumNr", "#Fragments", "RetentionTime", "Precursor_mz", 
                                 "Charge", "Delta_Mass", "Delta_Mass_Oxygen", "Base_Peak", "TIC")
  ########################################################################################
  ### looping
  SpectrumCount<-1
  for (SpectrumCount in c(1:length(Wt))) {
    MSMS_tempp<-Wt[[SpectrumCount]]
    MSMS_temp<-MSMS_tempp@mz
    Intensity_temp<-MSMS_tempp@intensity
    fragment_count<-1
    for (fragment_count in c(1:22)) {
      MSMS_temp_frag<-MSMS_temp-expected_fragements[fragment_count]
      MSMS_temp_frag_squared<-MSMS_temp_frag
      i<-1
      for (i in c(1:length(MSMS_temp))) {
        MSMS_temp_frag_squared[i]<-MSMS_temp_frag[i]*MSMS_temp_frag[i]
        i<-i+1
      }
      ifelse(
        min(MSMS_temp_frag_squared)<=FragTol,
        IonCheck[SpectrumCount,fragment_count]<-Intensity_temp[which(MSMS_temp_frag_squared==min(MSMS_temp_frag_squared))]/max(Intensity_temp),
        IonCheck[SpectrumCount,fragment_count]<-0
      )
      ifelse(
        min(MSMS_temp_frag_squared)<=FragTol,
        IonCheck_Absolute[SpectrumCount,fragment_count]<-Intensity_temp[which(MSMS_temp_frag_squared==min(MSMS_temp_frag_squared))],
        IonCheck_Absolute[SpectrumCount,fragment_count]<-0
      )
      fragment_count<-fragment_count+1
    }
    IonCheck_Absolute[SpectrumCount,30]<-max(Intensity_temp)
    IonCheck_Absolute[SpectrumCount,31]<-sum(Intensity_temp)  
    SpectrumCount<-SpectrumCount+1
  }
  j<-1
  for (j in c(1:length(Wt))) {
    IonCheck[j,23]<-j
    IonCheck_Absolute[j,23]<-j
    IonCheck[j,24]<-sum(IonCheck[j,c(1:22)]!=0)
    IonCheck_Absolute[j,24]<-sum(IonCheck[j,c(1:22)]!=0)
    MSMS_tempp<-Wt[[j]]
    MSMS_temp<-MSMS_tempp@rt
    IonCheck[j,25]<-MSMS_temp
    IonCheck_Absolute[j,25]<-MSMS_temp
    MSMS_temp<-MSMS_tempp@precursorMz
    IonCheck[j,26]<-MSMS_temp
    IonCheck_Absolute[j,26]<-MSMS_temp
    MSMS_temp<-MSMS_tempp@precursorCharge
    IonCheck[j,27]<-MSMS_temp
    IonCheck_Absolute[j,27]<-MSMS_temp
    IonCheck[j,28]<-IonCheck[j,26]*IonCheck[j,27]-IonCheck[j,27]*1.0072764-1568.8059
    IonCheck_Absolute[j,28]<-IonCheck[j,26]*IonCheck[j,27]-IonCheck[j,27]*1.0072764-1568.8059
    IonCheck[j,29]<-IonCheck[j,28]+15.99491463
    IonCheck_Absolute[j,29]<-IonCheck[j,28]+15.99491463
    j<-j+1
  }
  DoubleIonCheck<-list(IonCheck,IonCheck_Absolute)
  return(DoubleIonCheck)
}
########################################################################################
### Running all files through checkup tool
IonCheck_MR1_data<-MGF_to_IONcheck(MR1_data)
########################################################################################
### Apply filters
IonCheck_MR1_data_filtered<-IonCheck_MR1_data[[1]]
IonCheck_MR1_data_filtered<-IonCheck_MR1_data_filtered[IonCheck_MR1_data_filtered[,28]>=Min_mass_shift,]
IonCheck_MR1_data_filtered<-IonCheck_MR1_data_filtered[IonCheck_MR1_data_filtered[,24]>=Min_frag_number,]

IonCheck_MR1_data_filtered_absolute<-IonCheck_MR1_data[[2]]
IonCheck_MR1_data_filtered_absolute<-IonCheck_MR1_data_filtered_absolute[IonCheck_MR1_data_filtered_absolute[,28]>=Min_mass_shift,]
IonCheck_MR1_data_filtered_absolute<-IonCheck_MR1_data_filtered_absolute[IonCheck_MR1_data_filtered_absolute[,24]>=Min_frag_number,]

########################################################################################
### writing required Excel Tables
write.csv(IonCheck_MR1_data_filtered, file = "IonCheck_MR1_data.csv", quote = FALSE, row.names = FALSE, sep = ";")
write.csv(IonCheck_MR1_data_filtered_absolute, file = "IonCheck_MR1_data_absolute.csv", quote = FALSE, row.names = FALSE, sep = ";")
########################################################################################
