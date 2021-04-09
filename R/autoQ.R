#' Read files and quantify isotopologues in low- or high-resolution data
#' 
#' Read each file and quantify both Area and Maxo of each metabolite isotopologues, 
#' by evaluating the RT window of each metabolite to minimize background noise.
#' Low resolution data requires that \emph{massdiff} and \emph{mzerror} parameters are indicated to process the data.
#' High resolution data requires a limit of ppm (\emph{maxppm}) and that column \emph{Formula} from \emph{formulaTable}
#' includes the molecular formula of the derivatized metabolite for accurate isotopologue mass estimation by \link{enviPat} package.
#' @param SampleFiles List of mzML files to be processed (see \emph{base::list.files()})
#' @param formulaTable Dataframe including CompoundName, Formula, mz, RT and NumAtoms (see \link{buildFormulaTable})
#' @param RTwin Retention time window (in seconds) to extract the metabolite isotopologues
#' @param minwidth Minimum width (in seconds) of a detected EIC peak to be considered
#' @param maxwidth Maximum width (in seconds) of a detected EIC peak to be considered
#' @param minscans Minimum number of recorded scans to use to evaluate a peak.
#' @param fit.p Max linear model p-value (used to determine peak symmetry) for a peak to have its area quantified, otherwise only Maxo is returned.
#' @param SNR Minimum signal-to-noise ratio for EIC peak quantification.
#' @param resolution Set to c(1) in the case of low-resolution data. 
#' Set to a fixed average resolution value if qTOF instrument is used (for example, c(20000)). Set to 'R at MZ' as c(1e5,200) in Orbitrap instruments.
#' @param massdiff Stable isotope mass difference, for C13: 1.003355, necessary for low-resolution. 
#' @param mzerror Symmetric (+/-) absolute mass error to consider in low-resolution data.
#' @param maxppm Symmetric (+/-) maximum ppm mass error to consider in high-resolution data.
#' @param isotopes Dataframe with stable isotopes information as in  \link[enviPat]{isotopes}
#' @param labelatom Labeling atom used as in isotopes$isotope. For example, "13C".
#' @param thr Inherited from \link[enviPat]{isopattern}. Probability below which isotope peaks can be omitted. Default = 0.1. Recommended values between: 0.5 and 0.01. 
#' @return A data.frame including the Area and Maxo for each file and metabolite isotopologues. If peak was missing, values are NA.
#' @details In the case of high-resolution and additional column indicating the estimated mass accuracy error (in ppm) per each isotopologue and sample is included.
#' @export

autoQ <- function(SampleFiles=NULL, formulaTable=NULL, SNR=3, minscans = 6, RTwin = 5, fit.p = 0.05,
                  resolution = NA, minwidth=1,  maxwidth=4, thr=0.1,
                  mzerror = 0.1, massdiff = 1.003355, 
                  maxppm=5, isotopes, labelatom="13C"){
  
  
  if(maxwidth > RTwin){stop("Argument RTwin must be larger than maxwidth (check Help)",call. = help("autoQ"))}
  if(minwidth >= maxwidth){stop("Argument maxwidth must be larger than minwidth (check Help)",call. = help("autoQ"))}
  if(is.na(resolution[1])){stop("Please indicate the resolution of the data  as indicated (check Help)",call. = help("autoQ"))}
  
  if(resolution[1]==1){
    #nominal
    df <- isoQuant.LR(SampleFiles, formulaTable, SNR, minscans , RTwin, fit.p,
                      mzerror, massdiff, minwidth, maxwidth)
    
  }else{
    df <- isoQuant.HR(SampleFiles, formulaTable, SNR, minscans, RTwin, fit.p, maxppm,
                      minwidth,  maxwidth, isotopes, labelatom,  resolution, thr)
    
    # remove ppm abundance columns ?
  }
  df$Isotopologue <- paste0("M+",df$Isotopologue)
  # do not make CpdName factor, only in plot funcs
  return(df)
  
}