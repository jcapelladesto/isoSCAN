#' Read files and quantify isotopologues in low- or high-resolution data
#' 
#' Read each file and quantify both Area and Maxo of each metabolite isotopologues, 
#' by evaluating the RT window of each metabolite to minimize background noise.
#' Low resolution data requires that \emph{massdiff} and \emph{mzerr} parameters are indicated to process the data.
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
#' @param mzerr Symmetric (+/-) absolute mass error to consider in low-resolution data.
#' @param maxppm Symmetric (+/-) maximum ppm mass error to consider in high-resolution data.
#' @param isotopes Dataframe with stable isotopes information as in  \link{enviPat::isotopes}
#' @param labelatom Labeling atom used as in isotopes$isotope. For example, "13C".
#' @return A data.frame including the Area and Maxo for each file and metabolite isotopologues. If peak was missing, values are NA.
#' @export

autoQ <- function(SampleFiles=NULL, formulaTable=NULL, SNR=3, minscans = 6, RTwin = 5, fit.p = 0.05,
                  resolution = NA, minwidth=1,  maxwidth=4,
                  mzerr = 0.1, massdiff = 1.003355, 
                  maxppm=5, isotopes, labelatom="13C"){
  
  
  if( RTwin > maxwidth){stop("Argument RTwin must not be larger than maxwidth (check Help)",call. = help("autoQ"))}
  if(is.na(resolution)){stop("Please indicate the resolution of the data  as indicated (check Delp)",call. = help("autoQ"))}
  
  if(resolution==1){
    #nominal
    HR <- F
    df <- isoQuant.LR(SampleFiles, formulaTable, SNR, minscans , RTwin, fit.p,
                      mzerror, massdiff, minwidth, maxwidth, HR)
    
  }else{
    HR <-T
    df <- isoQuant.HR(SampleFiles, formulaTable, SNR, minscans, RTwin, fit.p, maxppm,
                      minwidth,  maxwidth, isotopes, labelatom,  resolution, HR)
    
    # remove ppm abundance columns ?
  }
  
  return(df)
  
}