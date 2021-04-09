#' @import mzR
#' @import ggplot2
#' @import reshape2
#' @importClassesFrom Rcpp "C++Object"
NULL
#' Reads and builds FormulaTable file into a data.frame
#' 
#' Reads text file and creates \code{NumAtoms} column by obtaining carbon number from compound formula
#' @param x File directory (as in \emph{utils::read.table})
#' @param atom.symbol Atom symbol of the labelling experiment, for C13 set it to \code{"C"}
#' @return formulaTable with NumAtoms column 
#' @export
buildFormulaTable <- 
	function(x=NULL,atom.symbol="C"){
		if(is.matrix(x)){
			stop("Object is not a data.frame but a matrix. Use as.data.frame() to transform it")
		}
		formulaTable <- x
		if (!all(c("CompoundName","Formula","mz","RT")%in%colnames(formulaTable))){
			stop("Column names should be: CompoundName, Formula, mz and RT")
		}
		expr1 <- paste0(atom.symbol,"([[:digit:]]+)")
		expr2 <- paste0(atom.symbol,"1")
		formulaTable$"NumAtoms" <- sapply(formulaTable$Formula, function(x){
			reg <- regmatches(x,gregexpr(expr1,x))[[1]]
			if (length(reg)==0){ reg <- expr2}
			res <- as.numeric(gsub(pattern=atom.symbol,replacement="",reg))
			return(res)
		})
		# if(!is.factor(formulaTable$CompoundName)){
		# formulaTable$CompoundName <- factor(formulaTable$CompoundName,levels=unique(formulaTable$CompoundName))
		# }
		return(formulaTable)
	}