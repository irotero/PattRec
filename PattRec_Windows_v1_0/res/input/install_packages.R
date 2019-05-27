        if(!"ShortRead"%in%installed.packages()[,"Package"]){
		source("http://bioconductor.org/biocLite.R")
		biocLite("ShortRead")
	}
