/*
 * Copyright (c) 2016 GRADIANT. All rights reserved.
 * This code cannot be used, copied, modified and/or distributed without the express permission of the authors.
 * This algorithm is protected under a Confidentiality Agreement between the UDTEMC (Fundación Ramón Domínguez) and GRADIANT.
 */
package org.gradiant.JRI;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;
import java.util.concurrent.TimeUnit;

import org.gradiant.utils.Commons;
import org.rosuda.JRI.REXP;
import org.rosuda.JRI.RMainLoopCallbacks;
import org.rosuda.JRI.Rengine;


public class RCNVScript implements RMainLoopCallbacks {
    final private String rCode;
    Boolean print;
    
    public RCNVScript() {
        
        //Use properties file
    	String memoryLimit;
        Properties prop = new Properties();
		InputStream is = null;
		try {
			is = new FileInputStream(Commons.resDir+"configuration.properties");
			prop.load(is);
		} catch(IOException e) {
			System.out.println(e.toString());
		}
		memoryLimit=prop.getProperty("memory.limit");
		
        this.print = false;
        this.rCode =  "scriptR <- function(per_dup, per_del, per_gene_dup, per_gene_del, cov_cont, fvcf, excel, genes, restr, fixed, down, do_plot, userDB, passwordDB, wd, name_test, name_bed, name_fasta, name_cont, polymorphic, poly_regions) {\n" +
        	        "setwd(wd)\n" +
        		"options(warn= -1)\n" +
        		"memory.limit("+Integer.parseInt(memoryLimit)+")\n"+
        		"set.seed(1)\n" +
        	        "list.of.packages <- c(\"outliers\", \"openxlsx\", \"RMySQL\", \"Rcpp\", \"remotes\")\n" +
        	        "new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,\"Package\"])]\n" +
        	        "if(length(new.packages)!=0) install.packages(new.packages, repos = \"http://cloud.r-project.org/\")\n" +
					"if(!\"poolR\"%in%installed.packages()[,\"Package\"]) remotes::install_github(\"ozancinar/poolR\")\n" +
        	        "library(outliers)\n" +
        	        "library(RMySQL)\n" +
        	        "library(openxlsx)\n" +
					"library(poolR)\n" +
			"if(!\"ShortRead\"%in%installed.packages()[,\"Package\"]){\n" +
			"	source(\"https://bioconductor.org/biocLite.R\")\n" +
			"	biocLite(\"ShortRead\", ask=FALSE)\n" +
			"}\n" +
        		"library(ShortRead)\n" +
        	        "#------------------------------------------------------------------------------\n" +
        	        "per_dup <- as.numeric(per_dup)\n" +
        	        "per_del <- as.numeric(per_del)\n" +
        	        "per_gene_dup <- as.numeric(per_gene_dup)\n" +
        	        "per_gene_del <- as.numeric(per_gene_del)\n" +
        	        "cov_cont <- as.numeric(cov_cont)\n" +
        	        "excel <- as.logical(excel)\n" +
        	        "genes <- as.logical(genes)\n" +
        	        "restr <- as.logical(restr)\n" +
        	        "fixed <- as.logical(fixed)\n" +
        	        "down <- as.logical(down)\n" +
        	        "do_plot <- as.logical(do_plot)\n" +
        	        "#------------------------------------------------------------------------------\n" +
        		"ppal <- function(fvcf, per, per_gene, cov_cont, name_test, name_cont, name_bed, name_fasta, excel, genes, polymorphic, poly_regions, do_plot, restr, fixed, down){\n" +
        		"	lectura <- rd(name_test, name_cont, name_bed, name_fasta, fvcf, polymorphic, poly_regions, fixed)\n" +
        		"	name_cont.old <- name_cont\n" +
        		"	name_test.old <- name_test\n" +
        		"	name_cont <- lectura$name_cont\n" +
        		"	name_test <- lectura$name_test\n" +
        		"	listgenes <- names(table(lectura$test[,4]))\n" +
        		"	sex <- gender(lectura$test, lectura$cont, lectura$n)\n" +
       			"	output_genes <- c()\n" +
        		"	if(genes==T) output_genes <- Rgenes(per_gene, lectura$n, listgenes, lectura$test, lectura$cont, lectura$gc, sex)\n" +
        		"	fun_mxs <- rgns(listgenes, lectura$n, lectura$test, lectura$cont, output_genes, cov_cont, genes)\n" +
        		"	mxs <- as.data.frame(fun_mxs$mxs, stringsAsFactors=F, row.names=F)\n" +
        		"	output_genes <- fun_mxs$output_genes\n" +
        		"	normalized <- nrml(listgenes, lectura$n, lectura$test, lectura$cont, mxs)\n" +
        		"	std <- sdgenes(listgenes, lectura$n, normalized$test, normalized$cont)\n" +
        		"	output <- tabla(lectura$n, listgenes, normalized$test, normalized$cont, lectura$gc, mxs, std, cov_cont)\n" +
        		"	if(fvcf!=\"NO\") output <- vcf_filter(output, lectura$vcf, cov_cont)\n" +
        		"	if(nrow(output)!=0) output <- concat(output,cov_cont,per)\n" +
        		"	nexons <- norm_exon(lectura$n, normalized$test, normalized$cont, listgenes, output)\n" +
        		"	if(nrow(output)!=0) residuos <- regression(lectura$n, nexons)\n" +
        		"	if(nrow(output)!=0) output1 <- filt(output, per, normalized$cont, lectura$n, residuos, sex, restr, cov_cont)\n" +
        		"	if(nrow(output)==0) output1 <- output\n" +
        		"	if(down==TRUE && nrow(output1)>2) output1 <- downsampling(name_test, name_cont, lectura$n, lectura$bed, output1, cov_cont, lectura$gc, per, sex, restr)\n" +
        		"	output <- output1\n" +
        		"	if(nrow(output1)!=0) output <- freq(output1, userDB, passwordDB)\n" +
        		"	if(nrow(output)>0){\n" +
        		"		if(lectura$n>1) output <- DCVar_exon(output, lectura$n, nexons)\n" +
        		"		output <- annot(output, userDB, passwordDB)\n" +
        		"	}\n" +
        		"	if(fvcf!=\"NO\" && class(output_genes)!=\"NULL\" && nrow(output_genes)>0) output_genes <- vcf_filter(output_genes, lectura$vcf, cov_cont)\n" +
        		"	info <- information(lectura$n, normalized$test, lectura$cont, name_test, name_cont, sex, nexons)\n" +
        		"	datos_dcvar <- as.numeric(info[(which(is.na(info[,1]))+1):(which(is.na(info[,1]))+2),2])\n" +
        		"	names(output) <- gsub(\"gen\",\"gene\",names(output))\n" +
        		"	names(output) <- gsub(\"Gen\",\"Gene\",names(output))\n" +
        		"	if(class(output_genes)!=\"NULL\") names(output_genes) <- gsub(\"Gen\",\"Gene\",names(output_genes))\n" +
        		"	if(excel==F) wrt(per, output, output_genes, name_test, name_cont, genes, info)\n" +
        		"	if(excel==T) write.excel(output, name_test, name_cont, datos_dcvar, per, output_genes, genes, info)\n" +
        		"	if(do_plot==T) graphs(output, name_test, name_cont, lectura$n, lectura$test, lectura$cont, mxs)\n" +
        		"}\n" +
        	        "#------------------------------------------------------------------------------\n" +
			"rd <- function(name_test, name_cont, name_bed, name_fasta, fvcf, polymorphic, poly_regions, fixed){\n" +
			"	bed <- read.table(name_bed, header=F, sep=\"\t\")\n" +
			"	bed[,1] <- as.character(bed[,1])\n" +
			"	bed[,4] <- as.character(bed[,4])\n" +
			"	cr <- names(table(bed[,1]))\n" +
			"	n <- length(name_cont)\n" +
			"	fasta.file <- scanFa(name_fasta, GRanges(seqname=bed[,1], IRanges(start=bed[,2]+1, end=bed[,3])))\n" +
			"	gc <- cbind(bed, letterFrequency(fasta.file,\"GC\")/letterFrequency(fasta.file,\"ATGC\"))\n" +
			"	names(gc) <- NULL\n" +
			"	gc[,1] <- toupper(gc[,1])\n" +
			"	gc[,1] <- gsub(\"CHR\",\"chr\",gc[,1])\n" +
			"	if(length(grep(\"chr\",gc[,1]))==0) gc[,1] <- paste0(\"chr\",gc[,1])\n" +
			"	gc[gc[,5]%in%c(Inf, -Inf),5] <- NA\n" +
			"	dat <- coverage(readGAlignments(file.path(name_test)))\n" +
			"	test <- c()\n" +
			"	for(crom in cr){\n" +
			"		dat1 <- as.numeric(dat[[which(names(dat)==crom)]])\n" +
			"		bed1 <- bed[bed[,1]==crom,]\n" +
			"		for(j in 1:nrow(bed1)){\n" +
			"			auxx <- as.data.frame(matrix(ncol=6,nrow=(bed1[j,3]-bed1[j,2])))\n" +
			"			auxx[,1] <- bed1[j,1]\n" +
			"			auxx[,4] <- bed1[j,4]\n" +
			"			auxx[,2] <- bed1[j,2]\n" +
			"			auxx[,3] <- bed1[j,3]\n" +
			"			auxx[,5] <- 1:nrow(auxx)\n" +
			"			auxx[,6] <- dat1[auxx[,2]+(1:nrow(auxx))]\n" +
			"			test <- rbind(test,auxx)\n" +
			"		}\n" +
			"	}\n" +
			"	names(test) <- c(\"chr\",\"start\",\"end\",\"gen\",\"pb\",\"cov\")\n" +
			"	test$gen <- toupper(test$gen)\n" +
			"	test$chr <- toupper(test$chr)\n" +
			"	test$chr <- gsub(\"CHR\",\"chr\",test$chr)\n" +
			"	if(length(grep(\"chr\",test$chr))==0) test$chr <- paste0(\"chr\",test$chr)\n" +
			"	cname <- file.path(name_cont)\n" +
			"	aux.cont <- as.data.frame(matrix(ncol=n,nrow=nrow(test)))\n" +
			"	for(k in 1:n){\n" +
			"		bm <- readGAlignments(cname[k])\n" +
			"		dat <- coverage(bm)\n" +
			"		cont <- c()\n" +
			"		for(crom in cr){\n" +
			"			dat1 <- as.numeric(dat[[which(names(dat)==crom)]])\n" +
			"			bed1 <- bed[bed[,1]==crom,]\n" +
			"			for(j in 1:nrow(bed1)){\n" +
			"				auxx <- as.data.frame(matrix(ncol=6,nrow=(bed1[j,3]-bed1[j,2])))\n" +
			"				auxx[,1] <- bed1[j,1]\n" +
			"				auxx[,4] <- bed1[j,4]\n" +
			"				auxx[,2] <- bed1[j,2]\n" +
			"				auxx[,3] <- bed1[j,3]\n" +
			"				auxx[,5] <- 1:nrow(auxx)\n" +
			"				auxx[,6] <- dat1[auxx[,2]+(1:nrow(auxx))]\n" +
			"				cont <- rbind(cont,auxx)\n" +
			"			}\n" +
			"		}\n" +
			"		names(cont) <- c(\"chr\",\"start\",\"end\",\"gen\",\"pb\",\"cov\")\n" +
			"		aux.cont[,k] <- cont[,6]\n" +
			"	}\n" +
			"	if(n>=2 & fixed==F){\n" +
			"		fnx <- function(x){cor(test$cov,x)}\n" +
			"		aux.cor <- sapply(aux.cont, \"fnx\")\n" +
			"		while(n>=2 & min(aux.cor)<0.95){\n" +
			"			l1 <- which.min(aux.cor)\n" +
			"			if(aux.cor[l1]<0.95){\n" +
			"				name_cont <- name_cont[-l1]\n" +
			"				n <- length(name_cont)\n" +
			"				aux.cont <- aux.cont[,-l1]\n" +
			"				aux.cor <- aux.cor[-l1]\n" +
			"			}\n" +
			"		}\n" +
			"	}\n" +
			"	cont <- as.data.frame(matrix(ncol=4+n,nrow=nrow(test)))\n" +
			"	cont[,1:4] <- test[,1:4]\n" +
			"	cont[,5:(4+n)] <- aux.cont[,1:n]\n" +
			"	name_cont <- gsub(\"\\\\\\\\\", \"/\", name_cont)\n" +
			"	name_test <- gsub(\"\\\\\\\\\", \"/\", name_test)\n" +
			"	fn.name <- function(x){gsub(\".bam\", \"\", unlist(strsplit(x,\"/\"))[length(unlist(strsplit(x,\"/\")))])}\n" +
			"	name_cont <- unlist(lapply(name_cont, fn.name))\n" +
			"	name_test <- unlist(lapply(name_test, fn.name))\n" +
			"	if(n>1){\n" +
			"		cont[,(5+n)] <- apply(cont[,5:(4+n)],1,mean)\n" +
			"		names(cont) <- c(\"chr\", \"start\", \"end\", \"gen\", name_cont, \"mean\")\n" +
			"	}\n" +
			"	if(n==1){\n" +
			"		cont$mean <- cont[,5]\n" +
			"		names(cont)[1:5] <- c(\"chr\", \"start\", \"end\", \"gen\", name_cont)\n" +
			"	}\n" +
			"	cont$gen <- toupper(cont$gen)\n" +
			"	cont$chr <- toupper(cont$chr)\n" +
			"	cont$chr <- gsub(\"CHR\",\"chr\",cont$chr)\n" +
			"	if(length(grep(\"chr\",cont$chr))==0) cont$chr <- paste0(\"chr\",cont$chr)\n" +
			"	if(polymorphic==\"SI\"){\n" +
			"		plm <- read.table(poly_regions, header=F, sep=\"\t\", stringsAsFactors=F)\n" +
			"		plm[,1] <- toupper(plm[,1])\n" +
			"		plm[,1] <- gsub(\"CHR\",\"chr\",plm[,1])\n" +
			"		if(length(grep(\"chr\",plm[,1]))==0) plm[,1]<-paste0(\"chr\",plm[,1])\n" +
			"		for(j in 1:nrow(plm)){\n" +
			"			m <- which(test$chr==plm[j,1] & test$end>=plm[j,2] & test$start<=plm[j,3])\n" +
			"			if(length(m)>0) test <- test[-m,]\n" +
			"			m <- which(cont$chr==plm[j,1] & cont$end>=plm[j,2] & cont$start<=plm[j,3])\n" +
			"			if(length(m)>0) cont <- cont[-m,]\n" +
			"			m <- which(gc[,1]==plm[j,1] & gc[,3]>=plm[j,2] & gc[,2]<=plm[j,3])\n" +
			"			if(length(m)>0) gc <- gc[-m,]\n" +
			"		}\n" +
			"	}\n" +
			"	if(fvcf!=\"NO\"){\n" +
			"		name_vcf <- dir()[grepl(\".vcf\",dir())]\n"+
			"		if(length(name_vcf)==0) vcf <- as.data.frame(matrix(ncol=10, nrow=0))\n" +
			"		if(length(name_vcf)!=0) vcf <- read.table(name_vcf, header=F, sep=\"\t\", stringsAsFactors=F)\n" +
			"		names(vcf) <- c(\"CHROM\", \"POS\", \"ID\", \"REF\", \"ALT\", \"QUAL\", \"FILTER\", \"INFO\", \"FORMAT\", \"SAMPLE\")\n" +
			"		vcf <- vcf[vcf$QUAL>=200,]\n" +
			"		vcf <- vcf[grepl(\"0/1\",vcf$SAMPLE),]\n" +
			"		vcf <- vcf[which(grepl(\"INDEL\",vcf$INFO)==F),]\n" +
			"		vcf$CHROM <- toupper(vcf$CHROM)\n" +
			"		vcf$CHROM <- gsub(\"CHR\",\"chr\",vcf$CHROM)\n" +
			"		return(list(\"n\"=n, \"test\"=test, \"cont\"=cont, \"gc\"=gc, \"vcf\"=vcf, \"name_cont\"=name_cont, \"name_test\"=name_test, \"bed\"=bed))\n" +
			"	}\n" +
			"	if(fvcf==\"NO\") return(list(\"n\"=n, \"test\"=test, \"cont\"=cont, \"gc\"=gc, \"name_cont\"=name_cont, \"name_test\"=name_test, \"bed\"=bed))\n" +
			"}\n" +
        	        "#------------------------------------------------------------------------------\n" +
        		"gender <- function(test,cont,n){\n" +
        		"	sex_test <- ifelse(mean(test[grepl(\"X\",test$chr),6]) <= 0.6*mean(test[which(grepl(\"X\",test$chr)==F),6]),\"M\",\"F\")\n" +
        		"	aux <- as.data.frame(matrix(ncol=n,nrow=nrow(test)))\n" +
        		"	sex_cont <- c()\n" +
        		"	for(i in 1:n){\n" +
        		"		cnt <- cont[,c(1,4+i)]\n" +
        		"		sex_cont <- c(sex_cont, ifelse(mean(cnt[grepl(\"X\",cnt$chr),2]) <= 0.6*mean(cnt[which(grepl(\"X\",cnt$chr)==F),2]), \"M\", \"F\"))\n" +
        		"	}\n" +
        		"	sex <- c(sex_test,sex_cont)\n" +
        		"	return(sex)\n" +
        		"}\n" +
        	        "#------------------------------------------------------------------------------\n" +
        		"Rgenes <- function(per, n, listgenes, test, cont, gc, sex){\n" +
        		"	output <- c()\n" +
        		"	out_genes <- as.data.frame(matrix(ncol=12,nrow=length(listgenes)))\n" +
        		"	names(out_genes) <- c(\"Chr\", \"Start\", \"End\", \"Gen\", \"pvalue\", \"Type\", \"Test_mean\", \"Cont_mean\", \"Cont_sd\", \"CV\", \"GC\", \"% decrease/increase\")\n" +
        		"	test[which(test==0,arr.ind=T)] <- 1.e-20\n" +
        		"	cont[which(cont==0,arr.ind=T)] <- 1.e-20\n" +
        		"	normtest <- test\n" +
        		"	if(sex[1]==\"F\" || is.na(sex[1])) normtest$cov <- normtest$cov/mean(normtest$cov,na.rm=T)\n" +
        		"	if(sex[1]==\"M\" && !is.na(sex[1])){\n" +
        		"		x <- which(grepl(\"X\",normtest$chr))\n" +
        		"		normtest$cov[x] <- normtest$cov[x]/mean(normtest$cov[x], na.rm=T)\n" +
        		"		normtest$cov[-x] <- normtest$cov[-x]/mean(normtest$cov[-x], na.rm=T)\n" +
        		"	}\n" +
        		"	normcont <- as.data.frame(cbind(cont$chr, cont$start, cont$end, cont$gen, cont$mean, cont[,5:(ncol(cont)-1)]))\n" +
        		"	names(normcont)[1:5] <- c(\"chr\",\"start\",\"end\",\"gen\",\"mean\")\n" +
        		"	normcont$mean <- as.numeric(as.character(normcont$mean))\n" +
        		"	if(n>1){\n" +
        		"		for(i in 2:(n+1)){\n" +
        		"			if(sex[i]==\"F\" || is.na(sex[i])) normcont[,4+i] <- normcont[,4+i]/mean(normcont[,4+i],na.rm=T)\n" +
        		"			if(sex[i]==\"M\" && !is.na(sex[i])){\n" +
        		"				x <- which(grepl(\"X\",normcont$chr))\n" +
        		"				normcont[x,4+i] <- normcont[x,4+i]/mean(normcont[x,4+i], na.rm=T)\n" +
        		"				normcont[-x,4+i] <- normcont[-x,4+i]/mean(normcont[-x,4+i], na.rm=T)\n" +
        		"			}\n" +
        		"		}\n" +
        		"		normcont[,5] <- rowSums(normcont[,6:(5+n)])/n\n" +
        		"	}\n" +
        		"	if(n==1){\n" +
        		"		if(sex[2]==\"F\" || is.na(sex[2])) normcont$mean <- normcont$mean/mean(normcont$mean,na.rm=T)\n" +
        		"		if(sex[2]==\"M\" && !is.na(sex[2])){\n" +
        		"			x <- which(grepl(\"X\",normcont$chr))\n" +
        		"			normcont$mean[x] <- normcont$mean[x]/mean(normcont$mean[x], na.rm=T)\n" +
        		"			normcont$mean[-x] <- normcont$mean[-x]/mean(normcont$mean[-x], na.rm=T)\n" +
        		"		}\n" +
        		"	}\n" +
        		"	for(k in 1:length(listgenes)){\n" +
        		"		auxtest <- test[test[,4]==listgenes[k],]\n" +
        		"		auxcont <- cont[cont[,4]==listgenes[k],]\n" +
        		"		regions <- as.numeric(names(table(auxtest[,2])))\n" +
        		"		mxs <- c()\n" +
        		"		for(j in regions){\n" +
        		"			cr <- auxcont[auxcont$start==j,]\n" +
        		"			tr <- auxtest[auxtest$start==j,]\n" +
        		"			mn <- 0.05*sum(cr$mean)\n" +
        		"			mx <- 0.95*sum(cr$mean)\n" +
        		"			a <- cumsum(cr$mean)\n" +
        		"			cr <- cr[which(a>mn & a<mx),]\n" +
        		"			tr <- tr[which(a>mn & a<mx),]\n" +
        		"			if(n>1){\n" +
        		"				mmm <- apply(cr[cr[,2]==j,6:(5+n)], FUN=\"max\", MARGIN=2)\n" +
        		"				mmm <- mean(mmm,na.rm=T)\n" +
        		"			}\n" +
        		"			if(n==1) mmm <- max(cr[cr[,2]==j,5], na.rm=T)\n" +
        		"			aux <- cbind(max(tr[tr[,2]==j,6], na.rm=T), mmm)\n" +
        		"			mxs <- rbind(mxs,aux)\n" +
        		"		}\n" +
        		"		out_genes$Chr[k] <- auxtest$chr[1]\n" +
        		"		out_genes$Start[k] <- min(auxtest$start,na.rm=T)\n" +
        		"		out_genes$End[k] <- max(auxtest$end,na.rm=T)\n" +
        		"		out_genes$Gen[k] <- listgenes[k]\n" +
        		"		out_genes$Test_mean[k] <- mean(mxs[,1])\n" +
        		"		out_genes$Cont_mean[k] <- mean(mxs[,2])\n" +
        		"		if(n>1){\n" +
        		"			out_genes$Cont_sd[k] <- sd(c(colSums(auxcont[,6:(5+n)])/nrow(auxcont)))\n" +
        		"			if(out_genes$Cont_mean[k]!=0) out_genes$CV[k] <- out_genes$Cont_sd[k]/out_genes$Cont_mean[k]\n" +
        		"			if(out_genes$Cont_mean[k]==0) out_genes$CV[k] <- 1\n" +
        		"		}\n" +
        		"		if(n==1){\n" +
        		"			out_genes$Cont_sd[k] <- 0\n" +
        		"			out_genes$CV[k] <- 0\n" +
        		"		}\n" +
        		"		auxtest <- normtest[normtest[,4]==listgenes[k],]\n" +
        		"		auxcont <- normcont[normcont[,4]==listgenes[k],]\n" +
        		"		gen_aux <- log(auxtest[,6]/auxcont[,5])\n" +
        		"		med_gen <- mean(gen_aux, na.rm=T)\n" +
        		"		std_gen <- sd(gen_aux, na.rm=T)\n" +
        		"		st_t_gen <- (med_gen/std_gen)\n" +
        		"		out_genes$pvalue[k] <- 2*pnorm(abs(st_t_gen), lower.tail=F)\n" +
        		"		x1 <- max(mean(auxtest[,6],na.rm=T),mean(auxcont[,5],na.rm=T))\n" +
        		"		x2 <- min(mean(auxtest[,6],na.rm=T),mean(auxcont[,5],na.rm=T))\n" +
        		"		out_genes[k,ncol(out_genes)] <- 1 - x2/x1\n" +
        		"		out_genes$Type[k] <- \"deletion\"\n" +
        		"		if(max(mxs[,1],na.rm=T)>max(mxs[,2],na.rm=T)) out_genes$Type[k] <- \"duplication\"\n" +
        		"		out_genes[k,ncol(out_genes)] <- round(100*out_genes[k,ncol(out_genes)],2)\n" +
        		" 		out_genes$GC[k] <- mean(gc[gc[,4]==listgenes[k],5])\n" +
        		"	}\n" +
        		"	out_genes <- out_genes[order(out_genes$pvalue),]\n" +
        		"	if(per[1]>=0.3 && per[2]>=0.26){\n" +
        		"		m <- c(which(out_genes[,ncol(out_genes)]>100*0.3 & out_genes$Type==\"deletion\"),which(out_genes[, ncol(out_genes)]>100*0.26 & out_genes$Type==\"duplication\"))\n" +
        		"		out_genes <- out_genes[m,]\n" +
        		"	}\n" +
        		"	ln <- length(which(out_genes$pvalue<0.05))\n" +
        		"	if(ln>15 && ln<50) out_genes$pvalue <- out_genes$pvalue*nrow(out_genes)/c(1:nrow(out_genes))\n" +
        		"	if(ln>=50) out_genes$pvalue <- out_genes$pvalue*nrow(out_genes)\n" +
        		"	if(nrow(out_genes)>0) out_genes <- out_genes[out_genes$pvalue<0.05,]\n" +
        		"	m <- c(which(out_genes[,ncol(out_genes)]>100*per[1] & out_genes$Type==\"deletion\"),which(out_genes[, ncol(out_genes)]>100*per[2] & out_genes$Type==\"duplication\"))\n" +
        		"	out_genes <- out_genes[m,]\n" +
        		"	m <- which(out_genes$Cont_mean>=cov_cont)\n" +
        		"	out_genes <- out_genes[m,]\n" +
        		"	out_genes[,c(7:11)] <- round(out_genes[,c(7:11)],2)\n" +
        		"	return(out_genes)\n" +
        		"}\n" +
        	        "#------------------------------------------------------------------------------\n" +
        		"rgns <- function(listgenes, n, test, cont, output_genes, cov_cont, genes){\n" +
        		"	listmx <- c()\n" +
        		"	for(k in 1:length(listgenes)){\n" +
        		"		auxtest <- test[test[,4]==listgenes[k],]\n" +
        		"		auxtest[which(auxtest==0,arr.ind=T)] <- 1.e-20\n" +
        		"		auxcont <- cont[cont[,4]==listgenes[k],]\n" +
        		"		auxcont[which(auxcont==0,arr.ind=T)] <- 1.e-20\n" +
        		"		regions <- as.numeric(names(table(auxtest[,2])))\n" +
        		"		exons_test <- c()\n" +
        		"		exons_cont <- c()\n" +
        		"		for(i in 1:length(regions)){\n" +
        		"			exons_test <- rbind(exons_test, c(regions[i],max(auxtest$cov[auxtest$start==regions[i]])))\n" +
        		"			if(n==1) exons_cont <- rbind(exons_cont, c(regions[i],max(auxcont[auxcont$start==regions[i],5])))\n" +
        		"			if(n>1) exons_cont <- rbind(exons_cont, c(regions[i],apply(auxcont[auxcont$start==regions[i],5:(4+n)],2,max)))\n" +
        		"		}\n" +
        		"		it <- which(exons_test[,2]<cov_cont)\n" +
        		"		if(n==1) ic <- which(exons_cont[,2]<cov_cont)\n" +
        		"		if(n>1 && nrow(exons_cont)>1) ic <- which(apply(exons_cont[,2:(n+1)], FUN=\"mean\", MARGIN=1)<cov_cont)\n" +
        		"		if(n>1 && nrow(exons_cont)==1) ic <- which(mean(exons_cont[,2:(n+1)])<cov_cont)\n" +
        		"		itc <- union(it,ic)\n" +
        		"		if(length(itc)>0 && length(itc)!=length(regions)){\n" +
        		"			exons_test <- exons_test[-itc,]\n" +
        		"			exons_cont <- exons_cont[-itc,]\n" +
        		"			regions <- regions[-itc]\n" +
        		"		}\n" +
        		"		if(length(ic)>0 && length(itc)==length(regions) && length(ic)!=length(regions)){\n" +
        		"			exons_test <- exons_test[-ic,]\n" +
        		"			exons_cont <- exons_cont[-ic,]\n" +
        		"			regions <- regions[-ic]\n" +
        		"		}\n" +
        		"		if(class(exons_test)!=\"matrix\") exons_test <- as.matrix(t(exons_test))\n" +
        		"		if(class(exons_cont)!=\"matrix\") exons_cont <- as.matrix(t(exons_cont))\n" +
        		"		exons_cont1 <- c()\n" +
        		"		nrw <- nrow(exons_cont)\n" +
        		"		if(n==1) exons_cont1 <- exons_cont\n" +
        		"		if(n>1 && nrw>1) exons_cont1 <- cbind(exons_cont[,1],apply(exons_cont[,2:(n+1)],1,mean))\n" +
        		"		if(n>1 && nrw<=1) exons_cont1 <- cbind(exons_cont[,1],mean(exons_cont[,2:(n+1)]))\n" +
        		"		dg <- c()\n" +
        		"		if(genes==T){\n" +
        		"			if(nrow(output_genes)>0){\n" +
        		"				dg <- which(output_genes$Gen==listgenes[k] & output_genes$Type==\"deletion\")\n" +
        		"				if(length(dg)>0){\n" +
        		"					ratio_mxs <- exons_test[,2]/exons_cont1[,2]\n" +
        		"					n_mxs <- which(ratio_mxs>=0.8)\n" +
        		"					if(length(n_mxs)>0){\n" +
        		"						if(length(n_mxs)>1) n_mxs <- which.max(ratio_mxs)\n" +
        		"						mxtest <- exons_test[n_mxs,2]\n" +
        		"						if(n==1) mxcont <- exons_cont1[n_mxs,2]\n" +
        		"						if(n>1) mxcont <- exons_cont[n_mxs,2:(n+1)]\n" +
        		"						output_genes <- output_genes[-dg,]\n" +
        		"					}\n" +
        		"					if(length(n_mxs)==0) dg <- c()\n" +
        		"				}\n" +
        		"			}\n" +
        		"		}\n" +
        		"		if(genes==F || length(dg)==0){\n" +
        		"			maximos <- cvp(n, exons_test, exons_cont1, regions)\n" +
        		"			mxtest <- maximos[1]\n" +
        		"			if(n==1) mxcont <- maximos[2]\n" +
        		"			if(n>1 && class(exons_cont)==\"matrix\"){\n" +
        		"				nnn <- which(exons_cont1[,2]==maximos[2])\n" +
        		"				if(length(nnn)==1) mxcont <- exons_cont[exons_cont[,1]==exons_cont1[nnn,1],2:(n+1)]\n" +
        		"				if(length(nnn)>1) mxcont <- apply(exons_cont[exons_cont[,1]%in%exons_cont1[nnn,1],2:(n+1)],FUN=\"max\",MARGIN=2)\n" +
        		"			}\n" +
        		"			if(n>1 && class(exons_cont)!=\"matrix\") mxcont <- exons_cont[2:(n+1)]\n" +
        		"		}\n" +
        		"		listmx <- rbind(listmx,cbind(listgenes[k],mxtest,rbind(mxcont)))\n" +
        		"	}\n" +
        		"	return(list(\"mxs\"=listmx,\"output_genes\"=output_genes))\n" +
        		"}\n" +
        	        "#------------------------------------------------------------------------------\n" +
        		"cvp <- function(n, exons_test, exons_cont,regions){\n" +
        		"	regions1 <- regions\n" +
        		"	ratio_mxs <- exons_test[,2]/exons_cont[,2]\n" +
        		"	tt <- c()\n" +
        		"	if(length(regions1)>=6 && length(ratio_mxs)>=6 && length(table(ratio_mxs))>1){\n" +
        		"		x <- matrix(ncol=5,nrow=length(ratio_mxs))\n" +
        		"		for(i in 1:5){\n" +
        		"			x[,i] <- as.character(kmeans(ratio_mxs,2)$cluster)\n" +
        		"			x[x[,i]==x[1,i],i] <- \"g1\"\n" +
        		"			x[x[,i]!=x[1,i],i] <- \"g2\"\n" +
        		"		}\n" +
        		"		a <- p <- 1\n" +
        		"		x1 <- x2 <- c()\n" +
        		"		if(class(apply(x, FUN=\"table\", MARGIN=1))==\"list\") a <- 0\n" +
        		"		if(a!=0){\n" +
        		"			x1 <- ratio_mxs[which(x[,1]==\"g1\")]\n" +
        		"			x2 <- ratio_mxs[which(x[,1]== \"g2\")]\n" +
        		"			p <- wilcox.test(x1,x2)$p.value\n" +
        		"			std <- sd(ratio_mxs)\n" +
        		"			if(p<0.025 && !is.na(p)){\n" +
        		"				ifelse(abs(mean(x1,na.rm=T)-1)<abs(mean(x2,na.rm=T)-1), tt <- x1, tt <- x2)\n" +
        		"			}\n" +
        		"			if(p>=0.025 && !is.na(p) && (length(x1)==2 || length(x2)==2)){\n" +
        		"				if(length(x2)==2){\n" +
        		"					tst <- c(chisq.out.test(c(x1,x2[1]),variance=var(c(x1,x2[1])))$p.value, chisq.out.test(c(x1,x2[2]),variance=var(c(x1,x2[2])))$p.value)\n" +
        		"					if(tst[1]<0.05 && tst[2]<0.05) ifelse(abs(mean(x1,na.rm=T)-1) < abs(mean(x2,na.rm=T)-1), tt <- x1, tt <- x2)\n" +
        		"				}\n" +
        		"				if(length(x1)==2){\n" +
        		"					tst <- c(chisq.out.test(c(x2,x1[1]),variance=var(c(x2,x1[1])))$p.value, chisq.out.test(c(x2,x1[2]),variance=var(c(x2,x1[2])))$p.value)\n" +
        		"					if(tst[1]<0.05 && tst[2]<0.05) ifelse(abs(mean(x1,na.rm=T)-1) < abs(mean(x2,na.rm=T)-1), tt <- x1, tt <- x2)\n" +
        		"				}\n" +
        		"			}\n" +
        		"		}\n" +
        		"		if(a==0 || (p>=0.025 && length(x1)!=2 && length(x2)!=2)  || is.na(p)){\n" +
        		"			xx <- c()\n" +
        		"			tt <- ratio_mxs\n" +
        		"			tst <- chisq.out.test(tt,variance=var(tt))\n" +
        		"			while(tst$p.value<0.05 && length(tt)>=3 & !is.na(tst$p.value)){\n" +
        		"				xx <- c(xx,which(ratio_mxs==outlier(tt)))\n" +
        		"				tt <- ratio_mxs[-xx]\n" +
        		"				tst <- chisq.out.test(tt,variance=var(tt))\n" +
        		"			}\n" +
        		"		}\n" +
        		"	}\n" +
        		"	if(length(regions1)<6 && length(regions1)>3 && length(ratio_mxs)<6 && length(ratio_mxs)>3 && length(table(ratio_mxs))>1){\n" +
        		"		ratio_mxs1 <- rm.outlier(ratio_mxs)\n" +
        		"		ratio_mxs2 <- outlier(ratio_mxs)\n" +
        		"		if(length(ratio_mxs1)>0){\n" +
        		"			if(length(rm.outlier(ratio_mxs1)>0) && abs(outlier(ratio_mxs1)-mean(rm.outlier(ratio_mxs1)))>abs(outlier(ratio_mxs1)-mean(ratio_mxs2))){\n" +
        		"				ratio_mxs2 <- c(ratio_mxs2, outlier(ratio_mxs1))\n" +
        		"				ratio_mxs1 <- rm.outlier(ratio_mxs1)\n" +
        		"				if(length(rm.outlier(ratio_mxs1)>0) && abs(outlier(ratio_mxs1)-mean(rm.outlier(ratio_mxs1)))>abs(outlier(ratio_mxs1)-mean(ratio_mxs2))){\n" +
        		"					ratio_mxs2 <- c(ratio_mxs2, outlier(ratio_mxs1))\n" +
        		"					ratio_mxs1 <- rm.outlier(ratio_mxs1)\n" +
        		"				}\n" +
        		"			}\n" +
        		"		        if(mean(ratio_mxs2)>(1.5*mean(ratio_mxs1,na.rm=T)) || mean(ratio_mxs2)<(mean(ratio_mxs1,na.rm=T)/2)) tt <- ratio_mxs[which.min(abs(1-ratio_mxs))]\n" +
        	        "		}\n" +
        		"	}\n" +
        		"	if(length(tt)>0) regions1 <- regions[which(ratio_mxs%in%tt)]\n" +
        		"	mxs <- c(max(exons_test[exons_test[,1]%in%regions1,2]), max(exons_cont[exons_cont[,1]%in%regions1,2]))\n" +
        		"	return(mxs)\n" +
        		"}\n" +
        	        "#------------------------------------------------------------------------------\n" +
        		"nrml <- function(listgenes, n, test, cont, mxs){\n" +
        		"	test$normcov <- NA\n" +
        		"	if(n>1) x <- apply(cont[,5:(n+4)],1,\"sd\")\n" +
        		"	if(n==1) x <- rep(0,nrow(cont))\n" +
        		"	for(i in listgenes){\n" +
        		"		k <- which(mxs[,1]==i)\n" +
        		"		m <- which(test[,4]==i)\n" +
        		"		test[m,7] <- test[m,6]/as.numeric(mxs[k,2])\n" +
        		"		m <- which(cont[,4]==i)\n" +
        		"		if(n==1) cont[m,5] <- cont[m,5]/as.numeric(mxs[k,3])\n" +
        		"		if(n>1) for(j in 3:(n+2)) cont[m,2+j] <- cont[m,2+j]/as.numeric(mxs[k,j])\n" +
        		"		if(length(k)!=0) mxs <- mxs[-k, ]\n" +
        		"	}\n" +
        		"	if(n==1) names(cont)[5]<-\"normmean\"\n" +
        		"	if(n>1) cont$normmean<-apply(cont[,5:(n+4)],1,\"mean\")\n" +
        		"	cont$covsd <- x\n" +
        		"	cont$covmin <- cont$covmax <- 0\n" +
        		"	if(n>1){\n" +
        		"		cont$covmin <- apply(cont[,5:(n+4)],1,\"min\")\n" +
        		"		cont$covmax <- apply(cont[,5:(n+4)],1,\"max\")\n" +
        		"	}\n" +
        		"	return(list(\"test\"=test,\"cont\"=cont))\n" +
        		"}\n" +
        	        "#------------------------------------------------------------------------------\n" +
        		"sdgenes <- function(listgenes, n, test, cont){\n" +
        		"	std <- c()\n" +
        		"	for(i in listgenes){\n" +
        		"		auxtest <- test[test[,4]==i,]\n" +
        		"		auxcont <- cont[cont[,4]==i,]\n" +
        		"		regions <- names(table(auxtest$start))\n" +
        		"		x <- c()\n" +
        		"		for(j in regions){\n" +
        		"			auxxt <- auxtest[auxtest[,2]==j,]\n" +
        		"			auxxc <- auxcont[auxcont[,2]==j,]\n" +
        		"			xxx <- auxxc$normmean\n" +
        		"			ratios <- auxxt$normcov/xxx\n" +
        		"			ratios <- ratios[!is.na(ratios) & ratios!=Inf]\n" +
        		"			ratios <- log(ratios)\n" +
        		"			x <- c(x, sd(ratios,na.rm=T))\n" +
        		"		}\n" +
        		"		std <- rbind(std,c(i,mean(x,na.rm=T)))\n" +
        		"	}\n" +
        		"	std <- as.data.frame(std, stringsAsFactors=F)\n" +
        		"	std[,2] <- as.numeric(as.character(gsub(\",\",\".\",std[,2])))\n" +
        		"	return(std)\n" +
        		"}\n" +
        	        "#------------------------------------------------------------------------------\n" +
        		"tabla <- function(n, listgenes, test, cont, gc, mxs, std, cov_cont){\n" +
        		"	output1 <- c()\n" +
        		"	for(k in listgenes){\n" +
        		"		m1 <- which(mxs[,1]==k)[1]\n" +
        		"		std1 <- std[std[,1]==k,2]\n" +
        		"		auxtest <- test[test[,4]==k,]\n" +
        		"		auxcont <- cont[cont[,4]==k,]\n" +
        		"		regions <- names(table(auxtest[,2]))\n" +
        		"		out <- as.data.frame(matrix(ncol=15,nrow=length(regions)))\n" +
        		"		out[,2] <- regions\n" +
        		"		out[,11] <- 0\n" +
        		"		out[,12] <- 0\n" +
        		"		for(i in regions){\n" +
        		"			cr <- auxcont[auxcont[,2]==i,]\n" +
        		"			tr <- auxtest[auxtest[,2]==i,]\n" +
        		"			mn <- 0.05*sum(cr$mean)\n" +
        		"			mx <- 0.95*sum(cr$mean)\n" +
        		"			a <- cumsum(cr$mean)\n" +
        		"			cr <- cr[which(a>mn & a<mx),]\n" +
        		"			tr <- tr[which(a>mn & a<mx),]\n" +
        		"			if(nrow(cr)!=0 && nrow(tr)!=0){\n" +
        		"				cr[which(cr==0,arr.ind=T)] <- 1.e-20\n" +
        		"				tr[which(tr==0,arr.ind=T)] <- 1.e-20\n" +
        		"				logratios <- log(tr$normcov/cr$normmean)\n" +
        		"				clusters <- picos(logratios,i)\n" +
        		"				if(length(clusters$regions)==1){\n" +
        		"					mi <- which(out[,2]==i)\n" +
        		"					rslt <- calcul(tr,cr,std1,out[mi,],gc[gc[,1]==tr[1,1],],n)\n" +
        		"					out[mi,] <- rslt\n" +
        		"				}\n" +
        		"				if(length(clusters$regions)>1){\n" +
        		"					mi <- which(out[,2]==i)\n" +
        		"					out <- out[-mi,]\n" +
        		"					mataux <- as.data.frame(matrix(ncol=ncol(out), nrow=length(clusters$regions)))\n" +
        		"					mataux[,2] <- clusters$regions\n" +
        		"					mataux[,11:12] <- 0\n" +
        		"					if(mi==1) out <- rbind(mataux,out)\n" +
        		"					if(mi>1 & mi<(nrow(out)+1)) out <- rbind(out[1:(mi-1),], mataux, out[mi:nrow(out),])\n" +
        		"					if(mi==nrow(out)+1) out <- rbind(out,mataux)\n" +
        		"					index <- which(out[,2]%in%clusters$regions)\n" +
        		"					rslt <- c()\n" +
        		"					for(j in 1:length(index)){\n" +
        		"						tr_j <- tr[which(clusters$grupos==names(table(clusters$grupos))[j]),]\n" +
        		"						cr_j <- cr[which(clusters$grupos==names(table(clusters$grupos))[j]),]\n" +
        		"						if(nrow(tr_j)!=0 && nrow(cr_j)!=0){\n" +
        		"							rslt <- rbind(rslt,calcul(tr_j,cr_j,std1,out[index[j],],gc[gc[,1]==tr_j[1,1],],n))\n" +
        		"						}\n" +
        		"					}\n" +
        		"					for(j in 1:(length(index)-1)) rslt[j,3] <- as.numeric(rslt[j+1,2])-1\n" +
        		"					indx.aux <- which(table(out[index,2])>1)\n" +
        		"					if(length(indx.aux)>0) index <- index[-which(out[index,2]%in%names(indx.aux))[2:length(which(out[index,2]%in%names(indx.aux)))]]\n" +
        		"					if(length(which(is.na(rslt[,2])))>0) rslt <- rslt[-which(is.na(rslt[,2])),]\n" +
        		"					out[index,] <- rslt\n" +
        		"				}\n" +
        		"			}\n" +
        		"			if(nrow(cr)==0 || nrow(tr)==0){\n" +
        		"				mi <- which(out[,2]==i)\n" +
        		"				out[mi,9] <- 0\n" +
        		"			}\n" +
        		"		}\n" +
        		"		out[,1] <- auxtest[1,1]\n" +
        		"		out[,4] <- k\n" +
        		"		out[,8] <- max(auxtest$cov,na.rm=T)\n" +
        		"		out[,10] <- max(auxcont$mean,na.rm=T)\n" +
        		"		out <- out[out[,9]>cov_cont,]\n" +
        		"		if(nrow(out)>0) output1 <- rbind(output1,out)\n" +
        		"	}\n" +
        		"	if(is.null(output1)) output1 <- as.data.frame(matrix(ncol=15, nrow=0))\n" +
        		"	names(output1) <- c(\"Chr\", \"Start\", \"End\", \"Gen\", \"Type\", \"pvalue\", \"Test_mean\", \"Test_gen_max\", \"Cont_mean\", \"Cont_gen_max\", \"Cont_sd\", \"CV\", \"CV_norm\", \"GC\", \"% decrease/increase\")\n" +
        		"	if(nrow(output1)!=0) output1 <- output1[which((as.numeric(output1$End)-as.numeric(output1$Start))>0),]\n" +
        		"	return(output1)\n" +
        		"}\n" +
        	        "#------------------------------------------------------------------------------\n" +
        		"picos <- function(logratios,region){\n" +
        		"	rgs <- region\n" +
        		"	if(length(table(logratios))<=2) grupos <- 0\n" +
        		"	if(length(table(logratios))>2){\n" +
        		"		m1 <- 1\n" +
        		"		m <- list()\n" +
        		"		x1 <- logratios[1]\n" +
        		"		x2 <- which(logratios*x1<0)[1]\n" +
        		"		if(is.na(x2)) x2 <- 0\n" +
        		"		pp <- 1\n" +
        		"		if(x2!=0) m[[pp]] <- 1:(x2-1)\n" +
        		"		if(x2==0) m1 <- m <- grupos <- 0\n" +
        		"		while(x2!=length(logratios) && x2!=0){\n" +
        		"			pp <- pp+1\n" +
        		"			x1 <- logratios[x2]\n" +
        		"			x <- x2 + which(logratios[(x2+1):length(logratios)]*x1<0)[1]\n" +
        		"			if(!is.na(x)){\n" +
        		"				m[[pp]] <- x2:(x-1)\n" +
        		"				x2 <- x\n" +
        		"			}\n" +
        		"			if(is.na(x)){\n" +
        		"				m[[pp]] <- x2:length(logratios)\n" +
        		"				x2 <- 0\n" +
        		"			}\n" +
        		"		}\n" +
        		"		if(m1!=0){\n" +
        		"			mdias <- c()\n" +
        		"			for(jk in 1:length(m)) mdias <- c(mdias,mean(logratios[m[[jk]]]))\n" +
        		"			difer <- abs(mdias[1:(length(mdias)-1)]-mdias[2:length(mdias)])\n" +
        		"			lll <- which(difer>=2*sd(logratios))\n" +
        		"			if(length(lll)==0){\n" +
        		"				m1 <- 0\n" +
        		"				grupos <- 0\n" +
        		"			}\n" +
        		"			if(length(lll)==1){\n" +
        		"				m1 <- list()\n" +
        		"				m1[[1]] <- unlist(m[1:lll])\n" +
        		"				m1[[2]] <- unlist(m[(lll+1):length(m)])\n" +
        		"				if(abs(mean(logratios[m1[[1]]])-mean(logratios[m1[[2]]]))<2*sd(logratios)){\n" +
        		"					m1 <- 0	\n" +
        		"					grupos <- 0\n" +
        		"				}\n" +
        		"			}\n" +
        		"			if(length(lll)>1){\n" +
        		"				m1 <- list()\n" +
        		"				m1[[1]] <- unlist(m[1:lll[1]])\n" +
        		"				for(jk in 2:length(lll)) m1[[jk]] <- unlist(m[(lll[jk-1]+1):lll[jk]])\n" +
        		"				m1[[jk+1]] <- unlist(m[(lll[jk]+1):length(m)])\n" +
        		"				mdias <- c()\n" +
        		"				for(jk in 1:length(m1)) mdias <- c(mdias,mean(logratios[m1[[jk]]]))\n" +
        		"				difer <- abs(mdias[1:(length(mdias)-1)]-mdias[2:length(mdias)])\n" +
        		"				if(length(which(difer>2*sd(logratios)))==0){\n" +
        		"					m1 <- 0\n" +
        		"					grupos <- 0\n" +
        		"				}\n" +
        		"			}\n" +
        		"		}\n" +
        		"		if(length(m1)>1){\n" +
        		"			grupos <- rep(0,length(logratios))\n" +
        		"			for(jk in 1:length(m1)){\n" +
        		"				grupos[m1[[jk]]] <- jk\n" +
        		"				rgs <- c(rgs,as.numeric(rgs[jk-1])+length(which(grupos==(jk))))\n" +
        		"			}\n" +
        		"		}\n" +
        		"	}\n" +
        		"	return(list(\"regions\"=rgs[!duplicated(rgs)], \"grupos\"=grupos))\n" +
        		"}\n" +
        	        "#------------------------------------------------------------------------------\n" +
        		"calcul <- function(tr,cr,std1,aux,gc,n){\n" +
        		"	logratios <- log(tr$normcov/cr$normmean)\n" +
			"	if(length(logratios)>5 & !is.nan(std1) & !is.na(std1) & std1>sd(logratios, na.rm=T)) std1 <- sd(logratios, na.rm=T)\n" +
        		"	media <- mean(logratios, na.rm=T)\n" +
        		"	st_t <- media/std1\n" +
        		"	aux[,6] <- 2*pnorm(abs(st_t),lower.tail=F)\n" +
        		"	aux[,3] <- cr[1,3]\n" +
        		"	aux[,7] <- mean(tr$cov)\n" +
        		"	aux[,9] <- mean(cr$mean)\n" +
        		"	aux[,11] <- mean(cr$covsd)\n" +
        		"	if(aux[,9]!=0 && n>1) aux[,12] <- aux[,11]/aux[,9]\n" +
        		"	aux[,5] <- \"deletion\"\n" +
        		"	if(max(tr$normcov,na.rm=T)>max(cr$normmean,na.rm=T)) aux[,5] <- \"duplication\"\n" +
        		"	x1 <- max(mean(tr$normcov,na.rm=T),mean(cr$normmean,na.rm=T))\n" +
        		"	x2 <- min(mean(tr$normcov,na.rm=T),mean(cr$normmean,na.rm=T))\n" +
        		"	aux[,15] <- 1 - x2/x1\n" +
        		"	r <- which(gc[,3]>=as.numeric(aux[,2]) & gc[,2]<=as.numeric(aux[,3]))\n" +
        		"	if(length(r)>0) aux[,14] <- gc[r[1],5]\n" +
        		"	aux[,13] <- 0\n" +
        		"	if(aux[,9]!=0 && n>1) aux[,13] <- sd(cr$normmean)/mean(cr$normmean)\n" +
        		"	return(aux)\n" +
        		"}\n" +
        	        "#------------------------------------------------------------------------------\n" +
        		"vcf_filter <- function(output,vcf,cov_cont){\n" +
        		"	num <- ncol(output)\n" +
        		"	output[,(num+1):(num+2)] <- NA\n" +
        		"	names(output)[(num+1):(num+2)] <- c(\"DUP_SNPs\", \"Normal_SNPs\")\n" +
        		"	el <- c()\n" +
        		"	for(i in 1:nrow(output)){\n" +
        		"		if(!is.na(output[i,1])){\n" +
        		"			vcf_aux <- vcf[vcf$CHROM==output$Chr[i] & vcf$POS>=(as.numeric(output$Start[i])-5) & vcf$POS<=(as.numeric(output$End[i])+5),]\n" +
        		"			rts <- c()\n" +
        		"			ll <- c()\n" +
        		"			if(nrow(vcf_aux)>0){\n" +
        		"				for(j in 1:nrow(vcf_aux)){\n" +
        		"					xx <- strsplit(vcf_aux$INFO[j],\";\")[[1]]\n" +
        		"					x1 <- as.numeric(gsub(\"DP=\",\"\",xx[grepl(\"DP=\",xx)]))\n" +
        		"					if(x1>cov_cont){\n" +
        		"						xx <- xx[grepl(\"DP4\",xx)]\n" +
        		"						if(length(xx)!=0){\n" +
        		"							xx <- strsplit(xx,\"\")[[1]]\n" +
        		"							a <- as.numeric(paste0(xx[(which(xx==\"=\")+1):(which(xx==\",\")[1]-1)],collapse=\"\")) + as.numeric(paste0(xx[(which(xx==\",\")[[1]]+1):(which(xx==\",\")[2]-1)],collapse=\"\"))\n" +
        		"							b <- as.numeric(paste0(xx[(which(xx==\",\")[[2]]+1):(which(xx==\",\")[3]-1)],collapse=\"\")) + as.numeric(paste0(xx[(which(xx==\",\")[[3]]+1):length(xx)],collapse=\"\"))\n" +
        		"						}\n" +
        		"						if(length(xx)==0 & T%in%grepl(\"GT:AD:\",vcf_aux$FORMAT)){\n" +
        		"							yy <- strsplit(vcf_aux$SAMPLE[j],\":\")[[1]]\n" +
        		"							yy <- yy[which(grepl(\"/\",yy))+1]\n" +
        		"							a <- as.numeric(unlist(strsplit(yy,split=\",\"))[1])\n" +
        		"							b <- as.numeric(unlist(strsplit(yy,split=\",\"))[2])\n" +
        		"						}\n" +
        		"						if(length(xx)==0 & !T%in%grepl(\"GT:AD:\",vcf_aux$FORMAT)){\n" +
        		"							a <- b <- NA\n" +
        		"						}\n" +
        		"						rts <- c(rts,a/b)\n" +
        		"					}\n" +
        		"				}\n" +
        		"				vcf_aux <- vcf_aux[ll,]\n" +
        		"				if(nrow(vcf_aux)!=0){\n" +
        		"					if(output$Type[i]==\"deletion\") el <- c(el,i)\n" +
        		"					if(output$Type[i]==\"duplication\"){\n" +
        		"						output[i,num+1] <- length(which(rts<=0.75 | rts>=1.75))\n" +
        		"						output[i,num+2] <- length(which(rts>0.75 & rts<1.75))\n" +
        		"					}\n" +
        		"				}\n" +
        		"			}\n" +
        		"		}\n" +
        		"	}\n" +
        		"	if(length(el)>0) output$pvalue[el] <- 1\n" +
        		"	output[which(output$Type==\"deletion\"),(num+1):(num+2)] <- \"-\"\n" +
        		"	output[which(output$Type==\"duplication\" & (is.na(output[,num+1]) | is.na(output[,num+2]))),(num+1):(num+2)] <- 0\n" +
        		"	return(output)\n" +
        		"}\n" +
        	        "#------------------------------------------------------------------------------\n" +
        		"concat <- function(output,cov_cont,per){\n" +
        		"	output$Start <- as.numeric(output$Start)\n" +
        		"	output$End <- as.numeric(output$End)\n" +
        		"	gns <- names(table(output$Gen))\n" +
        		"	output <- as.data.frame(cbind(output[,1:6],1,output[,7:ncol(output)]))\n" +
        		"	names(output)[7] <- \"n_exons\"\n" +
        		"	colpos <- which(grepl(\"decrease\",names(output)))\n" +
        		"	for(i in gns){\n" +
        		"		output_aux <- output[output$Gen==i,]\n" +
        		"		output_aux <- output_aux[order(output_aux[,2], decreasing=F),]\n" +
        		"		j <- 1\n" +
        		"		while(j<nrow(output_aux)){\n" +
        		"			if(!is.na(output_aux$pvalue[j]) && output_aux$pvalue[j]<0.05 && output_aux$Cont_mean[j]>=2*cov_cont && ((output_aux[j,colpos]>per[1] && output_aux$Type[j]==\"deletion\") || (output_aux[j,colpos]>per[2] && output_aux$Type[j]==\"duplication\"))){\n" +
        		"				tipo <- output_aux$Type[j]\n" +
        		"				concat <- c()\n" +
        		"				while(!is.na(output_aux$pvalue[j+1]) && output_aux$Type[j+1]==tipo && output_aux$pvalue[j+1]<0.05 && ((output_aux[j+1,colpos]>0.2 && tipo==\"deletion\") || (output_aux[j+1,colpos]>0.1 && tipo==\"duplication\")) && (j+1)<=nrow(output_aux) && output_aux[j+1,colpos]>=output_aux[j,colpos]/2){\n" +
        		"					concat <- c(concat,j,j+1)\n" +
        		"					j <- j+1\n" +
        		"				}\n" +
        		"				concat <- concat[!duplicated(concat)]\n" +
        		"				index <- which(output[,2]%in%output_aux[concat,2])\n" +
        		"				min_del <- length(which(output[index,colpos]<per[1] & output$Type[index]==\"deletion\"))\n" +
        		"				min_dup <- length(which(output[index,colpos]<per[2] & output$Type[index]==\"duplication\"))\n" +
        		"				if(length(index)>0 && (min_del+min_dup)/length(index)<0.2){\n" +
        		"					lll <- which(output$Cont_mean[index]>2*cov_cont)\n" +
        		"					if(length(lll)==length(index)) lll <- which(output$Cont_mean[index]>cov_cont)\n" +
        		"					if(length(lll)==length(index)) lll <- 1:length(index)\n" +
        		"					output$End[index[1]] <- output$End[index[length(index)]]\n" +
        		"					output$n_exons[index[1]] <- length(index)\n" +
        		"					output$Test_mean[index[1]] <- mean(output$Test_mean[index[lll]], na.rm=T)\n" +
        		"					output$Cont_mean[index[1]] <- mean(output$Cont_mean[index[lll]], na.rm=T)\n" +
        		"					output$Cont_sd[index[1]] <- mean(output$Cont_sd[index[lll]], na.rm=T)\n" +
        		"					output$CV[index[1]] <- mean(output$CV[index[lll]], na.rm=T)\n" +
        		"					output$CV_norm[index[1]] <- mean(output$CV_norm[index[lll]], na.rm=T)\n" +
        		"					d <- which(!is.na(output$GC[index]))\n" +
        		"					output$GC[index[1]] <- sum((output$End[index][d]-output$Start[index][d])*output$GC[index][d])/sum(output$End[index][d]-output$Start[index][d])\n" +
        		"					output$pvalue[index[1]] <- fisher(p=output$pvalue[index[lll]],adjust=\"generalized\", R=csconv(p = output$pvalue[index[lll]]))\n" +
        		"					output[index[1],colpos] <- max(output[index[lll],colpos],na.rm=T)\n" +
        		"					output<-output[-index[2:length(index)],]\n" +
        		"				}\n" +
        		"			}\n" +
        		"			j <- j+1\n" +
        		"		}\n" +
        		"		k <- 1\n" +
        		"		while(k <= (nrow(output[output$Gen==i,])-2)){\n" +
        		"			cond1 <- output$n_exons[output$Gen==i][k]>1 && output$n_exons[output$Gen==i][k+1]==1 && output$n_exons[output$Gen==i][k+2]>1 && output$Type[output$Gen==i][k]==output$Type[output$Gen==i][k+2] && !is.na(output$pvalue[output$Gen==i][k+1]) && output$pvalue[output$Gen==i][k+1]!=1\n" +
        		"			cond2 <- !is.na(output$pvalue[output$Gen==i][k]) && !is.na(output$pvalue[output$Gen==i][k+2]) && output$pvalue[output$Gen==i][k]<0.01 && output$pvalue[output$Gen==i][k+2]<0.01 && output$Type[output$Gen==i][k]==output$Type[output$Gen==i][k+2] && output$Type[output$Gen==i][k]==output$Type[output$Gen==i][k+1] && !is.na(output$pvalue[output$Gen==i][k+1]) && output$pvalue[output$Gen==i][k+1]!=1 && output$Cont_mean[output$Gen==i][k]>=2*cov_cont && output$Cont_mean[output$Gen==i][k+2]>=2*cov_cont && ((output[output$Gen==i,colpos][k]>per[1] && output$Type[output$Gen==i][k]==\"deletion\") || (output[output$Gen==i,colpos][k]>per[2] && output$Type[output$Gen==i][k]==\"duplication\")) && ((output[output$Gen==i,colpos][k+2]>per[1] && output$Type[output$Gen==i][k+2]==\"deletion\") || (output[output$Gen==i,colpos][k+2]>per[2] && output$Type[output$Gen==i][k+2]==\"duplication\"))\n" +
        		"			if(cond1 || cond2){\n" +
        		"				output$End[output$Gen==i][k] <- output$End[output$Gen==i][k+2]\n" +
        		"				output$n_exons[output$Gen==i][k] <- sum(output$n_exons[output$Gen==i][k:(k+2)])\n" +
        		"				output$Test_mean[output$Gen==i][k] <- mean(output$Test_mean[output$Gen==i][c(k,k+2)], na.rm=T)\n" +
        		"				output$Cont_mean[output$Gen==i][k] <- mean(output$Cont_mean[output$Gen==i][c(k,k+2)], na.rm=T)\n" +
        		"				output$Cont_sd[output$Gen==i][k] <- mean(output$Cont_sd[output$Gen==i][c(k,k+2)], na.rm=T)\n" +
        		"				output$CV[output$Gen==i][k] <- mean(output$CV[output$Gen==i][c(k,k+2)], na.rm=T)\n" +
        		"				output$CV_norm[output$Gen==i][k] <- mean(output$CV_norm[output$Gen==i][c(k,k+2)], na.rm=T)\n" +
        		"				output$GC[output$Gen==i][k] <- sum((output$End[output$Gen==i][k:(k+2)]-output$Start[output$Gen==i][k:(k+2)])*output$GC[output$Gen==i][k:(k+2)])/sum(output$End[output$Gen==i][k:(k+2)]-output$Start[output$Gen==i][k:(k+2)])\n" +
        		"				output$pvalue[output$Gen==i][k] <- fisher(p=output$pvalue[output$Gen==i][k:(k+2)], adjust=\"generalized\", R=csconv(p=output$pvalue[output$Gen==i][k:(k+2)]))\n" +
        		"				output[output$Gen==i,16][k] <- max(output[output$Gen==i,16][k:(k+2)],na.rm=T)\n" +
        		"				output <- output[-which(output[,2]%in%output[output$Gen==i,2][c(k+1,k+2)] & output$Gen==i),]\n" +
        		"			}\n" +
        		"			k <- k+1\n" +
        		"		}\n" +
        		"	}\n" +
        		"	return(output)\n" +
        		"}\n" +
        	        "#------------------------------------------------------------------------------\n" +
        		"norm_exon <- function(n, normtest, normcont, listgenes, output){\n" +
        		"	output$Start <- as.numeric(output$Start)\n" +
        		"	output$End <- as.numeric(output$End)\n" +
        		"	z <- as.data.frame(matrix(ncol = (5+n), nrow = nrow(output)))\n" +
        		"	z[, 1:4] <- output[, 1:4]\n" +
        		"	names(z)[1:5] <- c(\"Chr\", \"Start\", \"End\", \"Gen\", \"test\")\n" +
        		"	normtest$start <- normtest$start + normtest$pb - 1\n" +
        		"	normtest$end <- normtest$start + 1\n" +
        		"	normcont$start <- normtest$start\n" +
        		"	normcont$end <- normtest$end\n" +
        		"	for(k in listgenes){\n" +
        		"		lll_tst <- c()\n" +
        		"		lll_cnt <- c()\n" +
        		"		auxtest_matches <- which(normtest$gen==k)\n" +
        		"		auxcont_matches <- which(normcont$gen==k)\n" +
        		"		auxtest <- normtest[auxtest_matches,]\n" +
        		"		auxcont <- normcont[auxcont_matches,]\n" +
        		"		lll <- which(z$Gen==k)\n" +
        		"		if(length(lll)>0){\n" +
        		"			for(i in 1:nrow(z[lll,])){\n" +
        		"				matches <- which(auxtest$start >= z$Start[lll][i] & auxtest$end <= z$End[lll][i])\n" +
        		"				z[lll,5][i] <- mean(auxtest$normcov[matches])\n" +
        		"				if(length(matches)!=0) lll_tst <- c(lll_tst,matches)\n" +
        		"				matches <- which(auxcont$start>=z$Start[lll][i] & auxcont$end<=z$End[lll][i])\n" +
        		"				if(n>1) z[lll,6:(5+n)][i,] <- apply(auxcont[matches,5:(n+4)],FUN=\"mean\", MARGIN=2)\n" +
        		"				if(n==1) z[lll,6][i] <- mean(auxcont[matches,5])\n" +
        		"				if(length(matches)!=0) lll_cnt <- c(lll_cnt,matches)\n" +
        		"			}\n" +
        		"		}\n" +
        		"		if(length(lll_tst)>0) auxtest <- auxtest[-lll_tst, ]\n" +
        		"		if(length(lll_cnt)>0) auxcont <- auxcont[-lll_cnt, ]\n" +
        		"		if(length(auxtest_matches) > 0) normtest <- normtest[-auxtest_matches, ]\n" +
        		"		if(length(auxcont_matches) > 0) normcont <- normcont[-auxcont_matches, ]\n" +
        		"	}\n" +
        		"	z <- z[which(duplicated(z[,1:3]) == F),]\n" +
        		"	return(z)\n" +
        		"}\n" +
        	        "#------------------------------------------------------------------------------\n" +
        		"regression <- function(n,data){\n" +
        		"	aux <- data[,1:6]\n" +
        		"	aux[,7] <- 0\n" +
        		"	aux[,8] <- 1\n" +
        		"	names(aux)[6:8] <- c(\"cont\",\"residual\",\"penalty\")\n" +
        		"	if(n>1) aux[,6] <- apply(data[,6:(n+5)],1,\"mean\")\n" +
        		"	nans <- as.numeric(names(table(which(is.na(aux),arr.ind=T)[,1])))\n" +
        		"	if(length(nans)>0) aux <- aux[-nans,]\n" +
        		"	mod <- lm(aux$test~aux$cont)\n" +
        		"	if(!is.na(summary(mod)$adj.r.squared) & summary(mod)$adj.r.squared>0.8){\n" +
        		"		aux[,7] <- mod$residuals\n" +
        		"		y1 <- mod$fitted.values-sd(mod$fitted.values)\n" +
        		"		y2 <- mod$fitted.values+sd(mod$fitted.values)\n" +
        		"		y3 <- mod$fitted.values-qnorm(0.95)*sd(mod$fitted.values)\n" +
        		"		y4 <- mod$fitted.values+qnorm(0.95)*sd(mod$fitted.values)\n" +
        		"		m <- which(aux$test<=y1 & aux$test>y3)\n" +
        		"		if(length(m)>0) aux[m,8][order(aux[m,7],decreasing=T)] <- 1/(abs(aux$residual[m])*5e2)\n" +
        		"		m <- which(aux$test>=y2 & aux$test<y4)\n" +
        		"		if(length(m)>0) aux[m,8][order(aux[m,7],decreasing=F)] <- 1/(abs(aux$residual[m])*5e2)\n" +
        		"		m <- which(aux$test<=y3)\n" +
        		"		if(length(m)>0) aux[m,8][order(aux[m,7],decreasing=T)] <- 1/(abs(aux$residual[m])*1e3)\n" +
        		"		m <- which(aux$test>=y4)\n" +
        		"		if(length(m)>0) aux[m,8][order(aux[m,7],decreasing=F)] <- 1/(abs(aux$residual[m])*1e3)\n" +
        		"		y11 <- mod$fitted.values-0.15*sd(mod$fitted.values)\n" +
        		"		y21 <- mod$fitted.values+0.15*sd(mod$fitted.values)\n" +
        		"		m <- which(aux$test>y11 & aux$test<=mod$fitted.values)\n" +
        		"		if(length(m)>0) aux[m,8] <- (1-abs(aux$residual[m]))*1.e3\n" +
        		"		m <- which(aux$test<y21 & aux$test>=mod$fitted.values)\n" +
        		"		if(length(m)>0) aux[m,8] <- (1-abs(aux$residual[m]))*1.e3\n" +
        		"		y12 <- mod$fitted.values-0.5*sd(mod$fitted.values)\n" +
        		"		y22 <- mod$fitted.values+0.5*sd(mod$fitted.values)\n" +
        		"		m <- which(aux$test>y12 & aux$test<=y11)\n" +
        		"		if(length(m)>0) aux[m,8] <- (1-abs(aux$residual[m]))*5.e2\n" +
        		"		m <- which(aux$test<y22 & aux$test>=y21)\n" +
        		"		if(length(m)>0) aux[m,8] <- (1-abs(aux$residual[m]))*5.e2\n" +
        		"	}\n" +
        		"	return(aux)\n" +
        		"}\n" +
        	        "#------------------------------------------------------------------------------\n" +
        		"filt <- function(output,per,normcont,n,residuos,sex,restr,cov_cont){\n" +
        		"	nn <- nrow(output)\n" +
        		"	output <- output[output$Test_gen_max>0,]\n" +
        		"	if(per[1]>=0.3 && per[2]>=0.26) output <- output[c(which(output[,16]>=0.3 & output$Type==\"deletion\"), which(output[,16]>=0.26 & output$Type==\"duplication\")),]\n" +
        		"	output[,16] <- 100*output[,16]\n" +
        		"	output[,c(8:16)] <- round(output[,c(8:16)],2)\n" +
        		"	output <- output[order(output[,6]),]\n" +
        		"	hlp <- 1\n" +
        		"	if(nrow(output[output[,6]<0.05,])>5 & nrow(output[output[,6]<0.05,])<=10){\n" +
        		"		hlp <- output[,6]*nn/c(1:nrow(output))\n" +
        		"		h <- \"benj\"\n" +
        		"	}\n" +
        		"	if(nrow(output[output[,6]<0.05,])>=10){\n" +
        		"		hlp <- output[,6]*nn\n" +
        		"		h <- \"bonf\"\n" +
        		"	}\n" +
        		"	if(length(which(hlp<0.05))!=0){\n" +
        		"		output[,6] <- hlp\n" +
        		"		names(output)[6] <- paste(names(output)[6],h,sep=\"_\")\n" +
        		"	}\n" +
        		"	for(i in 1:nrow(output)){\n" +
        		"		nr <- which(residuos[,1]==output$Chr[i] & residuos[,2]==output$Start[i] & residuos[,3]==output$End[i])\n" +
        		"		if(length(nr)>0) output[i,6] <- output[i,6]*residuos$penalty[nr]\n" +
        		"	}\n" +
        		"	output <- output[output[,6]<0.05,]\n" +
        		"	output <- output[c(which(output[,16]>=100*per[1] & output$Type==\"deletion\"), which(output[,16]>=100*per[2] & output$Type==\"duplication\")),]\n" +
        		"	aux <- rep(0,nrow(output))\n" +
        		"	if(n>1){\n" +
        		"		for(i in 1:length(aux)){\n" +
        		"			auxcont <- normcont[normcont$chr==output$Chr[i] & normcont$start<=output$End[i] & normcont$end>=output$Start[i],]\n" +
        		"			x <- c()\n" +
        		"			for(j in 1:(n-1)){\n" +
        		"				a1 <- mean(auxcont[,4+j],na.rm=T)\n" +
        		"				for(k in (j+1):n){\n" +
        		"					a2 <- mean(auxcont[,4+k],na.rm=T)\n" +
        		"					x <- c(x, 1-min(c(a1,a2))/max(c(a1,a2)))\n" +
        		"				}\n" +
        		"			}\n" +
        		"			aux[i] <- mean(x,na.rm=T)\n" +
        		"		}\n" +
        		"	}\n" +
        		"	aux <- round(100*aux,2)\n" +
        		"	output <- output[which(output[,16]>aux),]\n" +
        		"	if(restr){\n" +
        		"		if(per[1]<=0.4 & length(which(output$Type==\"deletion\"))!=0) output <- output[c(which(output[,16]>=40 & output$Type==\"deletion\"), which(output$Type==\"duplication\")),]\n" +
        		"		if(per[2]<=0.35 & length(which(output$Type==\"duplication\"))!=0) output <- output[c(which(output$Type==\"deletion\"), which(output[,16]>=35 & output$Type==\"duplication\")),]\n" +
        		"		if(cov_cont<=100) output <- output[output$Cont_mean>=100,]\n" +
        		"		if(nrow(output)>10 && nrow(output)<30) output <- output[output[,6]<median(output[,6], na.rm=T),]\n" +
        		"		if(nrow(output)>=30) output <- output[output[,6]<quantile(output[,6], probs=0.25, na.rm=T),]\n" +
        		"	}\n" +
        		"	output <- output[order(output[,6]),]\n" +
        		"	output <- output[order(output$Type),]\n" +
        		"	if(nrow(output[output$Type==\"deletion\",])>0) output <- rbind(output[output$Type==\"deletion\",], rep(NA,ncol(output)), output[output$Type==\"duplication\",])\n" +
        		"	return(output)\n" +
        		"}\n" +
        	        "#------------------------------------------------------------------------------\n" +
        		"downsampling <- function(name_test, name_cont, n, bed, output1, cov_cont, gc, per, sex, restr){\n" +
        		"	outputdwn <- list()\n" +
        		"	cr <- as.data.frame(matrix(ncol=2,nrow=length(table(bed[,1]))))\n" +
        		"	cr[,1] <- names(table(bed[,1]))\n" +
        		"	cr[,2] <- gsub(\"chr\",\"\",cr[,1])\n" +
        		"	if(length(which(cr[,2]==\"X\"))!=0) cr[which(cr[,2]==\"X\"),2] <- 23\n" +
        		"	if(length(which(cr[,2]==\"Y\"))!=0) cr[which(cr[,2]==\"Y\"),2] <- 24\n" +
        		"	cr[,2] <- as.numeric(cr[,2])\n" +
        		"	cname <- paste0(name_cont, \"_dwn.bam\")\n" +
        		"	for(k in 1:n){\n" +
        		"		dat <- coverage(readGAlignments(cname[k]))\n" +
        		"		aux.cont <- c()\n" +
        		"		for(crom in cr[,1]){\n" +
        		"			dat1 <- as.numeric(dat[[cr[cr[,1]==crom,2]]])\n" +
        		"			bed1 <- bed[bed[,1]==crom,]\n" +
        		"			for(j in 1:nrow(bed1)){\n" +
        		"				auxx <- as.data.frame(matrix(ncol=6,nrow=(bed1[j,3]-bed1[j,2])))\n" +
        		"				auxx[,1] <- bed1[j,1]\n" +
        		"				auxx[,4] <- bed1[j,4]\n" +
        		"				auxx[,2] <- bed1[j,2]\n" +
        		"				auxx[,3] <- bed1[j,3]\n" +
        		"				auxx[,5] <- 1:nrow(auxx)\n" +
        		"				auxx[,6] <- dat1[auxx[,2]+(1:nrow(auxx))]\n" +
        		"				aux.cont <- rbind(aux.cont,auxx)\n" +
        		"			}\n" +
        		"		}\n" +
        		"		if(k==1){\n" +
        		"			cont <- as.data.frame(matrix(ncol=4+n,nrow=nrow(aux.cont)))\n" +
        		"			cont[,1:4] <- aux.cont[,1:4]\n" +
        		"		}\n" +
        		"		cont[,4+k] <- aux.cont[,6]\n" +
        		"	}\n" +
        		"	if(n>1){\n" +
        		"		cont[,(5+n)] <- apply(cont[,5:(4+n)],1,mean)\n" +
        		"		names(cont) <- c(\"chr\", \"start\", \"end\", \"gen\", name_cont, \"mean\")\n" +
        		"	}\n" +
        		"	if(n==1){\n" +
        		"		cont$mean <- cont[,5]\n" +
        		"		names(cont)[1:5] <- c(\"chr\", \"start\", \"end\", \"gen\", name_cont)\n" +
        		"	}\n" +
        		"	cont$gen <- toupper(cont$gen)\n" +
        		"	cont$chr <- toupper(cont$chr)\n" +
        		"	cont$chr <- gsub(\"CHR\",\"chr\",cont$chr)\n" +
        		"	if(length(grep(\"chr\",cont$chr))==0) cont$chr <- paste0(\"chr\",cont$chr)\n" +
        		"	for(i in c(1,2,3,5,8)){\n" +
        		"		dat <- coverage(readGAlignments(paste(c(name_test, \"_dwn\", i, \".bam\"), collapse=\"\")))\n" +
        		"		test <- c()\n" +
        		"		for(crom in cr[,1]){\n" +
        		"			dat1 <- as.numeric(dat[[cr[cr[,1]==crom,2]]])\n" +
        		"			bed1 <- bed[bed[,1]==crom,]\n" +
        		"			for(j in 1:nrow(bed1)){\n" +
        		"				auxx <- as.data.frame(matrix(ncol=6,nrow=(bed1[j,3]-bed1[j,2])))\n" +
        		"				auxx[,1] <- bed1[j,1]\n" +
        		"				auxx[,4] <- bed1[j,4]\n" +
        		"				auxx[,2] <- bed1[j,2]\n" +
        		"				auxx[,3] <- bed1[j,3]\n" +
        		"				auxx[,5] <- 1:nrow(auxx)\n" +
        		"				auxx[,6] <- dat1[auxx[,2]+(1:nrow(auxx))]\n" +
        		"				test <- rbind(test,auxx)\n" +
        		"			}\n" +
        		"		}\n" +
        		"		names(test) <- c(\"chr\",\"start\",\"end\",\"gen\",\"pb\",\"cov\")\n" +
        		"		test$gen <- toupper(test$gen)\n" +
        		"		test$chr <- toupper(test$chr)\n" +
        		"		test$chr <- gsub(\"CHR\",\"chr\",test$chr)\n" +
        		"		if(length(grep(\"chr\",test$chr))==0) test$chr <- paste0(\"chr\",test$chr)\n" +
        		"		listgenes <- names(table(test[,4]))\n" +
        		"		fun_mxs <- rgns(listgenes, n, test, cont, c(), cov_cont, F)\n" +
        		"		mxs <- as.data.frame(fun_mxs$mxs, stringsAsFactors=F, row.names=F)\n" +
        		"		normalized <- nrml(listgenes, n, test, cont, mxs)\n" +
        		"		std <- sdgenes(listgenes, n, normalized$test, normalized$cont)\n" +
        		"		output_dwn <- tabla(n, listgenes, normalized$test, normalized$cont, gc, mxs, std, cov_cont)\n" +
        		"		output_dwn <- concat(output_dwn, cov_cont, per)\n" +
        		"		nexons <- norm_exon(n, normalized$test, normalized$cont, listgenes, output_dwn)\n" +
        		"		residuos <- regression(n, nexons)\n" +
        		"		output_dwn <- filt(output_dwn, per, normalized$cont, n, residuos, sex, restr, cov_cont)\n" +
        		"		names(output_dwn)[6] <- \"pvalue\"\n" +
        		"		outputdwn[[i]] <- output_dwn\n" +
        		"	}\n" +
        		"	outputdwn <- do.call(\"rbind\", outputdwn)\n" +
        		"	outputdwn <- outputdwn[!is.na(outputdwn[,1]),]\n" +
        		"	eldwn <- c()\n" +
        		"	for(i in 1:nrow(output1)){\n" +
        		"		if(!is.na(output1$Chr[i]) & length(which(outputdwn$Chr==output1$Chr[i] & outputdwn$Start<=output1$End[i] & outputdwn$End>=output1$Start[i]))<4){\n" +
        		"			l <- which(bed[,4]==output1$Gen[i] & (bed[,2]>output1$End[i] | bed[,3]<output1$Start[i]))\n" +
        		"			if(length(l)==0) eldwn <- c(eldwn,i)\n" +
        		"			if(length(l)!=0){\n" +
        		"				aux.dwn <- outputdwn[outputdwn$Gen==output1$Gen[i] & outputdwn$Start<=max(bed[l,3]) & outputdwn$End>=min(bed[l,2]),]\n" +
        		"				if(nrow(aux.dwn)==0) eldwn <- c(eldwn,i)\n" +
        		"				if(nrow(aux.dwn)!=0 && (table(aux.dwn$Start)<4 | names(table(aux.dwn$Type))[1]==output1$Type[i])) eldwn <- c(eldwn,i)\n" +
        		"			}\n" +
        		"		}\n" +
        		"	}\n" +
        		"	if(length(eldwn)!=0) output1 <- output1[-eldwn,]\n" +
        		"	if(is.na(output1[1,1])) output1 <- output1[-1,]\n" +
        		"	return(output1)\n" +
        		"}\n" +
        	        "#------------------------------------------------------------------------------\n" +
        	        "freq <- function(output, userDB, passwordDB){\n" +
        	        "	dbPrueba<-dbConnect(MySQL(), user=userDB, password=passwordDB, host=\"localhost\", dbname=\"CNV\")\n" +
        	        "	gr <- output[,1:5]\n" +
        	        "	gr$Start <- as.numeric(sub(\",\", \".\", gr$Start, fixed = TRUE))\n" +
        	        "	gr$End <- as.numeric(sub(\",\", \".\", gr$End, fixed = TRUE))\n" +
        	        "	m1 <- m2 <- m3 <- c()\n" +
        	        "	for(i in 1:nrow(gr)){\n" +
        	        "		if(is.na(gr[i,1])){\n" +
        	        "			m1 <- c(m1,NA)\n" +
        	        "			m2 <- c(m2,NA)\n" +
        	        "			m3 <- c(m3,NA)\n" +
        	        "		}\n" +
        	        "		if(!is.na(gr[i,1])){\n" +
        	        "			consultaString = paste0(\"select * from CNV_global where Type=\\\"\",gr$Type[i],\"\\\" and Start<=\", gr$End[i], \" and End>=\", gr$Start[i], \";\")\n" +
        	        "			consultaSql <- dbSendQuery(dbPrueba, consultaString)\n" +
        	        "			aux <-fetch(consultaSql, n=-1)\n" +
        	        "			if(nrow(aux)>0){\n" +
        	        "				aux$Chr <- toupper(aux$Chr)\n" +
        	        "				aux$Chr <- paste0(\"chr\",gsub(\"CHR\",\"\",aux$Chr))\n" +
        	        "				aux <- aux[aux$Chr==paste0(\"chr\",gsub(\"CHR\",\"\",toupper(gr$Chr[i]))),]\n" +
        	        "			}\n" +
        	        "			smpls <- names(table(aux$Sample))\n" +
        	        "			m1 <- c(m1,length(smpls))\n" +
        	        "			if(length(smpls)>0 && length(smpls)<10){\n" +
        	        "				m2 <- c(m2, paste(smpls,collapse=\";\"))\n" +
        	        "				per_mx <- c()\n" +
        	        "				for(j in smpls){\n" +
        	        "					per_mx <- c(per_mx, max(aux$Percentage[aux$Sample==j],na.rm=T))\n" +
        	        "				}\n" +
        	        "				m3 <- c(m3, paste(per_mx,collapse=\";\"))\n" +
        	        "			}\n" +
        	        "			if(length(smpls)==0){\n" +
        	        "				m2 <- c(m2,\"\")\n" +
        	        "				m3 <- c(m3,\"\")\n" +
        	        "			}\n" +
        	        "			if(length(smpls)>=10){\n" +
        	        "				m2 <- c(m2,\"-\")\n" +
        	        "				m3 <- c(m3,\"-\")\n" +
        	        "			}\n" +
        	        "		}\n" +
        	        "	}\n" +
        	        "	dbDisconnect(dbPrueba)\n" +
        	        "	output <- cbind(output,m1,m2,m3)\n" +
        	        "	names(output)[(ncol(output)-2):ncol(output)] <- c(\"Freq\", \"Samples\", \"Percentage_Samples\")\n" +
        	        "	return(output)\n" +
        	        "}\n" +
        	        "#------------------------------------------------------------------------------\n" +
        		"DCVar_exon <- function(output, n, data){\n" +
        		"	output[,ncol(output)+1] <- NA\n" +
        		"	names(output)[ncol(output)] <- \"DCVar\"\n" +
        		"	output$Start <- as.numeric(output$Start)\n" +
        		"	output$End <- as.numeric(output$End)\n" +
        		"	for(i in 1:nrow(output)){\n" +
        		"		if(!is.na(output[i,1])){\n" +
        		"			aux <- data[data$Chr==output$Chr[i] & data$Start<=output$End[i] & data$End>=output$Start[i],]\n" +
        		"			if(nrow(aux)>1) aux <- aux[which(aux[,2]==output[i,2]),]\n" +
        		"			if(nrow(aux)>1) aux <- aux[which(aux[,3]==output[i,3]),]\n" +
        		"			z1 <- mean(as.numeric(aux[,6:(n+5)]))\n" +
        		"			if(z1>0) vl <- sd(as.numeric(aux[,6:(n+5)]),na.rm=T)/z1\n" +
        		"			if(z1==0) vl <- 0\n" +
        		"			output$DCVar[i] <- vl\n" +
        		"		}\n" +
        		"	}\n" +
        		"	return(output)\n" +
        		"}\n" +
        	        "#------------------------------------------------------------------------------\n" +
        		"annot <- function(output, userDB, passwordDB){\n" +
        		"	dbPrueba<-dbConnect(MySQL(), user=userDB, password=passwordDB, host=\"localhost\", dbname=\"CNV\")\n" +
        		"	for(i in 1:nrow(output)){\n" +
        		"		if(!is.na(output[i,1])){\n" +
        		"			consultaString = paste0(\"select * from exones_list where start<=\", output$End[i], \" and end>=\", output$Start[i], \";\")\n" +
        		"			consultaSql <- dbSendQuery(dbPrueba, consultaString)\n" +
        		"			aux <- fetch(consultaSql, n=-1)\n" +
        		"			if(nrow(aux)>0){\n" +
        		"				aux$chr <- toupper(aux$chr)\n" +
        		"				aux$chr <- paste0(\"chr\",gsub(\"CHR\",\"\",aux$chr))\n" +
        		"				aux <- aux[aux$chr==paste0(\"chr\",gsub(\"CHR\",\"\",toupper(output$Chr[i]))),]\n" +
        		"				ifelse(nrow(aux)>0, output$Region_info[i] <- paste(paste(aux$exon_number,paste0(aux$ID,\")\"),sep=\"(\"),collapse=\";\"), output$Region_info[i] <- NA)\n" +
        		"			}\n" +
        		"			if(nrow(aux)==0) output$Region_info[i] <- NA\n" +
        		"		}\n" +
        		"		if(is.na(output[i,1])) output$Region_info[i] <- NA\n" +
        		"	}\n" +
        		"	dbDisconnect(dbPrueba)\n" +
        		"	return(output)\n" +
        		"}\n" +
        	        "#------------------------------------------------------------------------------\n" +
        		"information <- function(n, test, cont, name_test, name_cont, sex, nexons){\n" +
        		"	corr <- cor(cbind(test$cov,cont[,5:(4+n)]))\n" +
        		"	info_aux <- as.data.frame(matrix(ncol=6, nrow=length(name_cont)+1))\n" +
        		"	info_aux[,1] <- c(paste(\"Test\",name_test,sep=\"-\"),paste0(\"Cont\",paste(1:n,name_cont,sep=\"-\")))\n" +
        		"	info_aux[,2] <- sex\n" +
        		"	if(n==1){\n" +
        		"		info_aux[,3] <- c(mean(test$cov),mean(cont$mean))\n" +
        		"		info_aux[,4] <- c(sd(test$cov),sd(cont$mean))\n" +
        		"	}\n" +
        		"	if(n>1){\n" +
        		"		info_aux[,3] <- c(mean(test$cov),apply(cont[,5:(4+n)],MARGIN=2,FUN=\"mean\"))\n" +
        		"		info_aux[,4] <- c(sd(test$cov),apply(cont[,5:(4+n)],MARGIN=2,FUN=\"sd\"))\n" +
        		"	}\n" +
        		"	info_aux[,3] <- round(info_aux[,3],2)\n" +
        		"	info_aux[,4] <- round(info_aux[,4],2)\n" +
        		"	info_aux[,5] <- round(info_aux[,4]/info_aux[,3],3)\n" +
        		"	info_aux[,6] <- c(\"-\",round(corr[1,2:(1+n)],3))\n" +
        		"	names(info_aux) <- c(\"Sample\", \"Sex\", \"Cov_mean\", \"Cov_SD\", \"CV\",\"Correlation\")\n" +
        		"	info_aux[(nrow(info_aux)+1):(nrow(info_aux)+2),] <- NA\n" +
        		"	info_aux[nrow(info_aux),1] <- \"DCVar_global(Mean)\"\n" +
        		"	dcvar_aux <- apply(nexons[,5:(5+n)], FUN=\"sd\", MARGIN=1)/apply(nexons[,5:(5+n)],FUN=\"mean\",MARGIN=1)\n" +
        		"	info_aux[nrow(info_aux),2] <- round(mean(dcvar_aux, na.rm=T),4)\n" +
        		"	info_aux[nrow(info_aux)+1,1] <- \"DCVar_global(SD)\"\n" +
        		"	info_aux[nrow(info_aux),2] <- round(sd(dcvar_aux, na.rm=T),4)\n" +
        		"	return(info_aux)\n" +
        		"}\n" +
        	        "#------------------------------------------------------------------------------\n" +
        		"wrt <- function(per, output, output_genes, name_test, name_cont, genes, info){\n" +
        		"	if(genes==T){\n" +
        		"		write.table(output_genes, file = \"output_genes.txt\", append = FALSE, quote = TRUE, sep = \"\\t\",eol = \"\\n\", na = \"\", dec = \",\", row.names = F,col.names = TRUE, qmethod = c(\"escape\", \"double\"),fileEncoding = \"\")\n" +
        		"	}\n" +
        		"	write.table(info, file = \"output_information.txt\", append = FALSE, quote = TRUE, sep = \"\\t\",eol = \"\\n\", na = \"\", dec = \",\", row.names = F,col.names = TRUE, qmethod = c(\"escape\", \"double\"),fileEncoding = \"\")\n" +
        		"	write.table(output, file = \"output_regions.txt\", append = FALSE, quote = TRUE, sep = \"\\t\",eol = \"\\n\", na = \"\", dec = \",\", row.names = F,col.names = TRUE, qmethod = c(\"escape\", \"double\"),fileEncoding = \"\")\n" +
        		"}\n" +
        	        "#------------------------------------------------------------------------------\n" +
        		"write.excel <- function(output, name_test, name_cont, datos_dcvar, per, output_genes, genes, info){\n" +
        		"	for(i in 3:5) info[,i] <- as.character(info[,i])\n" +
        		"	m <- nrow(output)\n" +
        		"	if(m!=0){\n" +
        		"		if(is.na(output[m,1])){\n" +
        		"			m <- m-1\n" +
        		"			output <- output[1:m,]\n" +
        		"		}\n" +
        		"	}\n" +
        		"	dcvar_c <- which(names(output)==\"DCVar\")\n" +
        		"	freq_c <- which(names(output)==\"Freq\")\n" +
        		"	file <- \"output.xlsx\"\n" +
        		"	wb <- createWorkbook(file)\n" +
        		"	addWorksheet(wb, paste(\"regions\",paste0(paste(100*per,collapse=\"-\"),\"%\"),sep=\"_\"))\n" +
        		"	output1 <- output\n" +
        		"	for(i in 1:ncol(output1)) output1[,i] <- as.character(output1[,i])\n" +
        		"	writeData(wb, 1, output1)\n" +
        		"	cs1 <- createStyle(fgFill = \"#FF3030\")\n" +
        		"	cs2 <- createStyle(fgFill = \"#00CDCD\")\n" +
        		"	cs3 <- createStyle(fgFill = \"#FF7F00\")\n" +
        		"	if(m!=0){\n" +
        		"		highlightred <- which(output[,15]<=0.35 || output[,15]>=0.6)\n" +
        		"		addStyle(wb, 1, style = cs1, rows = highlightred+1, cols = 15, gridExpand = TRUE)\n" +
        		"		highlightred <- which(output[,8]<50 && output[,10]<50)\n" +
        		"		addStyle(wb, 1, style = cs1, rows = highlightred+1, cols = c(8,10), gridExpand = TRUE)\n" +
        		"		highlightred <- which(output[,10]<50)\n" +
        		"		addStyle(wb, 1, style = cs1, rows = highlightred+1, cols = 10, gridExpand = TRUE)\n" +
        		"		highlightorange <- which(output[,10]>=50 && output[,10]<100)\n" +
        		"		addStyle(wb, 1, style = cs3, rows = highlightorange+1, cols = 10, gridExpand = TRUE)\n" +
        		"		highlightblue <- which(output[,10]>=100)\n" +
        		"		addStyle(wb, 1, style = cs2, rows = highlightblue+1, cols = 10, gridExpand = TRUE)\n" +
        		"		highlightblue <- which(output[,16]>=40 && !is.na(output[,16]))\n" +
        		"		addStyle(wb, 1, style = cs2, rows = highlightblue+1, cols = 16, gridExpand = TRUE)\n" +
        		"		highlightblue <- which(output[,14]<=0.01 && !is.na(output[,14]))\n" +
        		"		addStyle(wb, 1, style = cs2, rows = highlightblue+1, cols = 14, gridExpand = TRUE)\n" +
        		"	}\n" +
        		"	if(length(freq_c)!=0){\n" +
        		"		highlightblue <- which(output[,freq_c]<=5 && !is.na(output[,freq_c]))\n" +
        		"		addStyle(wb, 1, style = cs2, rows = highlightblue+1, cols = freq_c, gridExpand = TRUE)\n" +
        		"		highlightred <- which(output[,freq_c]>=10 && !is.na(output[,freq_c]))\n" +
        		"		addStyle(wb, 1, style = cs1, rows = highlightred+1, cols = freq_c, gridExpand = TRUE)\n" +
        		"	}\n" +
        		"	if(length(dcvar_c)!=0){\n" +
        		"		datos_dcvar <- as.numeric(as.character(gsub(\",\",\".\",info[(which(is.na(info[,2]))[1]+1):nrow(info),2])))\n" +
        		"		highlightred <- which(output[,dcvar_c]>sum(datos_dcvar) && !is.na(output[,dcvar_c]))\n" +
        		"		addStyle(wb, 1, style = cs1, rows = highlightred+1, cols = dcvar_c, gridExpand = TRUE)\n" +
        		"		highlightblue <- which(output[,dcvar_c]<datos_dcvar[1] && !is.na(output[,dcvar_c]))\n" +
        		"		addStyle(wb, 1, style = cs2, rows = highlightblue+1, cols = dcvar_c, gridExpand = TRUE)\n" +
        		"	}\n" +
        		"	freezePane(wb, 1, firstActiveRow=2)\n" +
        		"	if(genes==T){\n" +
        		"		output_genes1 <- output_genes\n" +
        		"		for(i in 5:ncol(output_genes1)) output_genes1[,i] <- as.character(output_genes1[,i])\n" +
        		"		addWorksheet(wb, \"genes\")\n" +
        		"		writeData(wb, 2, output_genes1)\n" +
        		"		highlightred <- which(output_genes[,11]<=0.35 || output_genes[,11]>=0.6)\n" +
        		"		addStyle(wb, 1, style = cs1, rows = highlightred+1, cols = 11, gridExpand = TRUE)\n" +
        		"		highlightred <- which(output_genes[,7]<50 && output_genes[,8]<50)\n" +
        		"		addStyle(wb, 1, style = cs1, rows = highlightred+1, cols = 7:8, gridExpand = TRUE)\n" +
        		"		highlightred <- which(output_genes[,8]<50)\n" +
        		"		addStyle(wb, 1, style = cs1, rows = highlightred+1, cols = 8, gridExpand = TRUE)\n" +
        		"		highlightorange <- which(output_genes[,8]>=50 && output_genes[,8]<100)\n" +
        		"		addStyle(wb, 1, style = cs3, rows = highlightorange+1, cols = 8, gridExpand = TRUE)\n" +
        		"		highlightblue <- which(output_genes[,8]>=100)\n" +
        		"		addStyle(wb, 1, style = cs2, rows = highlightblue+1, cols = 8, gridExpand = TRUE)\n" +
        		"		highlightblue <- which(output_genes[,12]>=40)\n" +
        		"		addStyle(wb, 1, style = cs2, rows = highlightblue+1, cols = 12, gridExpand = TRUE)\n" +
        		"		freezePane(wb, 2, firstActiveRow=2)\n" +
        		"		addWorksheet(wb, \"information\")\n" +
        		"		writeData(wb, 3, info)\n" +
        		"		freezePane(wb, 3, firstActiveRow=2)\n" +
        		"	}\n" +
        		"	if(genes==F){\n" +
        		"		addWorksheet(wb, \"information\")\n" +
        		"		writeData(wb, 2, info)\n" +
        		"		freezePane(wb, 2, firstActiveRow=2)\n" +
        		"	}\n" +
        		"	saveWorkbook(wb, file, overwrite=T)\n" +
        		"}\n" +
        	        "#------------------------------------------------------------------------------\n" +
        		"graphs <- function(output, name_test, name_cont, n, test, cont, mxs){\n" +
        		"	listgns <- names(table(output[,4]))\n" +
        		"	for(k in listgns){\n" +
        		"		cairo_pdf(paste(c(\"Gene \", k, \".pdf\"), collapse=\"\"))\n" +
        		"		auxtest <- test[test[,4]==k,]\n" +
        		"		auxtest <- auxtest[order(auxtest[,2]),]\n" +
        		"		auxcont <- cont[cont[,4]==k,]\n" +
        		"		auxcont <- auxcont[order(auxcont[,2]),]\n" +
        		"		regions <- names(table(auxtest[,2]))\n" +
        		"		rgs <- table(auxtest$start)\n" +
        		"		rgs <- cumsum(rgs)\n" +
        		"		rgs1 <- rgs\n" +
        		"		for(i in 2:length(rgs1)) rgs1[i] <- rgs[i-1]+(rgs[i]-rgs[i-1])/2\n" +
        		"		rgs1[1] <- rgs1[1]/2\n" +
        		"		aux_gen <- auxtest[,6]\n" +
        		"		if(n>1) aux_gen <- cbind(aux_gen,auxcont[,5:(4+n)])\n" +
        		"		if(n==1) aux_gen <- cbind(aux_gen,auxcont[,5])\n" +
        		"		aux_gen <- cbind(1:nrow(aux_gen),aux_gen)\n" +
        		"		mx_gen <- max(aux_gen[,2:(2+n)])+50\n" +
        		"		par(mfrow=c(2,1))\n" +
        		"		if(n>1) color <- rainbow(n+1)\n" +
        		"		if(n==1) color <- c(\"red\", \"blue\")\n" +
        		"		plot(aux_gen[,2]~aux_gen[,1], ylab=\"raw coverage\", ylim=c(-50,mx_gen), type=\"l\", main=paste0(k, \" (raw coverage)\"), xlab=\"\", xlim=c(0,nrow(aux_gen)), col=\"white\", xaxt=\"n\", cex.axis=0.7)\n" +
        		"		for(i in 2:(2+n)){\n" +
        		"			points(aux_gen[,i]~aux_gen[,1], col=color[i-1], type=\"l\")\n" +
        		"		}\n" +
        		"		xg <- apply(aux_gen[,2:ncol(aux_gen)], MARGIN=1, FUN=\"max\")\n" +
        		"		ifelse(max(xg[1:round(length(xg)/4)])>max(xg[(3*round(length(xg)/4)):length(xg)]), psx <- \"topright\", psx <- \"topleft\")\n" +
        		"		legend(x=psx, legend=c(\"Test\", paste0(\"Control\",1:n)), col=color, lty=1, cex=0.7, bty=\"n\")\n" +
        		"		abline(v=c(0,rgs),col=\"darkgrey\",lty=2)\n" +
        		"		axis(1,at=rgs1[which(!names(rgs1)%in%output$Start)],labels=names(rgs1[which(!names(rgs1)%in%output$Start)]),las=2, cex.axis=0.6)\n" +
        		"		axis(1,at=rgs1[which(names(rgs1)%in%output$Start)], labels=names(rgs1[which(names(rgs1)%in%output$Start)]), las=2, cex.axis=0.6, col.axis=\"darkred\")\n" +
        		"		if(n>1) xxx <- apply(aux_gen[,3:(n+2)],1,mean)\n" +
        		"		if(n==1) xxx <- aux_gen[,3]\n" +
        		"		aux_gen[,2] <- aux_gen[,2]/as.numeric(mxs[which(mxs[,1]==k),2])\n" +
        		"		if(n>1){\n" +
        		"			for(i in 3:(n+2)) aux_gen[,i] <- aux_gen[,i]/as.numeric(mxs[which(mxs[,1]==k),i])\n" +
        		"			xxx <- apply(aux_gen[,3:(n+2)],1,mean)\n" +
        		"		}\n" +
        		"		if(n==1) xxx <- aux_gen[,3]/as.numeric(mxs[which(mxs[,1]==k),3])\n" +
        		"		plot(aux_gen[,2]~aux_gen[,1], ylab=\"normalized coverage\", ylim=c(0,1.1*max(c(aux_gen[,2],xxx))), type=\"l\", main=paste0(k, \" (normalized coverage)\"), xlab=\"\", xlim=c(0,nrow(aux_gen)), col=\"red\", xaxt=\"n\", cex.axis=0.7)\n" +
        		"		points(xxx~aux_gen[,1], col=\"blue\", type=\"l\")\n" +
        		"		legend(x=psx, legend=c(\"Test\",\"Control(s)\"), col=c(\"red\", \"blue\"), lty=1, cex=0.7, bty=\"n\")\n" +
        		"		abline(v=c(0,rgs),col=\"darkgrey\",lty=2)\n" +
        		"		axis(1,at=rgs1[which(!names(rgs1)%in%output$Start)],labels=names(rgs1[which(!names(rgs1)%in%output$Start)]),las=2,cex.axis=0.6)\n" +
        		"		axis(1,at=rgs1[which(names(rgs1)%in%output$Start)], labels=names(rgs1[which(names(rgs1)%in%output$Start)]), las=2, cex.axis=0.6, col.axis=\"darkred\")\n" +
        		"		dev.off()\n" +
        		"	}\n" +
        		"}\n" +
        	        "#------------------------------------------------------------------------------\n" +
        	        "result <- ppal(fvcf, per=c(per_del,per_dup), per_gene=c(per_gene_del, per_gene_dup), cov_cont, name_test, name_cont, name_bed, name_fasta, excel, genes, polymorphic, poly_regions, do_plot, restr, fixed, down)\n" +
        	        "return(dim(result$output))\n" +
        	        "rm(list=ls(all=T))}";
    }
    
    public int execute (String per_dup, String per_del, String per_gene_dup, String per_gene_del, String cov_cont, 
            String fvcf, String excel, String genes, String restr, String fixed, String down, String do_plot, String userDB, String passwordDB, String wd, String name_test, String name_bed, String name_fasta, String[] name_cont, String polymorphic, String poly_regions) {
        // Start R session.
        
        Rengine re = Rengine.getMainEngine();
        if(re == null) {
            re = new Rengine(new String[] {"--vanilla"}, false, this);
        }
        
        // Check if the session is working.
        if (!re.waitForR()) {
            return -1;
        }
		
        // assign arguments to R 
        re.assign("per_dup", per_dup);
        re.assign("per_del", per_del);
        re.assign("per_gene_dup", per_gene_dup);
        re.assign("per_gene_del", per_gene_del);
        re.assign("cov_cont", cov_cont);
        re.assign("fvcf", fvcf);
        re.assign("excel", excel);
        re.assign("genes", genes);
        re.assign("restr", restr);
        re.assign("fixed", fixed);
        re.assign("down", down);
        re.assign("do_plot", do_plot);
        re.assign("userDB", userDB);
        re.assign("passwordDB", passwordDB);
        re.assign("wd", wd);
        re.assign("name_test", name_test);
        re.assign("name_bed", name_bed);
        re.assign("name_fasta", name_fasta);
        re.assign("name_cont", name_cont);
        re.assign("polymorphic", polymorphic);
        re.assign("poly_regions", poly_regions);
        re.eval(rCode);
        //print=true;
        REXP result = re.eval("scriptR(per_dup, per_del, per_gene_dup, per_gene_del, cov_cont, fvcf, excel, genes, restr, fixed, down, do_plot, userDB, passwordDB, wd, name_test, name_bed, name_fasta, name_cont, polymorphic, poly_regions)");
        re.jriFlushConsole();
        print=false;
        //re.end();
        
        if (result!=null) {
            return 0; //Exit code OK
        } else {
            return -1; //Exit code ERROR
        }
    }
    
    @Override
    public void rWriteConsole(Rengine rngn, String string, int i) {
        if (print) {
            System.out.println(string);
        }
    }

    @Override
    public void rBusy(Rengine rngn, int i) {
		//Don't generate exception --> Silence output messages when we are installing new packages
        //throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public String rReadConsole(Rengine rngn, String string, int i) {
        System.out.println("GRIDD rReadConsole = \n" + string);
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void rShowMessage(Rengine rngn, String string) {
        System.out.println(string);    }

    @Override
    public String rChooseFile(Rengine rngn, int i) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void rFlushConsole(Rengine rngn) {
        //throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void rSaveHistory(Rengine rngn, String string) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void rLoadHistory(Rengine rngn, String string) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
}
