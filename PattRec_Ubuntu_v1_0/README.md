PattRec
=============================

DESCRIPTION
------------
PattRec is a bioinformatics tool designed to detect rare copy number variants (CNVs) in targeted Next Generation Sequencing (tg-NGS) data. It is presented as a Java-based GUI, with its CNV detection algorithm implemented in R.
This tool was designed for use with target gene panels, sequenced in Illumina platforms.



SYSTEM REQUIREMENTS
-------------------
	- **Ubuntu 14.04 LTS 64bit**, **Ubuntu 16.04 LTS 64bit** or **Ubuntu 18.04 LTS 64bit** operating system
	- Java 8:
		https://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html
	- MySQL or MariaDB: 
		https://dev.mysql.com/downloads/
		https://mariadb.org/download/
	- R Project:
		https://www.r-project.org/
	- Bedtools:
		https://bedtools.readthedocs.io/en/latest/
	- SAMtools:
		https://sourceforge.net/projects/samtools/
	- 8GB RAM or greater is recommended.


RUN
------------
Some configurations have to be performed before running the program: 

1. If you haven't done it yet, install rJava:

    Terminal:            $ R CMD INSTALL lib/rJava_0.9-11.tar.gz
	R:					 $ install.packages("rJava")

2. Set configurations (change the path if needed):

		$ export R_HOME=/usr/lib/R
		$ export LD_LIBRARY_PATH=/usr/lib/R/site-library/rJava/jri/


3. Run the jar file:
	
		$ java -jar pattrec.jar
	



USER MANUAL
------------
1. Configure database: click "Database configuration" and enter user and password for MySQL or MariaDB database.

2. Select input files:
	- BAM test: BAM file of the test sample (patient).
	- BAM control: one or more BAM files to compare with the BAM test.
		
		*NOTE: It is highly recommended that all the BAM files were sequenced on the same run.*
	- BED file: file containing the target regions (sequenced regions) in BED format.
	- FASTA file: genome reference file (it must be the same used for the alignment of the BAM files).
		
		*NOTE: GRCh37-hg19 can be downloaded from: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/*
		*GRCh38-hg38 can be downloaded from: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/*

3. Click "RUN".

4. Configure parameters for the algorithm (if needed).

5. Click "Proceed".

6. You will be asked to save results in database. If you agree, that results will be used in the next execution of the algorithm.

Results will be placed in the user home, in the folder "PattRec".



