PattRec
=============================

DESCRIPTION
------------
PattRec is a bioinformatics tool designed to detect rare copy number variants (CNVs) in targeted Next Generation Sequencing (tg-NGS) data. It is presented as a Java-based GUI, with its CNV detection algorithm implemented in R.
This tool was designed for use with target gene panels, sequenced in Illumina platforms.



SYSTEM REQUIREMENTS
-------------------
	- **Windows 10 64bit** operating system
	- 8GB RAM or greater is recommended



RUN
------------
PattRec depends on R v3.4 (with rJava package), Java 8, MySQL or MariaDB, Visual Studio, and Perl 5.12. The executable allows the user to install all the dependencies, or just those needed.

1. If you already have all or some of the dependencies installed on your computer, please uncheck the corresponding box, and the installer will ask you to select the path to R (v3.4.3) and Java 8. Please make sure that you have installed the R package rJava; if you haven't, just open a R session and type: 'install.packages("rJava")'.

2. If you choose to install MySQL via the installer, please make sure to:
	- Select 'Typical' installation
	- Launch the 'Wizard' in order to configurate the database.
	- Select 'Detailed Configuration'
	- Add MySQL to the firewall
	- Create a password for 'root'


USER MANUAL
------------
1. Configure database: click "Database configuration" and enter user and password for MySQL or MariaDB database.

2. Select input files:
	- BAM test: BAM file of the test sample (patient).
	- BAM control: one or more BAM files to compare with the BAM test.

		*NOTE: It is highly recommended that all the BAM files were sequenced on the same run.*
	- BED file: file containing the target regions (sequenced regions) tab separated format.

		*IMPORTANT: the header of the BED file must be removed, so the first line in the file corresponds to a region*
		*Example: chr1	10000	20000	Gene_name* 
	- FASTA file: genome reference file (it must be the same used for the alignment of the BAM files).
		*NOTE: GRCh37-hg19 can be downloaded from: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/*
		*GRCh38-hg38 can be downloaded from: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/*
		

3. Click "RUN".

4. Configure parameters for the algorithm (if needed).

5. Click "Proceed".

6. You will be asked to save results in database. If you agree, that results will be used in the next execution of the algorithm.

Results will be placed in the user home, in the folder "PattRec".



