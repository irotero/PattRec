package org.gradiant.threads;

import org.gradiant.UI.RCNV;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.File;
import java.io.FilenameFilter;
import java.io.InputStreamReader;
import javax.swing.JOptionPane;
import javax.swing.JTextArea;
import org.gradiant.utils.Commons;
import org.gradiant.JRI.RCNVScript;
import org.gradiant.JRI.RDBScript;

/**
 *
 * @author mpacheco
 */
public class R_Thread extends NotifyingThread {

    private RCNV r_data;
    private String destFolder;
    private ProcessBuilder builder;
    private Process process;
    final String message = "Do you want to save in the DDBB the results of this simulation?";

    public R_Thread(RCNV R_data) {
        r_data = R_data;
	String s1=r_data.getFile_Test().getName();
 	String replaceString=s1.replace(".bam","");
        destFolder = Commons.baseDir + "output/" + replaceString;
    }

    public void doRun() {

        try {
            makeDirectories();
            
            // check bedtools version
            boolean btNewVersion = Commons.checkBedtoolsNewVersion();
            boolean stNewVersion = Commons.checkSamtoolsNewVersion();
            doCoverageControls(btNewVersion);
            doCoverageTest(btNewVersion);
            doGCcontent();

            if (r_data.isFvcf()) {
                if(r_data.getFilename_vcf() == null)
                    doSamTools(stNewVersion);
                else
                    copyVCF();
            }
            
            if (r_data.isDo_down()) {
                doDownsampling(btNewVersion);
            }
            
            if (r_data.isPolymorphic()) {
                if (Commons.checkFileExists(r_data.getFilename_poly())){
                    doRCNV();
                } else {
                    r_data.getOutput().append("Polymorphic file doesn't exist\n"); 
                    r_data.getOutput().append("The program has failed\n");
                }
            } else {
                doRCNV();
            }
            
            doRemoveTmp();
            
        } catch (IOException | InterruptedException ex) {
            ex.printStackTrace();
        }
    }

    public void makeDirectories() throws IOException, InterruptedException {        
        Runtime.getRuntime().exec("rm -rf " + destFolder).waitFor();
        Runtime.getRuntime().exec("mkdir -p " + destFolder).waitFor();
        Runtime.getRuntime().exec("mkdir -p " + destFolder + "/tmp").waitFor();
        Runtime.getRuntime().exec("mkdir -p " + destFolder + "/results").waitFor();
    }

    public void doCoverageControls(boolean btNewVersion) throws IOException, InterruptedException {
        r_data.getOutput().append("Executing BedTools for control samples...\n");
        File[] controls = r_data.getFiles_control();

        for (File control : controls) {
            if (btNewVersion == true){
                builder = new ProcessBuilder("bedtools", "coverage", "-a", r_data.getFilename_BED(), "-b", control.getAbsolutePath(), "-d");
            } else {
                builder = new ProcessBuilder("coverageBed", "-abam", control.getAbsolutePath(), "-b", r_data.getFilename_BED(), "-d");
            }
            builder.redirectOutput(new File(destFolder + "/tmp/" + control.getName() + ".txt"));
            builder.start().waitFor();
        }
    }

    public void doCoverageTest(boolean btNewVersion) throws IOException, InterruptedException {
        r_data.getOutput().append("Executing BedTools for test sample...\n");
        if (btNewVersion == true){
            builder = new ProcessBuilder("bedtools", "coverage", "-a", r_data.getFilename_BED(), "-b", r_data.getFile_Test().getAbsolutePath(), "-d");
        } else {
            builder = new ProcessBuilder("coverageBed", "-abam", r_data.getFile_Test().getAbsolutePath(), "-b", r_data.getFilename_BED(), "-d");

        }
        builder.redirectOutput(new File(destFolder + "/tmp/" + r_data.getFile_Test().getName() + ".txt"));
        builder.start().waitFor();
    }

    public void doGCcontent() throws IOException, InterruptedException {
        r_data.getOutput().append("Obtaining GC content...\n");
        builder = new ProcessBuilder("bedtools", "nuc", "-fi", r_data.getFilename_fasta(), "-bed", r_data.getFilename_BED());
        builder.redirectOutput(new File(destFolder + "/tmp/GC_TestSample.txt"));
        builder.start().waitFor();
    }

    public void doSamTools(boolean stNewVersion) throws IOException, InterruptedException {
        r_data.getOutput().append("Executing SamTools for test sample...\n");
        String[] cmd;
        if (stNewVersion == true) {
            cmd = new String[] {"/bin/sh", "-c", "samtools mpileup -ugf \"" + r_data.getFilename_fasta() + "\" \"" + r_data.getFile_Test().getAbsolutePath()
                + "\" | bcftools call -vc - > " + destFolder + "/tmp/" + r_data.getFile_Test().getName() + ".vcf"};

        } else {
            cmd = new String[] {"/bin/sh", "-c", "samtools mpileup -ugf \"" + r_data.getFilename_fasta() + "\" \"" + r_data.getFile_Test().getAbsolutePath()
                + "\" | bcftools view -vcg - > " + destFolder + "/tmp/" + r_data.getFile_Test().getName() + ".vcf"};

        }
        process = Runtime.getRuntime().exec(cmd);
        process.waitFor();
    }
    
    public void doDownsampling(boolean btNewVersion) throws IOException, InterruptedException {
        r_data.getOutput().append("Downsampling control samples...\n");
        File[] controls = r_data.getFiles_control();
        for (File control : controls) {
		String s1=control.getName();
 		String replaceString=s1.replace(".bam","");
                Runtime.getRuntime().exec("java " + "-jar " + Commons.libDir + "/picard.jar " + "DownsampleSam I=" + control.getAbsolutePath() + " O=" + destFolder + "/tmp/" + replaceString + "_dwn.bam " + "P=0.8 " + "R=1").waitFor();
		if (btNewVersion == true){
                	builder = new ProcessBuilder("bedtools", "coverage", "-a", r_data.getFilename_BED(), "-b", destFolder + "/tmp/" + replaceString + "_dwn.bam", "-d");
            	} else {
                	builder = new ProcessBuilder("coverageBed", "-abam", destFolder + "/tmp/" + replaceString + "_dwn.bam", "-b", r_data.getFilename_BED(), "-d");
            	}
            	builder.redirectOutput(new File(destFolder + "/tmp/" + replaceString + "_dwn.txt"));
            	builder.start().waitFor();
	}
        r_data.getOutput().append("Downsampling test sample...\n");
	String s1=r_data.getFile_Test().getName();
 	String replaceString=s1.replace(".bam","");
        Runtime.getRuntime().exec("java " + "-jar " + Commons.libDir + "/picard.jar " + "DownsampleSam I=" + r_data.getFile_Test().getAbsolutePath() + " O=" + destFolder + "/tmp/" + replaceString + "_dwn.bam " + "P=0.8 " + "R=1").waitFor();
       	if (btNewVersion == true){
		builder = new ProcessBuilder("bedtools", "coverage", "-a", r_data.getFilename_BED(), "-b", destFolder + "/tmp/" + replaceString + "_dwn.bam", "-d");
       	} else {
       		builder = new ProcessBuilder("coverageBed", "-abam", destFolder + "/tmp/" + replaceString + "_dwn.bam", "-b", r_data.getFilename_BED(), "-d");
	}
       	builder.redirectOutput(new File(destFolder + "/tmp/" + replaceString + "_dwn1.txt"));
	builder.start().waitFor();
        Runtime.getRuntime().exec("java " + "-jar " + Commons.libDir + "/picard.jar " + "DownsampleSam I=" + r_data.getFile_Test().getAbsolutePath() + " O=" + destFolder + "/tmp/" + replaceString + "_dwn.bam " + "P=0.8 " + "R=2").waitFor();
       	if (btNewVersion == true){
		builder = new ProcessBuilder("bedtools", "coverage", "-a", r_data.getFilename_BED(), "-b", destFolder + "/tmp/" + replaceString + "_dwn.bam", "-d");
       	} else {
       		builder = new ProcessBuilder("coverageBed", "-abam", destFolder + "/tmp/" + replaceString + "_dwn.bam", "-b", r_data.getFilename_BED(), "-d");
	}
       	builder.redirectOutput(new File(destFolder + "/tmp/" + replaceString + "_dwn2.txt"));
	builder.start().waitFor();
        Runtime.getRuntime().exec("java " + "-jar " + Commons.libDir + "/picard.jar " + "DownsampleSam I=" + r_data.getFile_Test().getAbsolutePath() + " O=" + destFolder + "/tmp/" + replaceString + "_dwn.bam " + "P=0.8 " + "R=3").waitFor();
       	if (btNewVersion == true){
		builder = new ProcessBuilder("bedtools", "coverage", "-a", r_data.getFilename_BED(), "-b", destFolder + "/tmp/" + replaceString + "_dwn.bam", "-d");
       	} else {
       		builder = new ProcessBuilder("coverageBed", "-abam", destFolder + "/tmp/" + replaceString + "_dwn.bam", "-b", r_data.getFilename_BED(), "-d");
	}
       	builder.redirectOutput(new File(destFolder + "/tmp/" + replaceString + "_dwn3.txt"));
	builder.start().waitFor();
        Runtime.getRuntime().exec("java " + "-jar " + Commons.libDir + "/picard.jar " + "DownsampleSam I=" + r_data.getFile_Test().getAbsolutePath() + " O=" + destFolder + "/tmp/" + replaceString + "_dwn.bam " + "P=0.8 " + "R=5").waitFor();
       	if (btNewVersion == true){
		builder = new ProcessBuilder("bedtools", "coverage", "-a", r_data.getFilename_BED(), "-b", destFolder + "/tmp/" + replaceString + "_dwn.bam", "-d");
       	} else {
       		builder = new ProcessBuilder("coverageBed", "-abam", destFolder + "/tmp/" + replaceString + "_dwn.bam", "-b", r_data.getFilename_BED(), "-d");
	}
       	builder.redirectOutput(new File(destFolder + "/tmp/" + replaceString + "_dwn5.txt"));
	builder.start().waitFor();
        Runtime.getRuntime().exec("java " + "-jar " + Commons.libDir + "/picard.jar " + "DownsampleSam I=" + r_data.getFile_Test().getAbsolutePath() + " O=" + destFolder + "/tmp/" + replaceString + "_dwn.bam " + "P=0.8 " + "R=8").waitFor();
       	if (btNewVersion == true){
		builder = new ProcessBuilder("bedtools", "coverage", "-a", r_data.getFilename_BED(), "-b", destFolder + "/tmp/" + replaceString + "_dwn.bam", "-d");
       	} else {
       		builder = new ProcessBuilder("coverageBed", "-abam", destFolder + "/tmp/" + replaceString + "_dwn.bam", "-b", r_data.getFilename_BED(), "-d");
	}
       	builder.redirectOutput(new File(destFolder + "/tmp/" + replaceString + "_dwn8.txt"));
	builder.start().waitFor();
    }

    public void doRCNV() throws IOException, InterruptedException {
        r_data.getOutput().append("Executing R script...\n");
        
        
        String[] controls = new String[r_data.getFiles_control().length];
        for (int i=0; i<controls.length ;i++) {
            controls[i] = r_data.getFiles_control()[i].getName(); //args[14-?] = control sample files
        }
        
        RCNVScript rcnv = new RCNVScript();
        String poly_regions = Commons.baseDir + "/input/polymorphic.bed"; 
        int exitValue = rcnv.execute(Double.toString(r_data.getPer_dup()), Double.toString(r_data.getPer_del()), Double.toString(r_data.getPer_gen_dup()), 
                Double.toString(r_data.getPer_gen_del()), Integer.toString(r_data.getCov_cont()), r_data.isFvcf() ? "SI" : "NO", 
                r_data.isDo_excel() ? "TRUE" : "FALSE", r_data.isGenes() ? "TRUE" : "FALSE", r_data.isDo_restr() ? "TRUE" : "FALSE", 
                r_data.isDo_fixed() ? "TRUE" : "FALSE", r_data.isDo_down() ? "TRUE" : "FALSE", r_data.isDo_plot() ? "TRUE" : "FALSE", 
                r_data.getCred().getUser(), r_data.getCred().getPassword(), destFolder + "/tmp", r_data.getFile_Test().getName(), r_data.getFilename_fasta(), 
                controls, r_data.isPolymorphic() ? "SI" : "NO", r_data.getFilename_poly());
        
        doMoveFiles();
        
        if(exitValue == 0) {
            r_data.getOutput().append("Done!\n");
            if(JOptionPane.showConfirmDialog(r_data, message, "Save results", JOptionPane.YES_NO_OPTION) == 0) {
                if(saveResults() == 0)
                    r_data.getOutput().append("Results successfully saved in the database\n");
                else
                    r_data.getOutput().append("There was an error storing the results in the database\n");
            }

        } else {
            r_data.getOutput().append("The program has failed\n");
        }
    }
    
    public void doRemoveTmp() throws IOException, InterruptedException {
        Runtime.getRuntime().exec(new String[]{"/bin/sh", "-c", "rm -r " + destFolder + "/tmp/"}).waitFor();
    }
    
    private static void printProcOutput(JTextArea area, Process proc) throws IOException {
        //BufferedReader reader = new BufferedReader(new InputStreamReader(proc.getInputStream()));
        BufferedReader reader = new BufferedReader(new InputStreamReader(proc.getErrorStream()));
        String line = null;
        while ( (line = reader.readLine()) != null) {
           area.append(line);
           area.append(System.getProperty("line.separator"));
        }
    }

    public void copyVCF() throws IOException, InterruptedException {
        r_data.getOutput().append("Using vcf file already existing...\n");
        
        String vcf = r_data.getFilename_vcf();
        Runtime.getRuntime().exec("cp " + vcf + " " + destFolder + "/tmp").waitFor();
        
        vcf = vcf.split("/")[vcf.split("/").length - 1]; // Get the final name only
        Runtime.getRuntime().exec("mv " + destFolder + "/tmp/" + vcf + " " + destFolder + "/tmp/" + r_data.getFile_Test().getName() + ".vcf").waitFor();
    }

    private int saveResults() throws IOException, InterruptedException {
        String sampleName = JOptionPane.showInputDialog(r_data, "Give an ID to save this sample:", null, JOptionPane.QUESTION_MESSAGE);
        if(sampleName != null && !sampleName.equals("")) {
            r_data.getOutput().append("Saving results in the database...\n");
            
            String results;
            if(r_data.isDo_excel())
                results = destFolder + "/results/" + getFilesByExtension(destFolder + "/results", ".xlsx")[0];
            else
                results = destFolder + "/results/" + getFilesByExtension(destFolder + "/results", "regions.txt")[0];
            
            RDBScript rdbs = new RDBScript();
            int exitValue = rdbs.execute(r_data.getCred().getUser(), r_data.getCred().getPassword(), 
                sampleName, results);

            return exitValue;
        }

        return -1;
    }
    
    private static String[] getFilesByExtension(String dir, String ext) {
        File folder = new File(dir);

        if(folder.isDirectory() == false)
            return null;

        return folder.list(new FilenameFilter() {
            public boolean accept(File dir, String name) {
                    return (name.endsWith(ext));
            }
        });
    }
    
    public void doMoveFiles() throws IOException, InterruptedException {
        Runtime.getRuntime().exec(new String[]{"/bin/sh", "-c", "mv " + destFolder + "/tmp/*.xlsx " + destFolder + "/results/"}).waitFor();
        Runtime.getRuntime().exec(new String[]{"/bin/sh", "-c", "mv " + destFolder + "/tmp/*.pdf " + destFolder + "/results/"}).waitFor();
        Runtime.getRuntime().exec(new String[]{"/bin/sh", "-c", "mv " + destFolder + "/tmp/*gen.txt " + destFolder + "/results/"}).waitFor();
        Runtime.getRuntime().exec(new String[]{"/bin/sh", "-c", "mv " + destFolder + "/tmp/*information.txt " + destFolder + "/results/"}).waitFor();
        Runtime.getRuntime().exec(new String[]{"/bin/sh", "-c", "mv " + destFolder + "/tmp/*regions.txt " + destFolder + "/results/"}).waitFor();
        
    }
}
