/*
 * Copyright (c) 2016 GRADIANT. All rights reserved.
 * This code cannot be used, copied, modified and/or distributed without the express permission of the authors.
 * This algorithm is protected under a Confidentiality Agreement between the UDTEMC (Fundaci√≥n Ram√≥n Dom√≠nguez) and GRADIANT.
 */
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

import com.google.common.io.Files;

import org.gradiant.JRI.RCNVScript;
import org.gradiant.JRI.RDBScript;

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
    //Modificado el cÛdigo para incluir File.separator en lugar de depender del sistema operativo
        destFolder = Commons.baseDir + "output" + File.separator + replaceString;
    }

    public void doRun() {

        try {
            makeDirectories();
            if (r_data.isFvcf()) {
                if(r_data.getFilename_vcf() == null) {
                    r_data.getOutput().append("The vcf file doesn't exist: ignoring vcf filtering...\n"); 
                } else {
                    copyVCF();
            }}

            if (r_data.isDo_down()) {
                    doDownsampling();
            }
            
            if (r_data.isPolymorphic()) {
                if (Commons.checkFileExists(r_data.getFilename_poly())){
                    doRCNV();
                } else {
                    r_data.getOutput().append("Polymorphic file doesn't exist: ignoring polymorphic filtering...\n");
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
        r_data.getOutput().append("Starting analysis...\n");
        File fileDestFolder = new File (destFolder);
        deleteDirectory(fileDestFolder);
        fileDestFolder.mkdirs();        
        File fileDestFolder1 = new File (destFolder + File.separator + "tmp");
        fileDestFolder1.mkdir();
        File fileDestFolder2 = new File (destFolder + File.separator + "results");
        fileDestFolder2.mkdir();
        //Runtime.getRuntime().exec("rm -rf " + destFolder).waitFor();
        //Runtime.getRuntime().exec("mkdir -p " + destFolder).waitFor();
        //Runtime.getRuntime().exec("mkdir -p " + destFolder + "/tmp").waitFor();
       // Runtime.getRuntime().exec("mkdir -p " + destFolder + "/results").waitFor();
        //System.out.println("Rscript " + Commons.resDir + "res" + File.separator + "input" + File.separator + "install_packages.R");
        
        //Inicialmente se descargaba aqui el paquete ShortRead pero tras un error que no supimos identificar se moviÛ al Script de R de la clase RCNVScript
        //el error era debido a que Rscript tiene un bug a la hora de interpretar directorios que tengan nombres con espacios en blanco
        //Runtime.getRuntime().exec("Rscript \"" + Commons.resDir + "res" + File.separator + "input" + File.separator + "install_packages.R\"").waitFor();
    }

    public void doDownsampling() throws IOException, InterruptedException {
        r_data.getOutput().append("Downsampling controls...\n");
        File[] controls = r_data.getFiles_control();
        for (File control : controls) {
        	String s1=control.getName();
        	String replaceString=s1.replace(".bam","");
        	System.out.println("java " + "-jar \"" + Commons.resDir + "lib" + File.separator + "picard.jar\" " + "DownsampleSam I=\"" + control.getAbsolutePath() + "\" O=\"" + destFolder + File.separator + "tmp" + File.separator + replaceString + "_dwn.bam\" " + "P=0.8 " + "R=1");
        	Runtime.getRuntime().exec("java " + "-jar \"" + Commons.resDir + "lib" + File.separator + "picard.jar\" " + "DownsampleSam I=\"" + control.getAbsolutePath() + "\" O=\"" + destFolder + File.separator + "tmp" + File.separator + replaceString + "_dwn.bam\" " + "P=0.8 " + "R=1").waitFor();
        }
        r_data.getOutput().append("Downsampling test...\n");
        String s1=r_data.getFile_Test().getName();
        String replaceString=s1.replace(".bam","");
        Runtime.getRuntime().exec("java " + "-jar \"" + Commons.resDir + "lib" + File.separator + "picard.jar\" " + "DownsampleSam I=\"" + r_data.getFile_Test().getAbsolutePath() + "\" O=\"" + destFolder + File.separator + "tmp" + File.separator + replaceString + "_dwn1.bam\" " + "P=0.8 " + "R=1").waitFor();
        Runtime.getRuntime().exec("java " + "-jar \"" + Commons.resDir + "lib" + File.separator + "picard.jar\" " + "DownsampleSam I=\"" + r_data.getFile_Test().getAbsolutePath() + "\" O=\"" + destFolder + File.separator + "tmp" + File.separator + replaceString + "_dwn2.bam\" " + "P=0.8 " + "R=2").waitFor();
        Runtime.getRuntime().exec("java " + "-jar \"" + Commons.resDir + "lib" + File.separator + "picard.jar\" " + "DownsampleSam I=\"" + r_data.getFile_Test().getAbsolutePath() + "\" O=\"" + destFolder + File.separator + "tmp" + File.separator + replaceString + "_dwn3.bam\" " + "P=0.8 " + "R=3").waitFor();
        Runtime.getRuntime().exec("java " + "-jar \"" + Commons.resDir + "lib" + File.separator + "picard.jar\" " + "DownsampleSam I=\"" + r_data.getFile_Test().getAbsolutePath() + "\" O=\"" + destFolder + File.separator + "tmp" + File.separator + replaceString + "_dwn5.bam\" " + "P=0.8 " + "R=5").waitFor();
        Runtime.getRuntime().exec("java " + "-jar \"" + Commons.resDir + "lib" + File.separator + "picard.jar\" " + "DownsampleSam I=\"" + r_data.getFile_Test().getAbsolutePath() + "\" O=\"" + destFolder + File.separator + "tmp" + File.separator + replaceString + "_dwn8.bam\" " + "P=0.8 " + "R=8").waitFor();
    }



    public void doRCNV() throws IOException, InterruptedException {
        r_data.getOutput().append("Executing R script...\n");
        
        
       String[] controls = new String[r_data.getFiles_control().length];
       for (int i=0; i<controls.length ;i++) {
           controls[i] = r_data.getFiles_control()[i].getAbsolutePath(); //args[14-?] = control sample files
        }
        
        RCNVScript rcnv = new RCNVScript();
        String poly_regions = Commons.baseDir + "input" + File.separator + "polymorphic.bed";
        int exitValue = rcnv.execute(Double.toString(r_data.getPer_dup()), Double.toString(r_data.getPer_del()), Double.toString(r_data.getPer_gen_dup()), 
                Double.toString(r_data.getPer_gen_del()), Integer.toString(r_data.getCov_cont()), r_data.isFvcf() ? "SI" : "NO", 
                r_data.isDo_excel() ? "TRUE" : "FALSE", r_data.isGenes() ? "TRUE" : "FALSE", r_data.isDo_restr() ? "TRUE" : "FALSE", r_data.isDo_fixed() ? "TRUE" : "FALSE", r_data.isDo_down() ? "TRUE" : "FALSE", r_data.isDo_plot() ? "TRUE" : "FALSE", r_data.getCred().getUser(), r_data.getCred().getPassword(), 
                destFolder + File.separator + "tmp", r_data.getFile_Test().getAbsolutePath(), r_data.getFilename_BED(), r_data.getFilename_fasta(), controls, r_data.isPolymorphic() ? "SI" : "NO", r_data.getFilename_poly());
        
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
    	//Cambiamos el cÛdigo para no ser dependientes del SO (Unix), para ello llamamos a la funcion deleteDirectory
        //Runtime.getRuntime().exec(new String[]{File.separator + "bin" + File.separator + "sh", "-c", "rm -r " + destFolder + File.separator + "tmp" + File.separator}).waitFor();
    	File deleteFile = new File(destFolder + File.separator + "tmp" + File.separator);
    	deleteDirectory(deleteFile);
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
        r_data.getOutput().append("Using vcf file...\n");
        
        String vcf = r_data.getFilename_vcf();         
        File originFolder= new File(vcf);
        
        String finalNameVcf = originFolder.getName();     
        File destinyFolder = new File(destFolder + File.separator + "tmp"+ File.separator + finalNameVcf);
        
    	if (originFolder.exists()) {
    		Files.copy(originFolder, destinyFolder);
    	}
    	
    	//Modificamos el cÛdigo de manera que no se esa dependiende del SO
        //Runtime.getRuntime().exec("cp " + vcf + " " + destFolder + "/tmp").waitFor();
        //Runtime.getRuntime().exec("mv " + destFolder + "/tmp/" + vcf + " " + destFolder + "/tmp/" + r_data.getFile_Test().getName() + ".vcf").waitFor();
    }

    private int saveResults() throws IOException, InterruptedException {
        String sampleName = JOptionPane.showInputDialog(r_data, "Give an ID to save this sample:", null, JOptionPane.QUESTION_MESSAGE);
        if(sampleName != null && !sampleName.equals("")) {
            r_data.getOutput().append("Saving results in the database...\n");
            
            String results;
            if(r_data.isDo_excel()) {
                results = destFolder + File.separator + "results" + File.separator + getFilesByExtension(destFolder + File.separator + "results", ".xlsx")[0];
//                System.out.println(results);
//                System.out.println("");
//                System.out.println("");
//                System.out.println(r_data.getCred().getUser());
//                System.out.println(r_data.getCred().getPassword());
//                System.out.println(sampleName);
            }
            else
                results = destFolder + File.separator + "results" + File.separator + getFilesByExtension(destFolder + File.separator + "results", "regions.txt")[0];
            
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
 
//FunciÛn modificada para no depender del SO (Unix)
   public void doMoveFiles() throws IOException, InterruptedException {
	  
	   File destDir = new File (destFolder + File.separator + "results");	   
	   File originFolder = new File (destFolder + File.separator + "tmp");
       File[] filesInDir = originFolder.listFiles();
       for(File file:filesInDir) {
           String name = file.getName();
           if (name.contains(".xlsx")||name.contains(".pdf")||name.contains("gen.txt")||name.contains("information.txt")||name.contains("regions.txt")) {
               //System.out.println(name);
        	   String newPath = destDir + File.separator + name;
        	   file.renameTo(new File(newPath));
           }
       }
        
        /*Runtime.getRuntime().exec(new String[]{"/bin/sh", "-c", "mv " + destFolder + "/tmp/*.xlsx " + destFolder + "/results/"}).waitFor();
        Runtime.getRuntime().exec(new String[]{"/bin/sh", "-c", "mv " + destFolder + "/tmp/*.pdf " + destFolder + "/results/"}).waitFor();
        Runtime.getRuntime().exec(new String[]{"/bin/sh", "-c", "mv " + destFolder + "/tmp/*gen.txt " + destFolder + "/results/"}).waitFor();
        Runtime.getRuntime().exec(new String[]{"/bin/sh", "-c", "mv " + destFolder + "/tmp/*information.txt " + destFolder + "/results/"}).waitFor();
        Runtime.getRuntime().exec(new String[]{"/bin/sh", "-c", "mv " + destFolder + "/tmp/*regions.txt " + destFolder + "/results/"}).waitFor();  */   
    }
    
 //Nuevo mÈtodo para desempeÒar la misma funciÛn que el programa rm de Unix y no depender por tanto del SO
    private boolean deleteDirectory(File directoryToBeDeleted) {
        File[] allContents = directoryToBeDeleted.listFiles();
        if (allContents != null) {
            for (File file : allContents) {
                deleteDirectory(file);
            }
        }
        return directoryToBeDeleted.delete();
    }
    
}
