/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.gradiant.utils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.io.File;

/**
 *
 * @author lgonzalez
 */
public class Commons {

    public final static String baseDir = System.getProperty("user.home") + "/PattRec/";
    public final static String defaultPolyFile = baseDir + "input/polymorphic.bed";
    
    public final static String defaultPicardJar = baseDir + "input/picard.jar";
    
    // path to the lib directory.
    // development: set this as current directory and place lib/ACNV/ and lib/RCNV/ folder inside "dist/"
    // release: set this in accordance to the structure of debian package. for example "/usr/lib/GRIDD/"
    public final static String libDir = System.getProperty("user.dir") + "/lib/"; // DEV
    //public static String libDir = "/usr/lib/GRIDD/lib"; // RELEASE
    
    // path to the res directory.
    // development: set this as current directory and place database/ folder inside "dist/"
    // release: set this in accordance to the structure of debian package. for example "/usr/share/gridd/"
    public final static String resDir = System.getProperty("user.dir") + "/"; // DEV
    //public static String resDir = "/usr/share/gridd/"; // RELEASE
    
        
    public static boolean checkBedtoolsNewVersion () {
        try {
            String[] args = {"bedtools", "--version"};
            String version = readBedtoolsVersion(args);
            if (version!=null){
                version = version.split("v")[1];
                String[] vnumber = version.split("\\.");
                if ( (Integer.parseInt(vnumber[0]) > 2) || ((Integer.parseInt(vnumber[0]) == 2 ) && (Integer.parseInt(vnumber[1])>23)) ) {
                    return true;
                } else {
                    return false;
                }
            } else {
                return false;
            }
        } catch (IOException ex) {
            Logger.getLogger(Commons.class.getName()).log(Level.SEVERE, null, ex);
        } catch (InterruptedException ex) {
            Logger.getLogger(Commons.class.getName()).log(Level.SEVERE, null, ex);
        }
        return false;
    }
    
    public static boolean checkSamtoolsNewVersion () {
        try {
            String[] args = {"samtools"};
            String version = readSamtoolsVersion(args);
            if (version!=null){
                version = version.split(" ")[1];
                String[] vnumber = version.split("-|\\.");
                if ( (Integer.parseInt(vnumber[0]) == 0) && (Integer.parseInt(vnumber[1]) <= 1 ) && (Integer.parseInt(vnumber[2])<=19) ) {
                    return false;
                } else {
                    return true;
                }
            } else {
                return false;
            }
        } catch (IOException ex) {
            Logger.getLogger(Commons.class.getName()).log(Level.SEVERE, null, ex);
        } catch (InterruptedException ex) {
            Logger.getLogger(Commons.class.getName()).log(Level.SEVERE, null, ex);
        } catch (ArrayIndexOutOfBoundsException ex) {
            
        }
        return false;
    }

    public static String readBedtoolsVersion(String args[]) throws IOException, InterruptedException {
        ProcessBuilder builder = new ProcessBuilder(args);
        Process p = builder.start();
        p.waitFor();
        BufferedReader br = null;
        try {
            br = new BufferedReader(new InputStreamReader(p.getInputStream()));
            String line = null;
            if ((line = br.readLine()) != null) {
                return line;
            }
        } catch (IOException ex) {
            Logger.getLogger(Commons.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            br.close();
        }
        return null;
    }
    
    public static String readSamtoolsVersion(String args[]) throws IOException, InterruptedException {
        ProcessBuilder builder = new ProcessBuilder(args);
        builder.redirectErrorStream(true);
        Process p = builder.start();
        p.waitFor();
        BufferedReader br = null;
        try {
            br = new BufferedReader(new InputStreamReader(p.getInputStream()));
            String line = null;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("Version")) {
                    return line.split(":")[1];
                }
            }
        } catch (IOException ex) {
            Logger.getLogger(Commons.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            br.close();
        }
        return null;
    }

    
    public static boolean checkFileExists(String filename) {
        File f = new File(filename);
        if(f.exists() && !f.isDirectory()) { 
            return true;
        } else {
            return false;
        }
    }
    
}
