/*
 * Copyright (c) 2016 GRADIANT. All rights reserved.
 * This code cannot be used, copied, modified and/or distributed without the express permission of the authors.
 * This algorithm is protected under a Confidentiality Agreement between the UDTEMC (Fundación Ramón Domínguez) and GRADIANT.
 */
package org.gradiant.utils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.io.File;

public class Commons {

    public final static String baseDir = System.getProperty("user.home") + File.separator + "PattRec" + File.separator;
    public final static String defaultPolyFile = baseDir + "input" + File.separator + "polymorphic.bed";
 //   public final static String defaultPicardJar = baseDir + "input/picard.jar";
 //   public final static String defaultInstallpackages = baseDir + "input/install_packages.R";
    
    // path to the lib directory.
    // development: set this as current directory and place lib/RCNV/ folder inside "dist/"
    // release: set this in accordance to the structure of debian package. for example "/usr/lib/GRIDD/"
     public final static String libDir = System.getProperty("user.dir") + File.separator; // DEV
    // public static String libDir = "/usr/lib/GRIDD/"; // RELEASE
    
    // path to the res directory.
    // development: set this as current directory and place database/ folder inside "dist/"
    // release: set this in accordance to the structure of debian package. for example "/usr/share/gridd/"
     public final static String resDir = System.getProperty("user.dir")  + File.separator; // DEV
   //  public static String resDir = "/usr/share/gridd/"; // RELEASE
    
    public static boolean checkBedtoolsNewVersion () throws IOException, InterruptedException {
        ProcessBuilder builder = new ProcessBuilder("bedtools", "--version");
        Process p = builder.start();
        p.waitFor();
        String[] version = readVersion (p.getInputStream()).split("\\.");
        if ( (Integer.parseInt(version[0]) > 2) || ((Integer.parseInt(version[0]) == 2 ) && (Integer.parseInt(version[1])>23)) ) {
            return true;
        } else {
            return false;
        }
    }

    public static String readVersion(InputStream inputStream) throws IOException {
        BufferedReader br = null;
        String vnumber = "0.0.0";
        try {
            br = new BufferedReader(new InputStreamReader(inputStream));
            String line = null;
            if ((line = br.readLine()) != null) {
                vnumber = line.split("v")[1];
            }
        } catch (IOException ex) {
            Logger.getLogger(Commons.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            br.close();
        }
        return vnumber;
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
