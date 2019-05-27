/*
 * Copyright (c) 2016 GRADIANT. All rights reserved.
 * This code cannot be used, copied, modified and/or distributed without the express permission of the authors.
 * This algorithm is protected under a Confidentiality Agreement between the UDTEMC (Fundación Ramón Domínguez) and GRADIANT.
 */
package org.gradiant;

import org.gradiant.UI.Interface;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Date;
import java.util.Properties;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.swing.JOptionPane;

import org.gradiant.utils.Commons;
import org.gradiant.utils.CurrentDate;
import org.gradiant.utils.OSUtil;

import com.google.common.io.Files;




public class main {
    public static void main(String args[]) {
        /* Set the Nimbus look and feel */
        //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
        /* If Nimbus (introduced in Java SE 6) is not available, stay with the default look and feel.
         * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html 
         */
        try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("Nimbus".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (ClassNotFoundException ex) {
            java.util.logging.Logger.getLogger(Interface.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (InstantiationException ex) {
            java.util.logging.Logger.getLogger(Interface.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            java.util.logging.Logger.getLogger(Interface.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(Interface.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>

        /* Create and display the form */
        java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {  
            	
            	try {
//            	    Date date = CurrentDate.getNTPDate();
//                	if(date==null) {
//                    	System.exit(1);
//                	}
//            	    CurrentDate.compareExpirationDate(date);
            	    if (OSUtil.isVM()) {
            	    	JOptionPane.showMessageDialog(null, "This program can't be executed in a virtual machine", "Can not launch application", JOptionPane.ERROR_MESSAGE);
                    	System.exit(1);
            	    }
            	} catch (Exception e) {
        			JOptionPane.showMessageDialog(null, "Unknown error", "Can not launch application", JOptionPane.ERROR_MESSAGE);
                	System.exit(1);
            	}
            	//Cambiado todo el c�digo de execs por manejo del objeto File para ser independientes del SO (Unix)
               try {
                	File folder = new File(Commons.baseDir + File.separator + "input");
                	folder.mkdirs(); 
                    //Runtime.getRuntime().exec("mkdir -p " + Commons.baseDir).waitFor();
                    //Runtime.getRuntime().exec("mkdir -p " + Commons.baseDir + "/input").waitFor();
                	
                    if (!Commons.checkFileExists(Commons.defaultPolyFile)) { //if not exists polymorphic in C:\Users\cinfo\PattRec\input
                    	File destinyfolder = new File(Commons.defaultPolyFile);
                    	File originfolder= new File(Commons.resDir + "res" + File.separator + "input" + File.separator + "polymorphic.bed");
                    	if (originfolder.exists()) { //if exists polymorphic in C:\Users\cinfo\Desktop\xxx\pattrec_windows\res\input
                    		Files.copy(originfolder, destinyfolder);
                    	}
                    	
                        //Runtime.getRuntime().exec("cp " + Commons.resDir + "res/input/polymorphic.bed " + Commons.defaultPolyFile).waitFor();
                    	//Runtime.getRuntime().exec("cp " + Commons.resDir + "res/input/picard.jar " + Commons.defaultPicardJar).waitFor();
                    	//Runtime.getRuntime().exec("cp " + Commons.resDir + "res/input/install_packages.R " + Commons.defaultInstallpackages).waitFor();
                    }
                } catch (IOException ex) {
                    Logger.getLogger(main.class.getName()).log(Level.SEVERE, null, ex);
                }// catch (InterruptedException ex) {
                 //  Logger.getLogger(main.class.getName()).log(Level.SEVERE, null, ex);
                 //}
                
                new Interface().setVisible(true);
                
        }
        });}
}
