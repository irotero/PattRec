/*
 * Copyright (c) 2016 GRADIANT. All rights reserved.
 * This code cannot be used, copied, modified and/or distributed without the express permission of the authors.
 * This algorithm is protected under a Confidentiality Agreement between the UDTEMC (Fundación Ramón Domínguez) and GRADIANT.
 */
package org.gradiant.UI;

import java.awt.Color;
import java.awt.Image;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Date;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.imageio.ImageIO;
import javax.swing.JFileChooser;
import javax.swing.JTextArea;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.filechooser.FileNameExtensionFilter;
import org.gradiant.utils.Credential;
import org.gradiant.utils.CurrentDate;

import com.google.common.io.Files;

import org.gradiant.threads.NotifyingThread;
import org.gradiant.threads.R_Thread;
import org.gradiant.threads.ThreadCompleteListener;
import org.gradiant.utils.Commons;

public class RCNV extends JFrame implements ThreadCompleteListener {
    
    private File file_test;
    private String filename_BED;
    private String filename_fasta;
    private String filename_vcf;
    private String filename_poly = Commons.defaultPolyFile;
    private Credential cred;
    
    private File[] files_control = new File[20];
    private double per_dup = 0.3, per_del = 0.35, per_gen_dup = 0.3, per_gen_del = 0.35;
    private int cov_cont = 50;
    private boolean genes, fvcf, polymorphic, do_excel, do_restr, do_fixed, do_down, do_plot;


    public double getPer_dup() {
        return per_dup;
    }

    public double getPer_del() {
        return per_del;
    }

    public double getPer_gen_dup() {
        return per_gen_dup;
    }

    public double getPer_gen_del() {
        return per_gen_del;
    }

    public int getCov_cont() {
        return cov_cont;
    }

    public boolean isGenes() {
        return genes;
    }

    public boolean isFvcf() {
        return fvcf;
    }

    public boolean isPolymorphic() {
        return polymorphic;
    }
    
    public boolean isDo_excel() {
        return do_excel;
    }
    
    public boolean isDo_restr() {
        return do_restr;
    }
    
    public boolean isDo_fixed() {
        return do_fixed;
    }
    
    public boolean isDo_down() {
        return do_down;
    }
    
    public boolean isDo_plot() {
        return do_plot;
    }

    public JTextArea getOutput() {
        return output;
    }

    public File getFile_Test() {
        return file_test;
    }

    public String getFilename_BED() {
        return filename_BED;
    }

    public File[] getFiles_control() {
        return files_control;
    }
    
    public String getFilename_fasta() {
        return filename_fasta;
    }

    public String getFilename_vcf() {
        return filename_vcf;
    }
    
    public String getFilename_poly() {
        return filename_poly;
    }

    public Credential getCred() {
        return cred;
    }

    public RCNV(File file_test, String filename_BED, File[] files_control, String filename_fasta, Credential cred) {

        String mostrar;
        initInterface();
        
        this.cred = cred;

        this.file_test = file_test;
        mostrar = "Test: \n" + file_test.getAbsolutePath() + "\n";

        this.filename_BED = filename_BED;
        mostrar += "BedFile: \n" + filename_BED + "\n";

        this.files_control = files_control;
        mostrar += "Control:\n";
        
        for (File controlFile : files_control) {
            mostrar += controlFile + "\n";
        }   
        
        this.filename_fasta = filename_fasta;
        mostrar += "Fasta:\n" + filename_fasta + "\n";

        files.setText(mostrar);
    }
    
    private RCNV() {
        
    }

    private void initInterface() {
        initComponents();
        getContentPane().setBackground(new Color(239,243,255));
        setLocationRelativeTo(null);
        setResizable(true);
        setTitle("PattRec");
        Image icon;
        try {
        	icon = ImageIO.read(new FileInputStream(Commons.resDir + "res" + File.separator + "icons" + File.separator + "gridd.png"));
            setIconImage(icon);
        } catch (IOException ex) {
            Logger.getLogger(Interface.class.getName()).log(Level.WARNING, null, ex);
        }
    }
    
    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jCheckBoxMenuItem1 = new javax.swing.JCheckBoxMenuItem();
        jCheckBoxMenuItem2 = new javax.swing.JCheckBoxMenuItem();
        jCheckBoxMenuItem3 = new javax.swing.JCheckBoxMenuItem();
        jScrollPane3 = new javax.swing.JScrollPane();
        jPanel1 = new javax.swing.JPanel();
        jLabel1 = new javax.swing.JLabel();
        jLabel2 = new javax.swing.JLabel();
        jScrollPane1 = new javax.swing.JScrollPane();
        files = new javax.swing.JTextArea();
        minDup = new javax.swing.JTextField();
        minDel = new javax.swing.JTextField();
        minGenDup = new javax.swing.JTextField();
        minGenDel = new javax.swing.JTextField();
        covCont = new javax.swing.JTextField();
        useVCF = new javax.swing.JCheckBox();
        genAnalysis = new javax.swing.JCheckBox();
        excel = new javax.swing.JCheckBox();
        restr = new javax.swing.JCheckBox();
        fixed = new javax.swing.JCheckBox();
        down = new javax.swing.JCheckBox();
        plot = new javax.swing.JCheckBox();
        correct = new javax.swing.JButton();
        proceed = new javax.swing.JButton();
        jScrollPane2 = new javax.swing.JScrollPane();
        output = new javax.swing.JTextArea();
        jLabel4 = new javax.swing.JLabel();
        jLabel5 = new javax.swing.JLabel();
        jLabel8 = new javax.swing.JLabel();
        jLabel6 = new javax.swing.JLabel();
        jLabel7 = new javax.swing.JLabel();
        fieldVCF = new javax.swing.JTextField();
        findVCF = new javax.swing.JButton();
        ignorePolymorphic = new javax.swing.JCheckBox();
        fieldPolymorphic = new javax.swing.JTextField();
        findPolymorphic = new javax.swing.JButton();

        jCheckBoxMenuItem1.setSelected(true);
        jCheckBoxMenuItem1.setText("jCheckBoxMenuItem1");

        jCheckBoxMenuItem2.setSelected(true);
        jCheckBoxMenuItem2.setText("jCheckBoxMenuItem2");

        jCheckBoxMenuItem3.setSelected(true);
        jCheckBoxMenuItem3.setText("jCheckBoxMenuItem3");

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
        setBackground(new java.awt.Color(243, 182, 120));

        jScrollPane3.setBorder(null);
        jScrollPane3.setAlignmentX(0.0F);
        jScrollPane3.setAlignmentY(0.0F);
        jScrollPane3.setPreferredSize(new java.awt.Dimension(800, 690));

        jPanel1.setBackground(new java.awt.Color(239, 243, 255));
        jPanel1.setBorder(null);
        jPanel1.setAlignmentX(0.0F);
        jPanel1.setAlignmentY(0.0F);
        jPanel1.setPreferredSize(new java.awt.Dimension(780, 660));

        jLabel1.setFont(new java.awt.Font("Ubuntu", 0, 20)); // NOI18N
        jLabel1.setForeground(new java.awt.Color(19, 21, 136));
        jLabel1.setText("Options");

        jLabel2.setFont(new java.awt.Font("Ubuntu", 1, 14)); // NOI18N
        jLabel2.setForeground(new java.awt.Color(63, 126, 168));
        jLabel2.setText("FILES");

        files.setEditable(false);
        files.setColumns(20);
        files.setRows(5);
        jScrollPane1.setViewportView(files);

        minDup.setText("0.3");

        minDel.setText("0.35");

        minGenDup.setText("0.3");

        minGenDel.setText("0.35");

        covCont.setText("50");

        useVCF.setText("USE VCF");
        useVCF.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                useVCFActionPerformed(evt);
            }
        });

        genAnalysis.setSelected(true);
        genAnalysis.setText("GENE ANALYSIS");
        genAnalysis.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                genAnalysisActionPerformed(evt);
            }
        });

        excel.setSelected(true);
        excel.setText("XLSX OUTPUT");

        restr.setSelected(true);
        restr.setText("RESTRICTIVE");

        fixed.setSelected(true);
        fixed.setText("FIXED CONTROL SAMPLES");

        down.setSelected(false);
        down.setText("DOWNSAMPLING");

        plot.setSelected(false);
        plot.setText("PLOT");

        correct.setText("Correct files");
        correct.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                correctActionPerformed(evt);
            }
        });

        proceed.setText("Proceed");
        proceed.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                proceedActionPerformed(evt);
            }
        });

        output.setEditable(false);
        output.setColumns(20);
        output.setRows(5);
        jScrollPane2.setViewportView(output);

        jLabel4.setText("MIN DUP (0-1)");

        jLabel5.setText("MIN DEL (0-1)");

        jLabel8.setText("MIN CONT COV (>0)");

        jLabel6.setText("MIN GENE DUP (0-1)");

        jLabel7.setText("MIN GENE DEL (0-1)");

        fieldVCF.setEditable(false);
        fieldVCF.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                fieldVCFActionPerformed(evt);
            }
        });

        findVCF.setText("find");
        findVCF.setEnabled(false);
        findVCF.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                findVCFActionPerformed(evt);
            }
        });

        ignorePolymorphic.setText("IGNORE POLYMORPHIC");
        ignorePolymorphic.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                ignorePolymorphicActionPerformed(evt);
            }
        });

        fieldPolymorphic.setEditable(false);
        fieldPolymorphic.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                fieldPolymorphicActionPerformed(evt);
            }
        });

        findPolymorphic.setText("find");
        findPolymorphic.setEnabled(false);
        findPolymorphic.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                findPolymorphicActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout jPanel1Layout = new javax.swing.GroupLayout(jPanel1);
        jPanel1.setLayout(jPanel1Layout);
        jPanel1Layout.setHorizontalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, jPanel1Layout.createSequentialGroup()
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addComponent(jLabel1)
                .addGap(350, 350, 350))
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addGroup(javax.swing.GroupLayout.Alignment.LEADING, jPanel1Layout.createSequentialGroup()
                        .addGap(20, 20, 20)
                        .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addComponent(jScrollPane2, javax.swing.GroupLayout.PREFERRED_SIZE, 675, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addGroup(jPanel1Layout.createSequentialGroup()
                                .addComponent(jLabel2)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                                    .addGroup(jPanel1Layout.createSequentialGroup()
                                        .addGap(32, 32, 32)
                                        .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                            .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                                .addComponent(jLabel6, javax.swing.GroupLayout.Alignment.TRAILING)
                                                .addComponent(jLabel7))
                                            .addComponent(jLabel8)
                                            .addComponent(jLabel5)
                                            .addComponent(jLabel4))
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                                            .addComponent(minGenDup, javax.swing.GroupLayout.Alignment.LEADING)
                                            .addComponent(minDel, javax.swing.GroupLayout.Alignment.LEADING)
                                            .addComponent(covCont, javax.swing.GroupLayout.Alignment.LEADING)
                                            .addComponent(minDup, javax.swing.GroupLayout.Alignment.LEADING)
                                            .addComponent(minGenDel, javax.swing.GroupLayout.PREFERRED_SIZE, 69, javax.swing.GroupLayout.PREFERRED_SIZE))
                                        .addGap(122, 122, 122)
                                        .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                            .addComponent(excel)
                                            .addComponent(genAnalysis)
                                            .addComponent(restr)
                                            .addComponent(fixed)
                                            .addComponent(down)
                                            .addComponent(plot)))
                                    .addGroup(jPanel1Layout.createSequentialGroup()
                                        .addGap(195, 195, 195)
                                        .addComponent(fieldVCF, javax.swing.GroupLayout.PREFERRED_SIZE, 410, javax.swing.GroupLayout.PREFERRED_SIZE)
                                        .addGap(18, 18, 18)
                                        .addComponent(findVCF, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                                    .addComponent(jScrollPane1, javax.swing.GroupLayout.PREFERRED_SIZE, 675, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addGroup(jPanel1Layout.createSequentialGroup()
                                        .addComponent(correct, javax.swing.GroupLayout.PREFERRED_SIZE, 130, javax.swing.GroupLayout.PREFERRED_SIZE)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                        .addComponent(proceed, javax.swing.GroupLayout.PREFERRED_SIZE, 130, javax.swing.GroupLayout.PREFERRED_SIZE))))))
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addContainerGap()
                        .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(useVCF)
                            .addGroup(jPanel1Layout.createSequentialGroup()
                                .addComponent(ignorePolymorphic)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                .addComponent(fieldPolymorphic, javax.swing.GroupLayout.PREFERRED_SIZE, 410, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addGap(18, 18, 18)
                                .addComponent(findPolymorphic, javax.swing.GroupLayout.PREFERRED_SIZE, 52, javax.swing.GroupLayout.PREFERRED_SIZE)))))
                .addContainerGap(55, Short.MAX_VALUE))
        );
        jPanel1Layout.setVerticalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addGap(24, 24, 24)
                .addComponent(jLabel1)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jLabel2)
                    .addComponent(jScrollPane1, javax.swing.GroupLayout.PREFERRED_SIZE, 120, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(fieldVCF, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(findVCF)
                    .addComponent(useVCF))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(fieldPolymorphic, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(findPolymorphic)
                    .addComponent(ignorePolymorphic))
                .addGap(18, 18, 18)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(minDup, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(jLabel4))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(minDel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(jLabel5))
                        .addGap(6, 6, 6)
                        .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(covCont, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(jLabel8))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(minGenDup, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(jLabel6))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(minGenDel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(jLabel7))
                        .addGap(21, 21, 21)
                        .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(correct)
                            .addComponent(proceed))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(jScrollPane2, javax.swing.GroupLayout.PREFERRED_SIZE, 141, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addComponent(excel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(genAnalysis)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(restr)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(fixed)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(down)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(plot)))
                .addContainerGap(33, Short.MAX_VALUE))
        );

        jScrollPane3.setViewportView(jPanel1);

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addComponent(jScrollPane3, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addGap(0, 0, Short.MAX_VALUE))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addComponent(jScrollPane3, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addGap(0, 0, Short.MAX_VALUE))
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    // Correct files button
    private void correctActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_correctActionPerformed
        Interface obj = new Interface(file_test, filename_BED, files_control, filename_fasta);
        obj.setVisible(true);
        dispose();
    }//GEN-LAST:event_correctActionPerformed

    // Proceed button
    private void proceedActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_proceedActionPerformed
    	try {
    	    Date date = CurrentDate.getNTPDate();
        	if(date==null) {
            	System.exit(1);
        	}
    	    CurrentDate.compareExpirationDate(date);
    	} catch (Exception e) {
			JOptionPane.showMessageDialog(null, "Unknown error", "Can not launch application", JOptionPane.ERROR_MESSAGE);
        	System.exit(1);
    	}
        genes = genAnalysis.isSelected();
        fvcf = useVCF.isSelected();
        polymorphic = ignorePolymorphic.isSelected();
        do_excel = excel.isSelected();
        do_restr = restr.isSelected();
        do_fixed = fixed.isSelected();
        do_down = down.isSelected();
        do_plot = plot.isSelected();

        if (checkValues()) {
            proceed.setEnabled(false);
            NotifyingThread thread1 = new R_Thread(this);
            thread1.addListener(this); // add ourselves as a listener
            (new Thread(thread1)).start();
        } else {
            output.append("ERROR: input values are not valid\n");
        }
    }//GEN-LAST:event_proceedActionPerformed

    private void findVCFActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_findVCFActionPerformed
        JFileChooser chooser = new JFileChooser();
        FileNameExtensionFilter filter = new FileNameExtensionFilter("VCF FILES", "vcf");
        chooser.setFileFilter(filter);
        chooser.showOpenDialog(null);

        if(chooser.getSelectedFile() != null) {
            filename_vcf = chooser.getSelectedFile().getAbsolutePath();
            
            if(filename_vcf.endsWith(".vcf") == false) {
                filename_vcf = null;
                fieldVCF.setText("ERROR: selected file is not a vcf file");
            } else {
                fieldVCF.setText(filename_vcf);
            }
        }
    }//GEN-LAST:event_findVCFActionPerformed

    private void useVCFActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_useVCFActionPerformed
        findVCF.setEnabled(useVCF.isSelected());
    }//GEN-LAST:event_useVCFActionPerformed

    //Funci�n modificada para ser independiente del SO
    private void ignorePolymorphicActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_ignorePolymorphicActionPerformed
        findPolymorphic.setEnabled(ignorePolymorphic.isSelected());
        if (fieldPolymorphic.getText() == null || fieldPolymorphic.getText().equals("")) {
            fieldPolymorphic.setText(Commons.defaultPolyFile);
        }
        
        try {
            if (!Commons.checkFileExists(Commons.defaultPolyFile)) {
            	File destinyfolder = new File(Commons.defaultPolyFile);
            	File originfolder= new File(Commons.resDir + "res" + File.separator + "input" + File.separator + "polymorphic.bed");
            	if (originfolder.exists()) {
            		Files.copy(originfolder, destinyfolder);
            	}
                //Runtime.getRuntime().exec("cp " + Commons.resDir + "res/input/polymorphic.bed " + Commons.defaultPolyFile).waitFor();
            }
        } catch (IOException ex) {
            Logger.getLogger(RCNV.class.getName()).log(Level.SEVERE, null, ex);
        } //catch (InterruptedException ex) {
            //Logger.getLogger(RCNV.class.getName()).log(Level.SEVERE, null, ex);
        //}
        
    }//GEN-LAST:event_ignorePolymorphicActionPerformed

    private void findPolymorphicActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_findPolymorphicActionPerformed
        JFileChooser chooser = new JFileChooser();
        FileNameExtensionFilter filter = new FileNameExtensionFilter("BED FILES", "bed");
        chooser.setFileFilter(filter);
        chooser.showOpenDialog(null);

        if(chooser.getSelectedFile() != null) {
            filename_poly = chooser.getSelectedFile().getAbsolutePath();
            
            if(filename_poly.endsWith(".bed") == false) {
                filename_poly = null;
                fieldPolymorphic.setText("ERROR: selected file is not a bed file");
            } else {
                fieldPolymorphic.setText(filename_poly);
            }
        }
    }//GEN-LAST:event_findPolymorphicActionPerformed

    private void fieldPolymorphicActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_fieldPolymorphicActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_fieldPolymorphicActionPerformed

    private void fieldVCFActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_fieldVCFActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_fieldVCFActionPerformed

    private void genAnalysisActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_genAnalysisActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_genAnalysisActionPerformed


    private boolean checkValues() {
        try {
            per_dup = Double.parseDouble(minDup.getText());
            per_del = Double.parseDouble(minDel.getText());
            per_gen_dup = Double.parseDouble(minGenDup.getText());
            per_gen_del = Double.parseDouble(minGenDel.getText());
            cov_cont = Integer.parseInt(covCont.getText());
            double pers[] = {per_dup, per_del, per_gen_dup, per_gen_del};

            if(!checkPercentages(pers))
                throw new Exception();
            if(cov_cont < 1)
                throw new Exception();

        } catch(Exception e) {
            per_dup = 0.3;
            minDup.setText("0.3");
            per_del = 0.35;
            minDel.setText("0.35");
            per_gen_dup = 0.3;
            minGenDup.setText("0.3");
            per_gen_del = 0.35;
            minGenDel.setText("0.35");
            cov_cont = 50;
            covCont.setText("50");

            return false;
        }

        return true;
    }

    private static boolean checkPercentages(double pers[]) {
        for(double per : pers) {
            if(!checkPercentage(per))
                return false;
        }

        return true;
    }

    private static boolean checkPercentage(double per) {
        return (per > 0.0 && per <= 1.0);
    }
    /**
     * @param args the command line arguments
     */
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
            java.util.logging.Logger.getLogger(RCNV.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (InstantiationException ex) {
            java.util.logging.Logger.getLogger(RCNV.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            java.util.logging.Logger.getLogger(RCNV.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(RCNV.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>
        //</editor-fold>

        /* Create and display the form */
        java.awt.EventQueue.invokeLater(() -> {
            new RCNV().setVisible(true);
        });
    }

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton correct;
    private javax.swing.JTextField covCont;
    private javax.swing.JCheckBox excel;
    private javax.swing.JCheckBox restr;
    private javax.swing.JCheckBox fixed;
    private javax.swing.JCheckBox down;
    private javax.swing.JCheckBox plot;
    private javax.swing.JTextField fieldPolymorphic;
    private javax.swing.JTextField fieldVCF;
    private javax.swing.JTextArea files;
    private javax.swing.JButton findPolymorphic;
    private javax.swing.JButton findVCF;
    private javax.swing.JCheckBox genAnalysis;
    private javax.swing.JCheckBox ignorePolymorphic;
    private javax.swing.JCheckBoxMenuItem jCheckBoxMenuItem1;
    private javax.swing.JCheckBoxMenuItem jCheckBoxMenuItem2;
    private javax.swing.JCheckBoxMenuItem jCheckBoxMenuItem3;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JLabel jLabel5;
    private javax.swing.JLabel jLabel6;
    private javax.swing.JLabel jLabel7;
    private javax.swing.JLabel jLabel8;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JScrollPane jScrollPane2;
    private javax.swing.JScrollPane jScrollPane3;
    private javax.swing.JTextField minDel;
    private javax.swing.JTextField minDup;
    private javax.swing.JTextField minGenDel;
    private javax.swing.JTextField minGenDup;
    private javax.swing.JTextArea output;
    private javax.swing.JButton proceed;
    private javax.swing.JCheckBox useVCF;
    // End of variables declaration//GEN-END:variables

    @Override
    public void notifyOfThreadComplete(Runnable runnable) {
        proceed.setEnabled(true);
    }
}
