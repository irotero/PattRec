/*
 * Copyright (c) 2016 GRADIANT. All rights reserved.
 * This code cannot be used, copied, modified and/or distributed without the express permission of the authors.
 * This algorithm is protected under a Confidentiality Agreement between the UDTEMC (Fundación Ramón Domínguez) and GRADIANT.
 */
package org.gradiant.JRI;


import org.rosuda.JRI.REXP;
import org.rosuda.JRI.RMainLoopCallbacks;
import org.rosuda.JRI.Rengine;

public class RDBScript  implements RMainLoopCallbacks {
    final private String rCode;
    Boolean print;
                
    public RDBScript () {

        this.print = false;
        this.rCode = "cnvs_database <- function(userDB, passwordDB, name_test, name_file) {\n" +
		"	list.of.packages <- c(\"RMySQL\", \"openxlsx\")\n" +
		"	new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,\"Package\"])]\n" +
		"	if(length(new.packages)) install.packages(new.packages)\n" +
		"	library(RMySQL)\n" +
		"	library(openxlsx)\n" +
		"	# Read data from result file:\n" +
		"	ifelse(grepl(\".xlsx\", name_file), data_aux <- read.xlsx(name_file, sheet=1), data_aux <- read.table(name_file, header = T, sep = \"\\t\", quote = \"\\\"\", dec = \",\"))\n" +
		"	#-------------------------------------------------------\n" +
		"	# Extract columns (chr, start, end, gen, type, %):\n" +
		"	if(nrow(data_aux)!=0){\n" +
		"		data_aux <- data_aux[data_aux$Chr!=\"\",]\n" +
		"		data <- data_aux[,1:5]\n" +
		"		data[,6] <- data_aux[,grepl(\"decrease\",names(data_aux))]\n" +
		"		# Add column for sample name:\n" +
		"		data <- as.data.frame(cbind(name_test,data))\n" +
		"		# Set column names of database:\n" +
		"		names(data) <- c(\"Sample\", \"Chr\", \"Start\", \"End\", \"Gene\", \"Type\", \"Percentage\")\n" +
		"		# Connect to mysql:\n" +
		"		dbPrueba<-dbConnect(MySQL(), user=userDB, password=passwordDB, host=\"localhost\", dbname=\"cnv\")\n" +
		"		consultaString = paste0(\"select * from cnv_global where Sample like \\\"\",name_test, \"\\\";\")\n" +
		"		consultaSql <- dbSendQuery(dbPrueba, consultaString)\n" +
		"		aux <-fetch(consultaSql, n=-1)\n" +
		"		if(nrow(aux)!=0){\n" +
		"			rpp <- do.call(paste0,data)%in%do.call(paste0,aux)\n" +
		"			if(TRUE%in%rpp) data <- data[rpp==FALSE,]\n" +
		"		}\n" +
		"		# Write to DB:\n" +
		"		if(nrow(data)!=0) dbWriteTable(dbPrueba, \"cnv_global\", data, append=TRUE, row.names=FALSE) \n" +
		"		dbDisconnect(dbPrueba)\n" +
		"	}\n" +
		"}";
	}
        
    public int execute (String userDB, String passwordDB, String name_test, String name_file) {
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
        re.assign("userDB", userDB);
        re.assign("passwordDB", passwordDB);
        re.assign("name_test", name_test);
        re.assign("name_file", name_file);
        re.eval(rCode);
        print=true;
        REXP result = re.eval("cnvs_database(userDB, passwordDB, name_test, name_file)");
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
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
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
