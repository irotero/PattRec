package org.gradiant.utils;

import java.io.IOException;
import java.net.InetAddress;
import java.util.Calendar;
import java.util.Date;
import java.util.GregorianCalendar;

import javax.swing.JOptionPane;

import org.apache.commons.net.ntp.NTPUDPClient;
import org.apache.commons.net.ntp.TimeInfo;

public class CurrentDate {

	 public static Date getNTPDate() {

	    	String[] hosts = new String[]{
	        	"hora.roa.es", "hora.rediris.es",
	        	"ntp.i2t.ehu.es"};

	    	NTPUDPClient client = new NTPUDPClient();
	    	// We want to timeout if a response takes longer than 5 seconds
	    	client.setDefaultTimeout(10000);

	    	for (String host : hosts) {
	        	try {
	            	//JOptionPane.showMessageDialog(null, host, "Probando con el servidor NTP:", JOptionPane.INFORMATION_MESSAGE);
	            	InetAddress hostAddr = InetAddress.getByName(host);
	            	TimeInfo info = client.getTime(hostAddr);          	
	            	Date date = new Date(info.getMessage().getTransmitTimeStamp().getTime());
	    	    	client.close();
	            	return date;
	        	}
	        	catch (IOException e) {
	            	e.printStackTrace();
	        	}
	    	}
			JOptionPane.showMessageDialog(null, "Cannot communicate with server. Please, check your Internet connection", "Can not launch application", JOptionPane.ERROR_MESSAGE);
	    	client.close();
	    	return null;

		}

	public static void compareExpirationDate(Date date) {
		Date expirationDate =new GregorianCalendar(2020, Calendar.DECEMBER, 31).getTime();
		if (date.after(expirationDate)) {
			JOptionPane.showMessageDialog(null, "Your license has expired. Contact us for renewal options", "Can not launch application", JOptionPane.ERROR_MESSAGE);
			System.exit(1);
		}
	}
}
