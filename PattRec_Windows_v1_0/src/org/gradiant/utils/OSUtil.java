package org.gradiant.utils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.swing.JOptionPane;

public class OSUtil {
	
	public enum Property {
		SERIALNUMBER,
		MODEL,
		MANUFACTURER,
		VERSION
	}
	
	public static boolean isVM() {
		Map<String, List<String>> wmiValues = getWMIValues();
		
		List<String> property = wmiValues.get(Property.MODEL.name());
		if (property!=null) {
			String model = property.get(1).toLowerCase();
			//JOptionPane.showMessageDialog(null, model, "El Model del HW es: ", JOptionPane.INFORMATION_MESSAGE);
			if (model.contains("virtual machine") || model.contains("vmware") || model.contains("virtualbox") || model.contains("hvm domu") || model.contains("kvm")) {
				return true;
			}
		}
		
		property = wmiValues.get(Property.MANUFACTURER.name());
		if (property!=null) {
			String manufacturer = property.get(1).toLowerCase();
			//JOptionPane.showMessageDialog(null, manufacturer, "El Manufacturer del HW es: ", JOptionPane.INFORMATION_MESSAGE);
			if (manufacturer.contains("microsoft") || manufacturer.contains("vmware") || manufacturer.contains("innotek gmbh") || manufacturer.contains("xen") || manufacturer.contains("red hat")) {
				return true;
			}
		}
		
		property = wmiValues.get(Property.SERIALNUMBER.name());
		if (property!=null) {
			String serialNumber = property.get(1).toLowerCase();
			//JOptionPane.showMessageDialog(null, serialNumber, "El Serial Number de la Bios es: ", JOptionPane.INFORMATION_MESSAGE);
			if (serialNumber.contains("vmware")) {
				return true;
			}
		}
		
		property = wmiValues.get(Property.VERSION.name());
		if (property!=null) {
			String version = property.get(1).toLowerCase();
			//JOptionPane.showMessageDialog(null, version, "El Version de la Bios es: ", JOptionPane.INFORMATION_MESSAGE);
			if (version.contains("vrtual") || version.contains("version") || version.contains("a m i") || version.contains("xen")) {
				return true;
			}
		}
		
		return false;
	}
	
	private static Map<String, List<String>> getWMIValues() {		
		Map<String, List<String>> wmiValueList = new HashMap<>();
		try {		
	        Process process =Runtime.getRuntime().exec("cmd /c WMIC BIOS GET " + Property.SERIALNUMBER.name());
	        List <String> propertyValue = getPropertyValue(Property.SERIALNUMBER.name(), process);
	        wmiValueList.put(Property.SERIALNUMBER.name(), propertyValue);
	        
	        process =Runtime.getRuntime().exec("cmd /c WMIC COMPUTERSYSTEM GET " + Property.MODEL.name());
	        propertyValue = getPropertyValue(Property.MODEL.name(), process);
	        wmiValueList.put(Property.MODEL.name(), propertyValue);
	        
	        process =Runtime.getRuntime().exec("cmd /c WMIC COMPUTERSYSTEM GET " + Property.MANUFACTURER.name());
	        propertyValue = getPropertyValue(Property.MANUFACTURER.name(), process);
	        wmiValueList.put(Property.MANUFACTURER.name(), propertyValue);
	        
	        process =Runtime.getRuntime().exec("cmd /c WMIC WMIC BIOS GET " + Property.VERSION.name());
	        propertyValue = getPropertyValue(Property.VERSION.name(), process);
	        wmiValueList.put(Property.VERSION.name(), propertyValue);
	        
	    } catch (Exception e) {
	        e.printStackTrace();
	    }
		//Values of return map can be null
		return wmiValueList;		
	}
	
	private static List<String> getPropertyValue(String property, Process process) throws InterruptedException, IOException {	
		List<String> result = new ArrayList<>();
		process.waitFor();
        BufferedReader reader=new BufferedReader(
                new InputStreamReader(process.getInputStream())
                );
        String line;
        while((line = reader.readLine()) != null) {
            if(line != null && !line.trim().isEmpty()) {
            	result.add(line);
            }
        }
        if (result.size() < 2) {
        	return null;
        }
        return result;		
	}
	
}
