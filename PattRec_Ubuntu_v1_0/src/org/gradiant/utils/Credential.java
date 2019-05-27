/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.gradiant.utils;

import java.io.Serializable;

/**
 *
 * @author mpacheco
 */
public class Credential implements Serializable {
    private String user;
    private String password;

    public Credential(String user, String password) {
        this.user = user;
        this.password = password;
    }
    
    public String getUser() {
        return user;
    }

    public String getPassword() {
        return password;
    }
}
