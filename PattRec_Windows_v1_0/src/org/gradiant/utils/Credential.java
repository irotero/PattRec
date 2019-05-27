/*
 * Copyright (c) 2016 GRADIANT. All rights reserved.
 * This code cannot be used, copied, modified and/or distributed without the express permission of the authors.
 * This algorithm is protected under a Confidentiality Agreement between the UDTEMC (Fundación Ramón Domínguez) and GRADIANT.
 */
package org.gradiant.utils;

import java.io.Serializable;

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
