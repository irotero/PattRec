/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.gradiant.threads;

/**
 *
 * @author lgonzalez
 */
public interface ThreadCompleteListener {
    void notifyOfThreadComplete(final Runnable runnable);
}
