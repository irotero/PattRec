/*
 * Copyright (c) 2016 GRADIANT. All rights reserved.
 * This code cannot be used, copied, modified and/or distributed without the express permission of the authors.
 * This algorithm is protected under a Confidentiality Agreement between the UDTEMC (Fundación Ramón Domínguez) and GRADIANT.
 */
package org.gradiant.threads;

import java.util.Set;
import java.util.concurrent.CopyOnWriteArraySet;

public abstract class NotifyingThread implements Runnable{
    
  private final Set<ThreadCompleteListener> listeners
                   = new CopyOnWriteArraySet<>();
  
  public final void addListener(final ThreadCompleteListener listener) {
    listeners.add(listener);
  }
  
  public final void removeListener(final ThreadCompleteListener listener) {
    listeners.remove(listener);
  }
  
  private final void notifyListeners() {
    for (ThreadCompleteListener listener : listeners) {
      listener.notifyOfThreadComplete(this);
    }
  }
  
  @Override
  public final void run() {
    try {
      doRun();
    } finally {
      notifyListeners();
    }
  }
  
  public abstract void doRun();
}
