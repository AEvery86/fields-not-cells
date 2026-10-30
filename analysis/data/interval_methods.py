# -*- coding: utf-8 -*-
import numpy as np
from scipy.stats import t as tdist

def wilson(y,n):
 z=1.959963984540054;p=y/n;den=1+z*z/n
 mid=(p+z*z/(2*n))/den;half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
 return np.clip(mid-half,0,1),np.clip(mid+half,0,1)

def region_interval(sample,M):
 k=len(sample);ns=sample.sum(1);phat=sample.sum(0)/ns.sum()
 residual=sample-ns[:,None]*phat
 var=(1-k/M)*np.sum(residual**2,axis=0)/(k-1)/(k*ns.mean()**2)
 half=tdist.ppf(.975,k-1)*np.sqrt(np.maximum(var,0))
 return np.maximum(phat-half,0),np.minimum(phat+half,1)
