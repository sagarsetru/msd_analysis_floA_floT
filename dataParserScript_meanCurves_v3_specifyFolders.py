#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 14:54:40 2020

@author: sagarsetru

plot MSD analysis for flo domains
"""

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import itertools
import glob
import pickle

#%%
def openFileList(fListTxtFile):
    with open(fListTxtFile) as f:
        txtFileList = f.readlines()
    txtFileList = [x.strip() for x in txtFileList]
    txtFileList.sort()
    return txtFileList
        
#%% compute average MSD with 95% bootstrap confidence intervals

def computeAverageAndBootstrap(fileList):
    # initialize lists to return values
    tausReturn = {}
    avgsReturn = {}
    btstrpciReturn = {}
    btstrpci2Return = {}
    stdReturn = {}
    nReturn = {}
    
    # loop through FloA, FloT files
    counter1 = 0
    for mainFile in fileList:
        counter1 += 1
        print('working on: '+mainFile)
        print('file',counter1,'out of',len(fileList))
        print(' ')
        
        # get directory in which to save 
        saveDir = os.path.split(mainFile)[0]
        baseFileName_w_ext = os.path.split(mainFile)[1]
        baseFileName = os.path.splitext(baseFileName_w_ext)[0]
        
        # get list of individual csv files
        csvFileList = openFileList(mainFile)
        
        # variable for max track length
        trackLengthMax = 0;
        
        # list for all deltaTaus
        deltaTaus = [];
    
        # list for all MSDs
        #msd_all = [];
    
        # loop through each individual csv file
        counter_nFiles = 0
        for f in csvFileList:
            
            counter_nFiles += 1
            
            # read the data into df
            tau_and_val = pd.read_csv(f,header=None)
            
            # get the shape of the df
            dataShape = tau_and_val.shape
            
            # get the length of this tracks msd
            trackLength = dataShape[0]
            
            # ignore any tracks that are only 1 entry long
            if trackLength == 1:
                print('This track is only one entry long. Ignoring.')
                print(f)
                print(' ')
                continue
            
            # get the max tracklength
            if trackLength > trackLengthMax:
#                print(trackLength)
#                print(f)
                trackLengthMax = trackLength
            #...
            
            # get the mean tau value for this data set
            tauMean = tau_and_val[0].diff().mean()
            
#            if np.isnan(tauMean):
#                print(f)
#                break
            # add to list of taus
            deltaTaus.append(tauMean)
            
             # get the msds
    #        msd_all.append(tau_and_msd.iloc[:,1])
        #...
        
        # check to make sure all deltaTaus are the same; first round
        deltaTaus = np.around(deltaTaus,decimals=1)
        if not np.all(deltaTaus==deltaTaus[0]):
            
            print('Not all delta tau values are the same!')
            break
        
        else:
            
            # use this value for tau (x axis)
            taus = np.cumsum(deltaTaus[0:trackLengthMax])
            
        #...
        
        # np array for all values
        val_all = np.zeros((counter_nFiles,trackLengthMax))
        
        # convert to nan
        val_all[:] = np.nan
        
        # add each msd to msd_all array (or vac)
        counter = -1
        for f in csvFileList:
            counter += 1
            
            # sparsify
#            if counter % 10 == 1:
            
            # read the msd and tau into df
            tau_and_val = pd.read_csv(f,header=None)
            
            # get the shape of the df
            dataShape = tau_and_val.shape
            
            # get the length of this tracks msd
            trackLength = dataShape[0]
            
            # ignore any tracks that are only 1 entry long
            if trackLength == 1:
                print('This track is only one entry long. Ignoring.')
                print(f)
                print(' ')
                continue
            
#            if trackLength == 17:
#                print(f)
#                test=tau_and_val[1]
#                break
            
            # get msd or vac
            val = tau_and_val[1]
            
            # remove spurious zeros
#            if 'vac' in mainFile:
            val=val[val!=0]
            #...
            
            # store values np array
            val_all[counter][0:len(val)] = val
        #...
        
        # get average value of msd at each tau
        val_mean = np.nanmean(val_all,axis=0)
        
        # remove remaining nan values
        val_mean = val_mean[~np.isnan(val_mean)]
        
        # adjust taus to compensate for removing nan values from mean
        taus = taus[0:np.shape(val_mean)[0]]
        
        # get the standard deviation
        val_std = np.nanstd(val_all,axis=0)
        
        # remove remaining nan values
        val_std = val_std[~np.isnan(val_std)]
        
        # get the number of samples (for standard error of mean)
        val_n = np.count_nonzero(~np.isnan(val_all),axis=0)
        
        # remove those which only have all nans, and so will have zero samples
        val_n = val_n[val_n!=0]
        
        # for bootstrap confidence intervals
        val_btstrpci = np.zeros((2,np.size(val_mean)))  
        
        # get 95% bootstrap confidence intervals
        for i in range(0,np.size(val_mean)):
            #print(i)
            
            # get the values at this value of tau
            val = val_all[:,i]
            
            # remove nan values
            val = [x for x in val if not np.isnan(x)]
            
            # do bootstrap
            
            # define bootstrap parameters 
            n = len(val)
            reps = 10000
            
            # define bootstrap distribution
            btstrpDist = np.random.choice(val,(n,reps))
            
            # calculate mean of bootstrap distribution
            btstrpMean = np.mean(btstrpDist,axis=0)
            btstrpMean.sort()
            
            btstrp_ci = np.percentile(btstrpMean,[2.5,97.5])
            val_btstrpci[:,i] = btstrp_ci
        #...
        
        # make 2d array with taus and the average values
        val_mean = np.vstack((taus,val_mean.T)).T
        
        # make 2d array with taus and the bootstrap ci's
#        val_btstrpci = np.vstack((taus,val_btstrpci[0][:],val_btstrpci[1][:])).T
        val_btstrpci = np.vstack((taus,val_btstrpci)).T

        
        # for bootstrap conf ints across trajectories
        val_btstrpci_2 = np.zeros((2,np.shape(val_mean)[0]))
        
        # bootstrap parameters
        n_trajectories = np.shape(val_all)[0]
#        n_choose = 200
#        n_choose = int(n_trajectories/2)
        n_choose = n_trajectories
        reps = 10000
        
#        btstrpDist_2 = np.zeros((reps,np.shape(val_all)[1]))
        btstrpDist_2 = np.zeros((reps,np.shape(val_mean)[0]))
                
        # build bootstrap distribution
        for i in range(0,reps):
#            print(i)
            idx = np.random.randint(n_trajectories,size=n_choose)
            
#            btstrpDist_temp = val_all[idx,:]
            btstrpDist_temp = val_all[idx,0:np.shape(val_mean)[0]]
            
            # get average value of msd at each tau
            val_mean_temp = np.nanmean(btstrpDist_temp,axis=0)
            
            # add to btstrp dist
            btstrpDist_2[i,:] = val_mean_temp
        #...
        
        # get the mean from btstrp dist
#        btstrpDist_2_mean = np.mean(btstrpDist_2,axis=0)
        
        # get 95% bootstrap confidence intervals across trajectories
        for i in range(0,np.shape(btstrpDist_2)[1]):
            #print(i)
            
            # get the values at this value of tau
            val = btstrpDist_2[:,i]
            
            # sort the values
            val.sort()
            
            btstrp_ci = np.nanpercentile(val,[2.5,97.5])
            val_btstrpci_2[:,i] = btstrp_ci
        #...
        
        # make 2d array with taus and the bootstrap ci's
#        val_btstrpci_2 = np.vstack((taus,val_btstrpci_2[0][:],val_btstrpci_2[1][:])).T
        val_btstrpci_2 = np.vstack((taus,val_btstrpci_2)).T
        
        # get the standard deviation from the bootstrap distribution
#        val_std_2 = np.nanstd(btstrpDist_2,axis=0)
        
#        plt.plot(val_mean[:,1])
#        plt.plot(val_std_2)
#        plt.plot(val_std)
##        plt.plot(val_std/val_n)
#        plt.plot(val_btstrpci[:,1])
#        plt.plot(val_btstrpci[:,2])
##        plt.plot(val_btstrpci_2[:,1])
##        plt.plot(val_btstrpci_2[:,2])
#        plt.plot(val_all.T)
#        plt.xscale('log')
#        plt.yscale('log')
#        matchAll={}
#        for a in idx:
#            matches = idx == a
#            matchAll[a] = np.sum(matches)
        
        
        # save taus, averaged , and 95% bootstrap conf ints to directory
        if 'msd' in mainFile:
            
            np.savetxt(saveDir+'/'+baseFileName+'_'+'avg_msd.csv',val_mean,delimiter=',')
            #np.savetxt(saveDir+'/'+baseFileName+'_'+'avg_msd_taus.csv',taus,delimiter=',')
            np.savetxt(saveDir+'/'+baseFileName+'_'+'btstrp_95ci_msd.csv',val_btstrpci,delimiter=',')
            np.savetxt(saveDir+'/'+baseFileName+'_'+'btstrp_95ci_acrossTraj_msd.csv',val_btstrpci_2,delimiter=',')
            np.savetxt(saveDir+'/'+baseFileName+'_'+'btstrp_std_msd.csv',val_std,delimiter=',')
            np.savetxt(saveDir+'/'+baseFileName+'_'+'btstrp_nSamples_msd.csv',val_n,delimiter=',')
            
        elif 'vac' in mainFile:
            
            np.savetxt(saveDir+'/'+baseFileName+'_'+'avg_vac.csv',val_mean,delimiter=',')
            #np.savetxt(saveDir+'/'+baseFileName+'_'+'avg_vac_taus.csv',taus,delimiter=',')
            np.savetxt(saveDir+'/'+baseFileName+'_'+'btstrp_95ci_vac.csv',val_btstrpci,delimiter=',')
            np.savetxt(saveDir+'/'+baseFileName+'_'+'btstrp_95ci_acrossTraj_msd.csv',val_btstrpci_2,delimiter=',')
            np.savetxt(saveDir+'/'+baseFileName+'_'+'btstrp_std_vac.csv',val_std,delimiter=',')
            np.savetxt(saveDir+'/'+baseFileName+'_'+'btstrp_nSamples_vac.csv',val_n,delimiter=',')
            
        else:
            
            print('neither "vac" nor "msd" is in text file name of individual msds. fix!')
            break
        
        #...
        
        # plot each individual track
        
#        # axis font size
#        afs = 24
#        
#        # tick font size
#        tfs = 18
#        
#        # line width
#        lw = 2
#        
#        plt.plot(val_all.T, linewidth=lw)
#        plt.xscale('log')
#        plt.yscale('log')
#        plt.xlabel(r'$\tau$ (s)',fontsize=afs)
#        plt.xticks(fontsize=tfs)
#        plt.yticks(fontsize=tfs)
#        
#        if 'msd' in mainFile:
#            plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)
#            plt.savefig(saveDir+'/'+baseFileName+'_allMSD.pdf',bbox_inches="tight")
#
#        elif 'vac' in mainFile:
#            plt.ylabel(r'MSD (($\mu$m/s)$^2$)',fontsize=afs)
#            plt.savefig(saveDir+'/'+baseFileName+'_allVAC.pdf',bbox_inches="tight")
            
        #...
        
        

        
        # add values to list
        #tausReturn[mainFile] = taus
        avgsReturn[baseFileName] = val_mean
        btstrpciReturn[baseFileName] = val_btstrpci
        btstrpci2Return[baseFileName] = val_btstrpci_2
        stdReturn[baseFileName] = val_std
#        stdReturn[baseFileName] = val_std_2
        nReturn[baseFileName] = val_n
        
    # return tuple of taus, averaged values, and bootstrap confidence intervals
    return avgsReturn, btstrpciReturn, btstrpci2Return, stdReturn, nReturn
#..  
    
#%%
   
    
def returnListForPlotting(dictToUse,keysToUse,indMain):
    # get the main dictionary
    mainDictList = dictToUse[keysToUse[indMain]]
    
    # retrun relevant values
    if len(keysToUse[indMain+1]) == 2:
        
        # get sub keys (e.g., for FloT vs FloA)
        subKey0 = keysToUse[indMain+1][0]
        subKey1 = keysToUse[indMain+1][1]
        
        # get the dictionaries
        dict_vals = mainDictList[0]
        dict_cis = mainDictList[1]
        dict_cis2 = mainDictList[2]
        dict_stds = mainDictList[3]
        dict_nSamples = mainDictList[4]
        
        # put into list
        returnList = [
                dict_vals[subKey0],
                dict_cis[subKey0],
                dict_cis2[subKey0],
                dict_stds[subKey0],
                dict_nSamples[subKey0],
                subKey0,
                dict_vals[subKey1],
                dict_cis[subKey1],
                dict_cis2[subKey1],
                dict_stds[subKey1],
                dict_nSamples[subKey1],
                subKey1,
                ]
        
        # return the relevant quantities, with keys to use for labelling axes
        return returnList
    
    elif len(keysToUse[indMain+1]) == 1:
        
        # get sub keys (e.g., for FloT vs FloA)
        subKey0 = keysToUse[indMain+1][0]
        
        # get the dictionaries
        dict_vals = mainDictList[0]
        dict_cis = mainDictList[1]
        dict_cis2 = mainDictList[2]
        dict_stds = mainDictList[3]
        dict_nSamples = mainDictList[4]
        
        returnList = [
                dict_vals[subKey0],
                dict_cis[subKey0],
                dict_cis2[subKey0],
                dict_stds[subKey0],
                dict_nSamples[subKey0],
                subKey0,
                ]
                      
        
        # return the relevant quantities, with keys to use for labelling axes
        return returnList
    
    # retrun relevant values
    elif len(keysToUse[indMain+1]) == 9:
        
        # get sub keys (e.g., for FloT vs FloA)
        subKey0 = keysToUse[indMain+1][0]
        subKey1 = keysToUse[indMain+1][1]
        subKey2 = keysToUse[indMain+1][2]
        subKey3 = keysToUse[indMain+1][3]
        subKey4 = keysToUse[indMain+1][4]
        subKey5 = keysToUse[indMain+1][5]
        subKey6 = keysToUse[indMain+1][6]
        subKey7 = keysToUse[indMain+1][7]
        subKey8 = keysToUse[indMain+1][8]
        
        # get the dictionaries
        dict_vals = mainDictList[0]
        dict_cis = mainDictList[1]
        dict_cis2 = mainDictList[2]
        dict_stds = mainDictList[3]
        dict_nSamples = mainDictList[4]
        
        # put into list
        returnList = [
                dict_vals[subKey0],
                dict_cis[subKey0],
                dict_cis2[subKey0],
                dict_stds[subKey0],
                dict_nSamples[subKey0],
                subKey0,
                dict_vals[subKey1],
                dict_cis[subKey1],
                dict_cis2[subKey1],
                dict_stds[subKey1],
                dict_nSamples[subKey1],
                subKey1,
                dict_vals[subKey2],
                dict_cis[subKey2],
                dict_cis2[subKey2],
                dict_stds[subKey2],
                dict_nSamples[subKey2],
                subKey2,
                dict_vals[subKey3],
                dict_cis[subKey3],
                dict_cis2[subKey3],
                dict_stds[subKey3],
                dict_nSamples[subKey3],
                subKey3,
                dict_vals[subKey4],
                dict_cis[subKey4],
                dict_cis2[subKey4],
                dict_stds[subKey4],
                dict_nSamples[subKey4],
                subKey4,
                dict_vals[subKey5],
                dict_cis[subKey5],
                dict_cis2[subKey5],
                dict_stds[subKey5],
                dict_nSamples[subKey5],
                subKey5,
                dict_vals[subKey6],
                dict_cis[subKey6],
                dict_cis2[subKey6],
                dict_stds[subKey6],
                dict_nSamples[subKey6],
                subKey6,
                dict_vals[subKey7],
                dict_cis[subKey7],
                dict_cis2[subKey7],
                dict_stds[subKey7],
                dict_nSamples[subKey7],
                subKey7,
                dict_vals[subKey8],
                dict_cis[subKey8],
                dict_cis2[subKey8],
                dict_stds[subKey8],
                dict_nSamples[subKey8],
                subKey8,
                ]
        
        # return the relevant quantities, with keys to use for labelling axes
        return returnList
    
    else:
        
        print('This dictionary list has more than 1 or 2 or 9 (eg for DltD, PBP3) elements.')
        print('This means more than two things were tracked in this experiment.')
        print('Investigate!')
    #...
#...
    
    
#%%
    
def doPlotCondition(returnList,col,leg,lw,alphaVal,errorFormat):
    
    # get x values (tau)
    x = returnList[0][:,0]
    
    # get y values (msd or vac)
    y = returnList[0][:,1]
    
    yMin = np.min(y)
    yMax = np.max(y)
    
    yExt = np.array([yMin,yMax])
    
    plt.plot(x,y,col,label=leg,linewidth=lw)
    
    if errorFormat == 'ci':
        
        # get lower confidence intervals
        ci_l = returnList[1][:,1]
        
        # get upper confidence intervals
        ci_u = returnList[1][:,2]
        
        plt.fill_between(x, ci_l, ci_u, alpha=alphaVal,facecolor=col)
        
    elif errorFormat == 'ci2':
        
        # get lower confidence intervals (across trajectories)
        ci_l = returnList[2][:,1]
        
        # get upper confidence intervals (across trajectories)
        ci_u = returnList[2][:,2]
        
        plt.fill_between(x, ci_l, ci_u, alpha=alphaVal,facecolor=col)
        
    elif errorFormat == 'std':
        
        # get the standard deviation
        std = returnList[3]
    
        plt.fill_between(x, y-std, y+std,alpha=alphaVal,facecolor=col)
        
    elif errorFormat == 'sem':
        
        # get the sem
        sem = returnList[3]/np.sqrt(returnList[4])
    
        plt.fill_between(x, y-sem, y+sem,alpha=alphaVal,facecolor=col)
        
    else:
        
        print('errorFormat string should be "ci" or "ci2" or "std" or "sem"!')
    
    #...
    
    return yExt,x,y
    

#%%
#dirs_all = [
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/m+LBO2',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/m-LBO2',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mAmp',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mBnz',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mCtrlBNZ',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mEPI_MAD',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mEPI_WU',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mFos',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mNis',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mProtoplasts',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mCtrlProtoplasts',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mTIRF100-0.2',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mTIRF35-0.3',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mTun',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mVal',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mVan',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/m∆dlt',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/m∆mreB',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/m∆mreB:BH:mbl',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/m∆pbpC',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mCtrl∆pbpC',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/sDMSO',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/sTunWTA',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/spH8.1',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/s∆dlt+Mg',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/s∆dltA',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/s∆flotillin',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/s∆rsgI',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/s∆tagU',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/s∆ugtP',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/FloAntTct',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/FloTntAct',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/MreB',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/PBP3',
#]

#dirs_all = [
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/m+LBO2',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/m-LBO2',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mAmp',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mCtrlAmp',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mBnz',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mCtrlBNZ',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mEPI_MAD',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mEPI_WU',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mFos',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mCtrlFos',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mNis',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mCtrlNis',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mProtoplasts',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mCtrlProtoplasts',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mTIRF100-0.2',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mTIRF35-0.3',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mTun',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mCtrlTun',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/sCtrlTunWTA',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mVal',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mCtrlVal',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mVan',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mCtrlVan',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/m∆dlt',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mCtrl∆dlt',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/m∆mreB',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/m∆mreB:BH:mbl',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/s∆rsgI',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/s∆rsgI-2',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/sCtrl∆rsgI-2',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/m∆pbpC',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/mCtrl∆pbpC',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/sDMSO',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/sCtrlDMSO',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/sTunWTA',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/spH8.1',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/s∆dlt+Mg',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/s∆dltA',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/s∆flotillin',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/sCtrl∆flotillin',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/s∆tagU',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/sCtrl∆tagU',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/s∆ugtP',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/sCtrl∆utgP',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/FloAntTct',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/FloTntAct',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/DltD',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/PBP3-2',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/MreB',
#]

dirs_all = [
'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/Amp10min',
'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/Amp30min',
'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/Amp90min',
'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/Fos10min',
'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/Fos30min',
'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/Fos90min',
'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/Nis10min',
'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/Nis30min',
'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/Nis90min',
'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/Tun10min',
'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/Tun30min',
'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/Tun90min',
'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/Van10min',
'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/Van30min',
'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/Van90min',
'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/AmpCtrl',
'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/FosCtrl',
'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/NisCtrl',
'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/TunCtrl',
'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/VanCtrl',
        ]


#mainDir = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/PBP3'
#mainDir = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins'

# where to save misc files
baseDir = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis'
os.chdir(baseDir)
#%%

# dictionary for msds
# each value in dictionary is a list of length two, keyed by the experimental condition
# each element of list is a dictionary

# list[0] = avg vals; list[1] = 95% conf. ints.; list[2] = std dev; list[3] = nSamples;
# that is,
# the zeroth element of each list is the avg vals keyed by condition (eg FloA or FloT)
# the first element of each list is the 95% conf. ints. keyed by condition (eg FloA or FloT)
# the second elemenet is the std dev
# the third element is the nSamples 

msd_dict = {}
keys_msd_dict = []

# dictionary for vacs, follows same structure
vac_dict = {}
keys_vac_dict = []

counter = 0
for mainDir in dirs_all:
    
    counter += 1
    print('WORKING ON directory',counter,'out of',len(dirs_all))
#    print(mainDir)
    # get key to use for dictionary
    dict_key = os.path.split(mainDir)[1]
    print(dict_key)
    
    ## FOR VACS
    
    # update list of keys
    keys_vac_dict.append(dict_key)
    print('vac')
    
    # get msd file list
    fileList = glob.glob(mainDir+'/*vacTracks.txt')
    
    # get avg vacs and cis and std and nsamples
    avgs, cis, cis_traj, stds, nSamples = computeAverageAndBootstrap(fileList)
    
    # load msd dictionary
    vac_dict[dict_key] = [avgs,cis,cis_traj,stds,nSamples]
    
    # update list of keys with keys of internal dictionary
    keys_vac_dict.append(list(avgs))
    print(list(avgs))
    
    print(' ')
    
    ## FOR MSDS
    
    # update list of keys
    keys_msd_dict.append(dict_key)
    
    # get msd file list
    fileList = glob.glob(mainDir+'/*msdTracks.txt')
    print('msd')
    
    # get avg msds and cis and std and nsamples
    avgs, cis, cis_traj, stds, nSamples = computeAverageAndBootstrap(fileList)
    
    # load msd dictionary
    msd_dict[dict_key] = [avgs,cis,cis_traj,stds,nSamples]
    
    # update list of keys with keys of internal dictionary
    keys_msd_dict.append(list(avgs))
    print(list(avgs))
    

#...
    
#%% get msd file list for all FloA
    
fileList = ['/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/FloA_ctrl_msdFiles_list.txt']

f = open(fileList[0],"r")
contents = f.read().splitlines()
all_FloA = []
for c in contents:
    f2 = open(c,"r")
    contents2 = f2.read().splitlines()
    all_FloA.append(contents2)
#...

fileToWrite_FloA = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/FloA_ctrl_msd.txt'
all_FloA = list(itertools.chain.from_iterable(all_FloA))
with open(fileToWrite_FloA, 'w') as f:
    for item in all_FloA:
        f.write("%s\n" % item)

#%% get msd file list for all FloT
    
fileList = ['/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/FloT_ctrl_msdFiles_list.txt']

f = open(fileList[0],"r")
contents = f.read().splitlines()
all_FloT = []
for c in contents:
    f2 = open(c,"r")
    contents2 = f2.read().splitlines()
    all_FloT.append(contents2)
#...

fileToWrite_FloT = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/FloT_ctrl_msd.txt'
all_FloT = list(itertools.chain.from_iterable(all_FloT))
with open(fileToWrite_FloT, 'w') as f:
    for item in all_FloT:
        f.write("%s\n" % item)


#%% for all floa control
    
dict_key = 'FloA_allCtrl'

# update list of keys
keys_msd_dict.append(dict_key)

fileList = ['/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/FloA_ctrl_msd.txt']

# get avg msds and cis and std and nsamples
avgs, cis, cis_traj, stds, nSamples = computeAverageAndBootstrap(fileList)

# load msd dictionary
msd_dict[dict_key] = [avgs,cis,cis_traj,stds,nSamples]

# update list of keys with keys of internal dictionary
keys_msd_dict.append(list(avgs))
print(list(avgs))
    
#%% for all flot control

dict_key = 'FloT_allCtrl'

# update list of keys
keys_msd_dict.append(dict_key)

fileList = ['/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/Flotillins/FloT_ctrl_msd.txt']

# get avg msds and cis and std and nsamples
avgs, cis, cis_traj, stds, nSamples = computeAverageAndBootstrap(fileList)

# load msd dictionary
msd_dict[dict_key] = [avgs,cis,cis_traj,stds,nSamples]

# update list of keys with keys of internal dictionary
keys_msd_dict.append(list(avgs))
print(list(avgs))
    

#counter = 0
#for mainDir in dirs_all:
#    
#    counter += 1
#    print('WORKING ON directory',counter,'out of',len(dirs_all))
##    print(mainDir)
#    # get key to use for dictionary
#    dict_key = os.path.split(mainDir)[1]
#    print(dict_key)
#    
#    ## FOR MSDS
#    
#    # update list of keys
#    keys_msd_dict.append(dict_key)
#    
#    # get msd file list
#    fileList = glob.glob(mainDir+'/*msdTracks.txt')
#    print(fileList)
    
#%% save dictionaries and lists of keys via pickle
pklFile = baseDir+'msd_vac_dicts_keys.pkl'

with open(pklFile, 'wb') as f:
    
    pickle.dump([msd_dict,keys_msd_dict,vac_dict,keys_vac_dict],f)
#...
    
## DUMMY VARIABLES
#dict_msd = msd_dict
#keys_dict_msd = keys_msd_dict
#dict_vac = vac_dict
#keys_dict_vac = keys_vac_dict
#%%
txtDictFile = baseDir+'msd_data_all.txt'

#import json
#
#with open(txtDictFile,'w') as file:
#    file.write(json.dumps(msd_dict))
    
    
with open(txtDictFile, 'w') as f:
    for key, value in msd_dict.items():
        f.write('%s:%s\n' % (key, value[4]))
    
    
#%% TO RELOAD DICTIONARIES AND LISTS OF KEYS
baseDir = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis'

pklFile = baseDir+'msd_vac_dicts_keys.pkl'

with open(pklFile, 'rb') as f:  # Python 3: open(..., 'rb')
#    f.open(pklFile,'rb')
    msd_dict, keys_msd_dict, vac_dict, keys_vac_dict = pickle.load(f)
#...
    
#%% save dictionaries for antibiotic time course via pickle
pklFile = baseDir+'msd_vac_dicts_keys_antibioticTimeCourse.pkl'

with open(pklFile, 'wb') as f:
    
    pickle.dump([msd_dict,keys_msd_dict,vac_dict,keys_vac_dict],f)
#...

#%% TO RELOAD DICTIONARIES AND LISTS OF KEYS for antibiotic time course
baseDir = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis'

pklFile = baseDir+'msd_vac_dicts_keys_antibioticTimeCourse.pkl'

with open(pklFile, 'rb') as f:  # Python 3: open(..., 'rb')
#    f.open(pklFile,'rb')
    msd_dict, keys_msd_dict, vac_dict, keys_vac_dict = pickle.load(f)
    
    
#%%
txtDictFile = baseDir+'/antibioticTimeCourse/msd_data_all.txt'


with open(txtDictFile, 'w') as f:
    for key, value in msd_dict.items():
        f.write('%s:%s\n' % (key, value[4]))
    
#%%
    
# choose plot colors
col1 = '#0004ff'
col2 = '#ff0000'
col3 = '#00fffb'
col4 = '#e600ff'

# choose plot colors
colCyan = '#00fffb' #col3 = '#00fffb'
colMagenta = '#e600ff' #col4 = '#e600ff'
colRed = '#ff0000'#col1 = '#0004ff' 
colBlue = '#0004ff'#col2 = '#ff0000'
colBlack = '#000000'

#%% for all simple plots (control plus experiment, floA, floT)
    
infoLists = [
        [[0,2],['FloA +LBO2','FloT +LBO2','FloA -LBO2','FloT -LBO2'],[colMagenta,colCyan,colRed,colBlue]],
        [[4,6],['FloA Amp','FloT Amp','FloT Control','FloA Control'],[colRed,colBlue,colCyan,colMagenta]],
        [[8,10],['FloT BNZ','FloA BNZ','FloA Control','FloT Control'],[colBlue,colRed,colMagenta,colCyan]],
        [[16,18],['FloA Fos','FloT Fos','FloT Control','FloA Control'],[colRed,colBlue,colCyan,colMagenta]],
        [[20,22],['FloA Nis','FloT Nis','FloT Control','FloA Control'],[colRed,colBlue,colCyan,colMagenta]],
        [[24,26],['FloT Protoplats','FloA Protoplasts','FloA Control','FloT Control'],[colBlue,colRed,colMagenta,colCyan]],
        [[32,34],['FloA Tun','FloT Tun','FloT Control','FloA Control'],[colRed,colBlue,colCyan,colMagenta]],
        [[38,40],['FloT Val','FloA Val','FloA Control','FloT Control'],[colBlue,colRed,colMagenta,colCyan]],
        [[42,44],['FloA Van','FloT Van','FloT Control','FloA Control'],[colRed,colBlue,colCyan,colMagenta]],
        [[46,48],['FloT ∆dlt','FloA ∆dlt','FloA Control','FloT Control'],[colBlue,colRed,colMagenta,colCyan]],
        [[52,54],['FloT ∆mreB ∆BH ∆mbl','FloA ∆mreB ∆BH ∆mbl','FloT ∆rsgI','FloA ∆rsgI'],[colBlue,colRed,colCyan,colMagenta]],
        [[56,58],['FloA ∆rsgI','FloT ∆rsgI','FloT Control','FloA Control'],[colRed,colBlue,colCyan,colMagenta]],
        [[60,62],['FloA ∆pbpC','FloT ∆pbpC','FloT Control','FloA Control'],[colRed,colBlue,colCyan,colMagenta]],
        [[64,66],['FloT DMSO','FloA DMSO','FloT Control','FloA Control'],[colBlue,colRed,colCyan,colMagenta]],
        [[68,36],['FloT TunWTA','FloA TunWTA','FloA Control','FloTControl'],[colBlue,colRed,colMagenta,colCyan]],
        [[76,78],['FloA ∆flotillin','FloT ∆flotillin','FloA Control','FloT Control'],[colRed,colBlue,colMagenta,colCyan]],
        [[80,82],['FloA ∆tagU','FloT ∆tagU','FloA Control','FloT Control'],[colRed,colBlue,colMagenta,colCyan]],
        [[84,86],['FloT ∆ugtP','FloA ∆ugtP','FloA Control','FloT Control'],[colBlue,colRed,colMagenta,colCyan]],
        ]

# TO CHECK:
# -keys_msd_dict[50] (m delta mreB)
# -[70], spH*.1
# -[72], s delta dlt+Mg
# -[74], s dleta dltA
# -[92], mreb
# -[94], PBP3

doSave = 1
saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'

# legend font size
lfs = 10
# axis font size
afs = 24
# tick font size
tfs = 18
# line width
lw = 2

# choose plot colors
col1 = '#0004ff'
col2 = '#ff0000'
col3 = '#00fffb'
col4 = '#e600ff'

# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict

# pick opacity for shade plot
alphaVal = 0.25

# pick error format
errorFormat = 'ci'
#errorFormat = 'sem'
#errorFormat = 'std'

counter = 0
for infoList in infoLists:
    
    counter += 1
    
    print(infoList)
    print('info list:')
    print(counter)
    print(' ')
    # get indices 
    indsMain = infoList[0]
    
    # get legend entires
    legEntries = infoList[1]
    
    # get plot colors
    col = infoList[2]
    
    # define filename as string of legend entries
    fname = legEntries[0].replace(" ","")+'_'\
    +legEntries[1].replace(" ","")+'_'\
    +legEntries[2].replace(" ","")+'_'\
    +legEntries[3].replace(" ","")

#    returnList = [
#                    dict_vals[subKey0],
#                    dict_cis[subKey0],
#                    dict_cis2[subKey0],
#                    dict_stds[subKey0],
#                    dict_nSamples[subKey0],
#                    subKey0,
#                    ]

    # get lists with plotting info
    returnList_1 = returnListForPlotting(dictToUse,keysToUse,indsMain[0])
    
    # get lists with plotting info
    returnList_2 = returnListForPlotting(dictToUse,keysToUse,indsMain[1])


    # do plotting
    
#    yExt1,x1,y1 = doPlotCondition(returnList_1[0:5],col1,legEntries[0],lw,alphaVal,errorFormat)
#        
#    yExt2,x2,y2 = doPlotCondition(returnList_1[6:11],col2,legEntries[1],lw,alphaVal,errorFormat)
#    
#    yExt3,x3,y3 = doPlotCondition(returnList_2[0:5],col3,legEntries[2],lw,alphaVal,errorFormat)
#        
#    yExt4,x4,y4 = doPlotCondition(returnList_2[6:11],col4,legEntries[3],lw,alphaVal,errorFormat)
    
    yExt1,x1,y1 = doPlotCondition(returnList_1[0:5],col[0],legEntries[0],lw,alphaVal,errorFormat)
        
    yExt2,x2,y2 = doPlotCondition(returnList_1[6:11],col[1],legEntries[1],lw,alphaVal,errorFormat)
    
    yExt3,x3,y3 = doPlotCondition(returnList_2[0:5],col[2],legEntries[2],lw,alphaVal,errorFormat)
        
    yExt4,x4,y4 = doPlotCondition(returnList_2[6:11],col[3],legEntries[3],lw,alphaVal,errorFormat)


    # get min and max y
    yExt = np.concatenate((yExt1,yExt2,yExt3,yExt4),axis=0)
    yMin = np.min(yExt)
    yMax = np.max(yExt)
    
    # get array for x values
    x = x1
    
    
    
    # for reference curves
    xMin = 0
    xMax = -1
    #multFactor = yMax
    
    # adjust start
    multFactorX = 1
    multFactorY1 = .07
    multFactorY2 = .3
    multFactorY3 = .032
    #multFactorY1 = .1
    #multFactorY2 = 1
    xDelta = .1
    
    # adjusts end of ref line
    xMaxOffset1 = 0
    xMaxOffset2 = .8
    xMaxOffset3 = 0
    
    xMinScale = .6
    
    xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
    xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
    xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX
    
    #b = np.min(x)+.1
    #b=10e-2
    #b=yMax
    #b=xtest1[0]-(xtest1[0]-.001)
    b=0
    #xtest2[0] = x[0]
    ytest1 = (xtest1+b)*multFactorY1
    #ytest1 = ytest1-ytest1[0]+y[0]+b
    ytest2 = (np.square(xtest2)+b)*multFactorY2
    #ytest2 = ytest2-ytest2[0]+y[0]+b
    ytest3 = (np.power(xtest3,0.5))*multFactorY3
    
    iSt = 0;
    plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
    plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
    plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)
    
#    plt.legend(loc=4,frameon=False,fontsize=lfs)
    
    plt.xscale('log')
    plt.yscale('log')
    
    plt.xlabel(r'$\tau$ (s)',fontsize=afs)
    plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)
    
    plt.xticks(fontsize=tfs)
    plt.yticks(fontsize=tfs)
    plt.xlim(x[0]-.05,x[-1]+1)
    plt.ylim(.001,.2)
    
    if doSave:
        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'.pdf',bbox_inches="tight")
    #...
        
    plt.legend(loc=4,frameon=False,fontsize=lfs)

    if doSave:
        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'_withLegend.pdf',bbox_inches="tight")
    #...
    
    plt.clf()
#...
    
#%% for antibiotic timecourse data

doSave = 1

# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict

# pick opacity for shade plot
alphaVal = 0.25

# pick error format
errorFormat = 'ci'
#errorFormat = 'sem'
#errorFormat = 'std'

counter = 0
for i in range(0,len(keys_msd_dict),2):
    
    counter += 1
    
    infoKey = keysToUse[i]
    print(infoKey)
    print('key:')
    print(counter)
    print(' ')
    
    # get lists with msd info
    returnList = returnListForPlotting(dictToUse,keysToUse,i)
    
    #    returnList = [
#                    dict_vals[subKey0],
#                    dict_cis[subKey0],
#                    dict_cis2[subKey0],
#                    dict_stds[subKey0],
#                    dict_nSamples[subKey0],
#                    subKey0,
#                    ]


    # get init x values (tau)
    x1 = returnList[0][0,0]
    
    # get init y values (msd)
    y1 = returnList[0][0,1]
    
    # get init D
    D1 = y1/(4*x1)
    
    # get sem
    D1_err = (returnList[3][0]/returnList[4][0])/(4*x1)
    
#    D1_err = (returnList[3][0]/returnList[4][0])/(4*x1)
    
    # get id of track for save string
    ID1 = returnList[5]
    
    print(D1)
    print(D1_err)
#    print(returnList[3][0])
#    print(returnList[4][0])
    print(' ')
    
    
    # get x values (tau)
    x2 = returnList[6][0,0]
    
    # get y values (msd or vac)
    y2 = returnList[6][0,1]
    
    # get init D
    D2 = y2/(4*x2)
    
    # get ID of track for save string
    ID2 = returnList[11]
    
    # get sem
    D2_err = (returnList[9][0]/returnList[10][0])/(4*x1)
    
    print(D2)
    print(D2_err)
#    print(returnList[9][0])
#    print(returnList[10][0])
    print(' ')
    
    
    if doSave:
    
        txtDictFile = baseDir+'/antibioticTimeCourse/'+ID1+'_D.txt'
        
        with open(txtDictFile, 'w') as f:
            f.write('%s:%s\n' % (ID1, D1))
            
            
        txtDictFile = baseDir+'/antibioticTimeCourse/'+ID2+'_D.txt'
        
        with open(txtDictFile, 'w') as f:
            f.write('%s:%s\n' % (ID2, D2))
        
            txtDictFile = baseDir+'/antibioticTimeCourse/'+ID1+'_D_err_sem.txt'
            
            
        
        with open(txtDictFile, 'w') as f:
            f.write('%s:%s\n' % (ID1, D1_err))
            
            
        txtDictFile = baseDir+'/antibioticTimeCourse/'+ID2+'_D_err_sem.txt'
        
        with open(txtDictFile, 'w') as f:
            f.write('%s:%s\n' % (ID2, D2_err))
        
#    # get lists with plotting info
#    returnList_2 = returnListForPlotting(dictToUse,keysToUse,i)



#%% additional plots (see email from Rabea)


#%% just floa, flot controls
infoLists  =  [
        [[98,100],['FloA','FloT']],
        ]

doSave = 0
saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'

# legend font size
lfs = 10
# axis font size
afs = 24
# tick font size
tfs = 18
# line width
lw = 2

# choose plot colors

col2 = '#00fffb'
col1 = '#e600ff'

# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict

# pick opacity for shade plot
alphaVal = 0.25

# pick error format
errorFormat = 'ci'
#errorFormat = 'sem'
#errorFormat = 'std'

counter = 0
for infoList in infoLists:
    
    counter += 1
    
    print(infoList)
    print('info list:')
    print(counter)
    print(' ')
    # get indices 
    indsMain = infoList[0]
    
    # get legend entires
    legEntries = infoList[1]
    
    # define filename as string of legend entries
    fname = legEntries[0].replace(" ","")+'_'\
    +legEntries[1].replace(" ","")


#    returnList = [
#                    dict_vals[subKey0],
#                    dict_cis[subKey0],
#                    dict_cis2[subKey0],
#                    dict_stds[subKey0],
#                    dict_nSamples[subKey0],
#                    subKey0,
#                    ]

    # get lists with plotting info
    returnList_1 = returnListForPlotting(dictToUse,keysToUse,indsMain[0])
    
    # get lists with plotting info
    returnList_2 = returnListForPlotting(dictToUse,keysToUse,indsMain[1])
    


    # do plotting
    
    yExt1,x1,y1 = doPlotCondition(returnList_1[0:5],col1,legEntries[0],lw,alphaVal,errorFormat)
            
    yExt2,x2,y2 = doPlotCondition(returnList_2[0:5],col2,legEntries[1],lw,alphaVal,errorFormat)

    # get min and max y
    yExt = np.concatenate((yExt1,yExt2),axis=0)
    yMin = np.min(yExt)
    yMax = np.max(yExt)
    
    # get array for x values
    x = x1
    
    # for reference curves
    xMin = 0
    xMax = -1
    #multFactor = yMax
    
    # adjust start
    multFactorX = 1
    multFactorY1 = .07
    multFactorY2 = .3
    multFactorY3 = .032
    #multFactorY1 = .1
    #multFactorY2 = 1
    xDelta = .1
    
    # adjusts end of ref line
    xMaxOffset1 = 0
    xMaxOffset2 = .8
    xMaxOffset3 = 0
    
    xMinScale = .6
    
    xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
    xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
    xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX
    
    #b = np.min(x)+.1
    #b=10e-2
    #b=yMax
    #b=xtest1[0]-(xtest1[0]-.001)
    b=0
    #xtest2[0] = x[0]
    ytest1 = (xtest1+b)*multFactorY1
    #ytest1 = ytest1-ytest1[0]+y[0]+b
    ytest2 = (np.square(xtest2)+b)*multFactorY2
    #ytest2 = ytest2-ytest2[0]+y[0]+b
    ytest3 = (np.power(xtest3,0.5))*multFactorY3
    
    iSt = 0;
    plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
    plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
    plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)
    
    
    plt.xscale('log')
    plt.yscale('log')
    
    plt.xlabel(r'$\tau$ (s)',fontsize=afs)
    plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)
    
    plt.xticks(fontsize=tfs)
    plt.yticks(fontsize=tfs)
    plt.xlim(x[0]-.05,x[-1]+1)
    plt.ylim(.001,.2)
    
    
    
    if doSave:
        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'.pdf',bbox_inches="tight")
    #...
    
    plt.legend(loc=4,frameon=False,fontsize=lfs)
    
    if doSave:
        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'_withLegend.pdf',bbox_inches="tight")
    #...
    
#    plt.clf()
    
#%% just floa, flot controls PLAYING AROUND
infoLists  =  [
        [[98,100],['FloA','FloT']],
        ]

doSave = 0
saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'

# legend font size
lfs = 10
# axis font size
afs = 24
# tick font size
tfs = 18
# line width
lw = 2

# choose plot colors

col2 = '#00fffb'
col1 = '#e600ff'

# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict

# pick opacity for shade plot
alphaVal = 0.25

# pick error format
errorFormat = 'ci'
#errorFormat = 'sem'
#errorFormat = 'std'

counter = 0
for infoList in infoLists:
    
    counter += 1
    
    print(infoList)
    print('info list:')
    print(counter)
    print(' ')
    # get indices 
    indsMain = infoList[0]
    
    # get legend entires
    legEntries = infoList[1]
    
    # define filename as string of legend entries
    fname = legEntries[0].replace(" ","")+'_'\
    +legEntries[1].replace(" ","")


#    returnList = [
#                    dict_vals[subKey0],
#                    dict_cis[subKey0],
#                    dict_cis2[subKey0],
#                    dict_stds[subKey0],
#                    dict_nSamples[subKey0],
#                    subKey0,
#                    ]

    # get lists with plotting info
    returnList_1 = returnListForPlotting(dictToUse,keysToUse,indsMain[0])
    
    # get lists with plotting info
    returnList_2 = returnListForPlotting(dictToUse,keysToUse,indsMain[1])
    


    # do plotting
    
    yExt1,x1,y1 = doPlotCondition(returnList_1[0:5],col1,legEntries[0],lw,alphaVal,errorFormat)
            
    yExt2,x2,y2 = doPlotCondition(returnList_2[0:5],col2,legEntries[1],lw,alphaVal,errorFormat)

    # get min and max y
    yExt = np.concatenate((yExt1,yExt2),axis=0)
    yMin = np.min(yExt)
    yMax = np.max(yExt)
    
    # get array for x values
    x = x1
    
    # for reference curves
    xMin = 0
    xMax = -1
    #multFactor = yMax
    
    # adjust start
    multFactorX = 1
    multFactorY1 = .07
    multFactorY2 = .3
    multFactorY3 = .032
    #multFactorY1 = .1
    #multFactorY2 = 1
    xDelta = .1
    
    # adjusts end of ref line
    xMaxOffset1 = 0
    xMaxOffset2 = .8
    xMaxOffset3 = 0
    
    xMinScale = .6
    
    xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
    xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
    xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX
    
    #b = np.min(x)+.1
    #b=10e-2
    #b=yMax
    #b=xtest1[0]-(xtest1[0]-.001)
    b=0
    #xtest2[0] = x[0]
    ytest1 = (xtest1+b)*multFactorY1
    #ytest1 = ytest1-ytest1[0]+y[0]+b
    ytest2 = (np.square(xtest2)+b)*multFactorY2
    #ytest2 = ytest2-ytest2[0]+y[0]+b
    exp = 0.5
    ytest3 = (np.power(xtest3,exp))*multFactorY3
    
    iSt = 0;
    plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
    plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
    plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)
    
    
#    plt.xscale('log')
#    plt.yscale('log')
    
    plt.xlabel(r'$\tau$ (s)',fontsize=afs)
    plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)
    
    plt.xticks(fontsize=tfs)
    plt.yticks(fontsize=tfs)
#    plt.xlim(.001,.2)
    plt.ylim(.001,.05)
    plt.xlim(x[0]-.05,x[-1]+1)
#    plt.ylim(.001,.2)
    
    
    if doSave:
        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'.pdf',bbox_inches="tight")
    #...
    
    plt.legend(loc=4,frameon=False,fontsize=lfs)
    
    if doSave:
        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'_withLegend.pdf',bbox_inches="tight")
    #...
    
#    plt.clf()
    

#%% chimera plots
infoLists  =  [
        [[90,88,100,98],['FloTntAct','FloAntTct','FloT','FloA']],
        ]

doSave = 0
saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'

# legend font size
lfs = 10
# axis font size
afs = 24
# tick font size
tfs = 18
# line width
lw = 2

# choose plot colors
col1 = '#0004ff'
col2 = '#ff0000'
col3 = '#00fffb'
col4 = '#e600ff'

# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict

# pick opacity for shade plot
alphaVal = 0.25

# pick error format
errorFormat = 'ci'
#errorFormat = 'sem'
#errorFormat = 'std'

counter = 0
for infoList in infoLists:
    
    counter += 1
    
    print(infoList)
    print('info list:')
    print(counter)
    print(' ')
    # get indices 
    indsMain = infoList[0]
    
    # get legend entires
    legEntries = infoList[1]
    
    # define filename as string of legend entries
    fname = legEntries[0].replace(" ","")+'_'\
    +legEntries[1].replace(" ","")+'_'\
    +legEntries[2].replace(" ","")+'_'\
    +legEntries[3].replace(" ","")


#    returnList = [
#                    dict_vals[subKey0],
#                    dict_cis[subKey0],
#                    dict_cis2[subKey0],
#                    dict_stds[subKey0],
#                    dict_nSamples[subKey0],
#                    subKey0,
#                    ]

    # get lists with plotting info
    returnList_1 = returnListForPlotting(dictToUse,keysToUse,indsMain[0])
    
    # get lists with plotting info
    returnList_2 = returnListForPlotting(dictToUse,keysToUse,indsMain[1])
    
    # get lists with plotting info
    returnList_3 = returnListForPlotting(dictToUse,keysToUse,indsMain[2])
    
    # get lists with plotting info
    returnList_4 = returnListForPlotting(dictToUse,keysToUse,indsMain[3])


    # do plotting
    
    yExt1,x1,y1 = doPlotCondition(returnList_1[0:5],col1,legEntries[0],lw,alphaVal,errorFormat)
            
    yExt2,x2,y2 = doPlotCondition(returnList_2[0:5],col2,legEntries[1],lw,alphaVal,errorFormat)
        
    yExt3,x3,y3 = doPlotCondition(returnList_3[0:5],col3,legEntries[2],lw,alphaVal,errorFormat)
    
    yExt4,x4,y4 = doPlotCondition(returnList_4[0:5],col4,legEntries[3],lw,alphaVal,errorFormat)

    # get min and max y
    yExt = np.concatenate((yExt1,yExt2),axis=0)
    yMin = np.min(yExt)
    yMax = np.max(yExt)
    
    # get array for x values
    x = x1
    
    # for reference curves
    xMin = 0
    xMax = -1
    #multFactor = yMax
    
    # adjust start
    multFactorX = 1
    multFactorY1 = .07
    multFactorY2 = .3
    multFactorY3 = .032
    #multFactorY1 = .1
    #multFactorY2 = 1
    xDelta = .1
    
    # adjusts end of ref line
    xMaxOffset1 = 0
    xMaxOffset2 = .8
    xMaxOffset3 = 0
    
    xMinScale = .6
    
    xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
    xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
    xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX
    
    #b = np.min(x)+.1
    #b=10e-2
    #b=yMax
    #b=xtest1[0]-(xtest1[0]-.001)
    b=0
    #xtest2[0] = x[0]
    ytest1 = (xtest1+b)*multFactorY1
    #ytest1 = ytest1-ytest1[0]+y[0]+b
    ytest2 = (np.square(xtest2)+b)*multFactorY2
    #ytest2 = ytest2-ytest2[0]+y[0]+b
    ytest3 = (np.power(xtest3,0.5))*multFactorY3
    
    iSt = 0;
    plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
    plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
    plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)
    
    
    plt.xscale('log')
    plt.yscale('log')
    
    plt.xlabel(r'$\tau$ (s)',fontsize=afs)
    plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)
    
    plt.xticks(fontsize=tfs)
    plt.yticks(fontsize=tfs)
    plt.xlim(x[0]-.05,x[-1]+1)
    plt.ylim(.001,.2)
    
    
    
    if doSave:
        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'.pdf',bbox_inches="tight")
    #...
    
    plt.legend(loc=4,frameon=False,fontsize=lfs)

    
    if doSave:
        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'_withLegend.pdf',bbox_inches="tight")
    #...
#    plt.clf()
    
#%% chimera plots
infoLists  =  [
        [[90,88,100,98],['FloTntAct','FloAntTct','FloT','FloA']],
        ]

doSave = 1
saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'

# legend font size
lfs = 10
# axis font size
afs = 24
# tick font size
tfs = 18
# line width
lw = 2

# choose plot colors
col2 = '#0004ff'
col1 = '#ff0000'
col3 = '#00fffb'
col4 = '#e600ff'

# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict

# pick opacity for shade plot
alphaVal = 0.25

# pick error format
errorFormat = 'ci'
#errorFormat = 'sem'
#errorFormat = 'std'

counter = 0
for infoList in infoLists:
    
    counter += 1
    
    print(infoList)
    print('info list:')
    print(counter)
    print(' ')
    # get indices 
    indsMain = infoList[0]
    
    # get legend entires
    legEntries = infoList[1]
    
    # define filename as string of legend entries
    fname = legEntries[0].replace(" ","")+'_'\
    +legEntries[1].replace(" ","")+'_'\
    +legEntries[2].replace(" ","")+'_'\
    +legEntries[3].replace(" ","")


#    returnList = [
#                    dict_vals[subKey0],
#                    dict_cis[subKey0],
#                    dict_cis2[subKey0],
#                    dict_stds[subKey0],
#                    dict_nSamples[subKey0],
#                    subKey0,
#                    ]

    # get lists with plotting info
    returnList_1 = returnListForPlotting(dictToUse,keysToUse,indsMain[0])
    
    # get lists with plotting info
    returnList_2 = returnListForPlotting(dictToUse,keysToUse,indsMain[1])
    
    # get lists with plotting info
    returnList_3 = returnListForPlotting(dictToUse,keysToUse,indsMain[2])
    
    # get lists with plotting info
    returnList_4 = returnListForPlotting(dictToUse,keysToUse,indsMain[3])


    # do plotting
    
    yExt1,x1,y1 = doPlotCondition(returnList_1[0:5],col1,legEntries[0],lw,alphaVal,errorFormat)
            
    yExt2,x2,y2 = doPlotCondition(returnList_2[0:5],col2,legEntries[1],lw,alphaVal,errorFormat)
        
    yExt3,x3,y3 = doPlotCondition(returnList_3[0:5],col3,legEntries[2],lw,alphaVal,errorFormat)
    
    yExt4,x4,y4 = doPlotCondition(returnList_4[0:5],col4,legEntries[3],lw,alphaVal,errorFormat)

    # get min and max y
    yExt = np.concatenate((yExt1,yExt2),axis=0)
    yMin = np.min(yExt)
    yMax = np.max(yExt)
    
    # get array for x values
    x = x1
    
    # for reference curves
    xMin = 0
    xMax = -1
    #multFactor = yMax
    
    # adjust start
    multFactorX = 1
    multFactorY1 = .07
    multFactorY2 = .3
    multFactorY3 = .032
    #multFactorY1 = .1
    #multFactorY2 = 1
    xDelta = .1
    
    # adjusts end of ref line
    xMaxOffset1 = 0
    xMaxOffset2 = .8
    xMaxOffset3 = 0
    
    xMinScale = .6
    
    xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
    xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
    xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX
    
    #b = np.min(x)+.1
    #b=10e-2
    #b=yMax
    #b=xtest1[0]-(xtest1[0]-.001)
    b=0
    #xtest2[0] = x[0]
    ytest1 = (xtest1+b)*multFactorY1
    #ytest1 = ytest1-ytest1[0]+y[0]+b
    ytest2 = (np.square(xtest2)+b)*multFactorY2
    #ytest2 = ytest2-ytest2[0]+y[0]+b
    ytest3 = (np.power(xtest3,0.5))*multFactorY3
    
    iSt = 0;
    plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
    plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
    plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)
    
    
    plt.xscale('log')
    plt.yscale('log')
    
    plt.xlabel(r'$\tau$ (s)',fontsize=afs)
    plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)
    
    plt.xticks(fontsize=tfs)
    plt.yticks(fontsize=tfs)
    plt.xlim(x[0]-.05,x[-1]+1)
    plt.ylim(.001,.2)
    
    
    
    if doSave:
        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'_colSwapped.pdf',bbox_inches="tight")
    #...
    
    plt.legend(loc=4,frameon=False,fontsize=lfs)

    
    if doSave:
        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'_colSwapped_withLegend.pdf',bbox_inches="tight")
    #...
#    plt.clf()

#%%

#
#
#
#
#%% pbp3 original mtirf35 DEPRECEATED
#infoLists  =  [
#        [[102,30],['PBP3','FloA','FloT']],
#        ]
#
#doSave = 0
#saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'
#
## legend font size
#lfs = 10
## axis font size
#afs = 24
## tick font size
#tfs = 18
## line width
#lw = 2
#
## choose plot colors
#col3 = '#00fffb'
#col2 = '#e600ff'
#col1 = '#000000'
#
## choose the dictionary and associated list of keys(msd or vac)
#dictToUse = msd_dict
#keysToUse = keys_msd_dict
#
## pick opacity for shade plot
#alphaVal = 0.25
#
## pick error format
#errorFormat = 'ci'
##errorFormat = 'sem'
##errorFormat = 'std'
#
#counter = 0
#for infoList in infoLists:
#    
#    counter += 1
#    
#    print(infoList)
#    print('info list:')
#    print(counter)
#    print(' ')
#    # get indices 
#    indsMain = infoList[0]
#    
#    # get legend entires
#    legEntries = infoList[1]
#    
#    # define filename as string of legend entries
#    fname = legEntries[0].replace(" ","")+'_'\
#    +legEntries[1].replace(" ","")+'_'\
#    +legEntries[2].replace(" ","")
#
#
##    returnList = [
##                    dict_vals[subKey0],
##                    dict_cis[subKey0],
##                    dict_cis2[subKey0],
##                    dict_stds[subKey0],
##                    dict_nSamples[subKey0],
##                    subKey0,
##                    ]
#
#    # get lists with plotting info
#    returnList_1 = returnListForPlotting(dictToUse,keysToUse,indsMain[0])
#    
#    # get lists with plotting info
#    returnList_2 = returnListForPlotting(dictToUse,keysToUse,indsMain[1])
#
#    # do plotting
#    
#    yExt1,x1,y1 = doPlotCondition(returnList_1[0:5],col1,legEntries[0],lw,alphaVal,errorFormat)
#            
#    yExt2,x2,y2 = doPlotCondition(returnList_2[0:5],col2,legEntries[1],lw,alphaVal,errorFormat)
#        
#    yExt3,x3,y3 = doPlotCondition(returnList_2[6:11],col3,legEntries[2],lw,alphaVal,errorFormat)
#    
#
#    # get min and max y
#    yExt = np.concatenate((yExt1,yExt2),axis=0)
#    yMin = np.min(yExt)
#    yMax = np.max(yExt)
#    
#    # get array for x values
#    x = x1
#    
#    # for reference curves
#    xMin = 0
#    xMax = -1
#    #multFactor = yMax
#    
#    # adjust start
#    multFactorX = 1
#    multFactorY1 = .07
#    multFactorY2 = .3
#    multFactorY3 = .032
#    #multFactorY1 = .1
#    #multFactorY2 = 1
#    xDelta = .1
#    
#    # adjusts end of ref line
#    xMaxOffset1 = 0
#    xMaxOffset2 = .8
#    xMaxOffset3 = 0
#    
#    xMinScale = .6
#    
#    xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
#    xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
#    xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX
#    
#    #b = np.min(x)+.1
#    #b=10e-2
#    #b=yMax
#    #b=xtest1[0]-(xtest1[0]-.001)
#    b=0
#    #xtest2[0] = x[0]
#    ytest1 = (xtest1+b)*multFactorY1
#    #ytest1 = ytest1-ytest1[0]+y[0]+b
#    ytest2 = (np.square(xtest2)+b)*multFactorY2
#    #ytest2 = ytest2-ytest2[0]+y[0]+b
#    ytest3 = (np.power(xtest3,0.5))*multFactorY3
#    
#    iSt = 0;
#    plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
#    plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
#    plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)
#        
#    plt.xscale('log')
#    plt.yscale('log')
#    
#    plt.xlabel(r'$\tau$ (s)',fontsize=afs)
#    plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)
#    
#    plt.xticks(fontsize=tfs)
#    plt.yticks(fontsize=tfs)
#    plt.xlim(x[0]-.05,x[-1]+1)
#    plt.ylim(.001,.2)
#    
#    
#    
#    if doSave:
#        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'_orig.pdf',bbox_inches="tight")
#    #...
#    
#    plt.legend(loc=4,frameon=False,fontsize=lfs)
#
#    if doSave:
#        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'_orig_withLegend.pdf',bbox_inches="tight")
#    #...
#    
##    plt.clf()
    
    
    
#%% pbp3 mtirf35
infoLists  =  [
        [[94,30],['PBP3','FloA','FloT']],
        ]

doSave = 0
saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'

# legend font size
lfs = 10
# axis font size
afs = 24
# tick font size
tfs = 18
# line width
lw = 2

# choose plot colors
col3 = '#00fffb'
col2 = '#e600ff'
col1 = '#000000'

# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict

# pick opacity for shade plot
alphaVal = 0.25

# pick error format
errorFormat = 'ci'
#errorFormat = 'sem'
#errorFormat = 'std'

counter = 0
for infoList in infoLists:
    
    counter += 1
    
    print(infoList)
    print('info list:')
    print(counter)
    print(' ')
    # get indices 
    indsMain = infoList[0]
    
    # get legend entires
    legEntries = infoList[1]
    
    # define filename as string of legend entries
    fname = legEntries[0].replace(" ","")+'_'\
    +legEntries[1].replace(" ","")+'_'\
    +legEntries[2].replace(" ","")


#    returnList = [
#                    dict_vals[subKey0],
#                    dict_cis[subKey0],
#                    dict_cis2[subKey0],
#                    dict_stds[subKey0],
#                    dict_nSamples[subKey0],
#                    subKey0,
#                    ]

    # get lists with plotting info
    returnList_1 = returnListForPlotting(dictToUse,keysToUse,indsMain[0])
    
    # get lists with plotting info
    returnList_2 = returnListForPlotting(dictToUse,keysToUse,indsMain[1])

    # do plotting
    
    yExt1,x1,y1 = doPlotCondition(returnList_1[47-5:47],col1,legEntries[0],lw,alphaVal,errorFormat)
            
    yExt2,x2,y2 = doPlotCondition(returnList_2[0:5],col2,legEntries[1],lw,alphaVal,errorFormat)
        
    yExt3,x3,y3 = doPlotCondition(returnList_2[6:11],col3,legEntries[2],lw,alphaVal,errorFormat)
    

    # get min and max y
    yExt = np.concatenate((yExt1,yExt2),axis=0)
    yMin = np.min(yExt)
    yMax = np.max(yExt)
    
    # get array for x values
    x = x1
    
    # for reference curves
    xMin = 0
    xMax = -1
    #multFactor = yMax
    
    # adjust start
    multFactorX = 1
    multFactorY1 = .07
    multFactorY2 = .3
    multFactorY3 = .032
    #multFactorY1 = .1
    #multFactorY2 = 1
    xDelta = .1
    
    # adjusts end of ref line
    xMaxOffset1 = 0
    xMaxOffset2 = .8
    xMaxOffset3 = 0
    
    xMinScale = .6
    
    xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
    xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
    xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX
    
    #b = np.min(x)+.1
    #b=10e-2
    #b=yMax
    #b=xtest1[0]-(xtest1[0]-.001)
    b=0
    #xtest2[0] = x[0]
    ytest1 = (xtest1+b)*multFactorY1
    #ytest1 = ytest1-ytest1[0]+y[0]+b
    ytest2 = (np.square(xtest2)+b)*multFactorY2
    #ytest2 = ytest2-ytest2[0]+y[0]+b
    ytest3 = (np.power(xtest3,0.5))*multFactorY3
    
    iSt = 0;
    plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
    plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
    plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)
        
    plt.xscale('log')
    plt.yscale('log')
    
    plt.xlabel(r'$\tau$ (s)',fontsize=afs)
    plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)
    
    plt.xticks(fontsize=tfs)
    plt.yticks(fontsize=tfs)
    plt.xlim(x[0]-.05,x[-1]+1)
    plt.ylim(.001,.2)
    
    
    
    if doSave:
        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'.pdf',bbox_inches="tight")
    #...
    
    plt.legend(loc=4,frameon=False,fontsize=lfs)

    if doSave:
        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'_withLegend.pdf',bbox_inches="tight")
    #...
    
#    plt.clf()
    
    
#%% pbp3 delta floa delta flot
infoLists  =  [
        [[94],['PBP3 Ctrl','∆floA','∆floT']],
        ]

doSave = 0
saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'

# legend font size
lfs = 10
# axis font size
afs = 24
# tick font size
tfs = 18
# line width
lw = 2

# choose plot colors
col3 = '#00fffb'
col2 = '#e600ff'
col1 = '#000000'

# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict

# pick opacity for shade plot
alphaVal = 0.25

# pick error format
errorFormat = 'ci'
#errorFormat = 'sem'
#errorFormat = 'std'

counter = 0
for infoList in infoLists:
    
    counter += 1
    
    print(infoList)
    print('info list:')
    print(counter)
    print(' ')
    # get indices 
    indsMain = infoList[0]
    
    # get legend entires
    legEntries = infoList[1]
    
    # define filename as string of legend entries
    fname = legEntries[0].replace(" ","")+'_'\
    +legEntries[1].replace(" ","")+'_'\
    +legEntries[2].replace(" ","")


#    returnList = [
#                    dict_vals[subKey0],
#                    dict_cis[subKey0],
#                    dict_cis2[subKey0],
#                    dict_stds[subKey0],
#                    dict_nSamples[subKey0],
#                    subKey0,
#                    ]

    # get lists with plotting info
    returnList_1 = returnListForPlotting(dictToUse,keysToUse,indsMain[0])
    


    # do plotting
    
    yExt1,x1,y1 = doPlotCondition(returnList_1[47-5:47],col1,legEntries[0],lw,alphaVal,errorFormat)
            
    yExt2,x2,y2 = doPlotCondition(returnList_1[35-5:35],col2,legEntries[1],lw,alphaVal,errorFormat)
        
    yExt3,x3,y3 = doPlotCondition(returnList_1[41-5:41],col3,legEntries[2],lw,alphaVal,errorFormat)
    

    # get min and max y
    yExt = np.concatenate((yExt1,yExt2),axis=0)
    yMin = np.min(yExt)
    yMax = np.max(yExt)
    
    # get array for x values
    x = x1
    
    # for reference curves
    xMin = 0
    xMax = -1
    #multFactor = yMax
    
    # adjust start
    multFactorX = 1
    multFactorY1 = .07
    multFactorY2 = .3
    multFactorY3 = .032
    #multFactorY1 = .1
    #multFactorY2 = 1
    xDelta = .1
    
    # adjusts end of ref line
    xMaxOffset1 = 0
    xMaxOffset2 = .8
    xMaxOffset3 = 0
    
    xMinScale = .6
    
    xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
    xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
    xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX
    
    #b = np.min(x)+.1
    #b=10e-2
    #b=yMax
    #b=xtest1[0]-(xtest1[0]-.001)
    b=0
    #xtest2[0] = x[0]
    ytest1 = (xtest1+b)*multFactorY1
    #ytest1 = ytest1-ytest1[0]+y[0]+b
    ytest2 = (np.square(xtest2)+b)*multFactorY2
    #ytest2 = ytest2-ytest2[0]+y[0]+b
    ytest3 = (np.power(xtest3,0.5))*multFactorY3
    
    iSt = 0;
    plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
    plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
    plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)
    
    
    plt.xscale('log')
    plt.yscale('log')
    
    plt.xlabel(r'$\tau$ (s)',fontsize=afs)
    plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)
    
    plt.xticks(fontsize=tfs)
    plt.yticks(fontsize=tfs)
    plt.xlim(x[0]-.05,x[-1]+1)
    plt.ylim(.001,.2)
    
    
    
    if doSave:
        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'.pdf',bbox_inches="tight")
    #...
    
    plt.legend(loc=4,frameon=False,fontsize=lfs)

    if doSave:
        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'_withLegend.pdf',bbox_inches="tight")
    #...
#    plt.clf()
    
  
    
#%% pbp3 antibiotic plots
infoLists  =  [
        [[94],['PBP3 Ctrl','Amp','Fos','Tun','Van']],
        ]

doSave = 0
saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'

# legend font size
lfs = 10
# axis font size
afs = 24
# tick font size
tfs = 18
# line width
lw = 2

# choose plot colors
col4 = '#0004ff'
col5 = '#ff0000'
col3 = '#00fffb'
col2 = '#e600ff'
col1 = '#000000'

# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict

# pick opacity for shade plot
alphaVal = 0.25

# pick error format
errorFormat = 'ci'
#errorFormat = 'sem'
#errorFormat = 'std'

counter = 0
for infoList in infoLists:
    
    counter += 1
    
    print(infoList)
    print('info list:')
    print(counter)
    print(' ')
    # get indices 
    indsMain = infoList[0]
    
    # get legend entires
    legEntries = infoList[1]
    
    # define filename as string of legend entries
    fname = legEntries[0].replace(" ","")+'_'\
    +legEntries[1].replace(" ","")+'_'\
    +legEntries[2].replace(" ","")+'_'\
    +legEntries[3].replace(" ","")+'_'\
    +legEntries[4].replace(" ","")
    
#    returnList = [
#                    dict_vals[subKey0],
#                    dict_cis[subKey0],
#                    dict_cis2[subKey0],
#                    dict_stds[subKey0],
#                    dict_nSamples[subKey0],
#                    subKey0,
#                    ]

    # get lists with plotting info
    returnList_1 = returnListForPlotting(dictToUse,keysToUse,indsMain[0])
    


    # do plotting
    #['PBP3 Ctrl','Amp','Fos','Tun','Van']
    yExt1,x1,y1 = doPlotCondition(returnList_1[47-5:47],col1,legEntries[0],lw,alphaVal,errorFormat)
            
    yExt2,x2,y2 = doPlotCondition(returnList_1[6:11],col2,legEntries[1],lw,alphaVal,errorFormat)
        
    yExt3,x3,y3 = doPlotCondition(returnList_1[24:29],col3,legEntries[2],lw,alphaVal,errorFormat)
    
    yExt4,x4,y4 = doPlotCondition(returnList_1[12:17],col4,legEntries[3],lw,alphaVal,errorFormat)
    
    yExt5,x5,y5 = doPlotCondition(returnList_1[0:5],col5,legEntries[4],lw,alphaVal,errorFormat)

    # get min and max y
    yExt = np.concatenate((yExt1,yExt2),axis=0)
    yMin = np.min(yExt)
    yMax = np.max(yExt)
    
    # get array for x values
    x = x1
    
    # for reference curves
    xMin = 0
    xMax = -1
    #multFactor = yMax
    
    # adjust start
    multFactorX = 1
    multFactorY1 = .07
    multFactorY2 = .3
    multFactorY3 = .032
    #multFactorY1 = .1
    #multFactorY2 = 1
    xDelta = .1
    
    # adjusts end of ref line
    xMaxOffset1 = 0
    xMaxOffset2 = .8
    xMaxOffset3 = 0
    
    xMinScale = .6
    
    xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
    xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
    xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX
    
    #b = np.min(x)+.1
    #b=10e-2
    #b=yMax
    #b=xtest1[0]-(xtest1[0]-.001)
    b=0
    #xtest2[0] = x[0]
    ytest1 = (xtest1+b)*multFactorY1
    #ytest1 = ytest1-ytest1[0]+y[0]+b
    ytest2 = (np.square(xtest2)+b)*multFactorY2
    #ytest2 = ytest2-ytest2[0]+y[0]+b
    ytest3 = (np.power(xtest3,0.5))*multFactorY3
    
    iSt = 0;
    plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
    plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
    plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)
    
    
    plt.xscale('log')
    plt.yscale('log')
    
    plt.xlabel(r'$\tau$ (s)',fontsize=afs)
    plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)
    
    plt.xticks(fontsize=tfs)
    plt.yticks(fontsize=tfs)
    plt.xlim(x[0]-.05,x[-1]+1)
    plt.ylim(.001,.2)
    
    
    
    if doSave:
        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'.pdf',bbox_inches="tight")
    #...
    
    plt.legend(loc=4,frameon=False,fontsize=lfs)

    if doSave:
        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'_withLegend.pdf',bbox_inches="tight")
    #...
#    plt.clf()
    
    
#%% pbp3 cells protoplasts
infoLists  =  [
        [[94],['PBP3 cells','PBP3 protoplasts']],
        ]

doSave = 0
saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'

# legend font size
lfs = 10
# axis font size
afs = 24
# tick font size
tfs = 18
# line width
lw = 2

# choose plot colors
col2 = '#0004ff'
col1 = '#ff0000'


# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict

# pick opacity for shade plot
alphaVal = 0.25

# pick error format
errorFormat = 'ci'
#errorFormat = 'sem'
#errorFormat = 'std'

counter = 0
for infoList in infoLists:
    
    counter += 1
    
    print(infoList)
    print('info list:')
    print(counter)
    print(' ')
    # get indices 
    indsMain = infoList[0]
    
    # get legend entires
    legEntries = infoList[1]
    
    # define filename as string of legend entries
    fname = legEntries[0].replace(" ","")+'_'\
    +legEntries[1].replace(" ","")
    
#    returnList = [
#                    dict_vals[subKey0],
#                    dict_cis[subKey0],
#                    dict_cis2[subKey0],
#                    dict_stds[subKey0],
#                    dict_nSamples[subKey0],
#                    subKey0,
#                    ]

    # get lists with plotting info
    returnList_1 = returnListForPlotting(dictToUse,keysToUse,indsMain[0])
    


    # do plotting
    
    yExt1,x1,y1 = doPlotCondition(returnList_1[48:53],col1,legEntries[0],lw,alphaVal,errorFormat)
            
    yExt2,x2,y2 = doPlotCondition(returnList_1[18:23],col2,legEntries[1],lw,alphaVal,errorFormat)

    # get min and max y
    yExt = np.concatenate((yExt1,yExt2),axis=0)
    yMin = np.min(yExt)
    yMax = np.max(yExt)
    
    # get array for x values
    x = x1
    
    # for reference curves
    xMin = 0
    xMax = -1
    #multFactor = yMax
    
    # adjust start
    multFactorX = 1
    multFactorY1 = .07
    multFactorY2 = .3
    multFactorY3 = .032
    #multFactorY1 = .1
    #multFactorY2 = 1
    xDelta = .1
    
    # adjusts end of ref line
    xMaxOffset1 = 0
    xMaxOffset2 = .8
    xMaxOffset3 = 0
    
    xMinScale = .6
    
    xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
    xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
    xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX
    
    #b = np.min(x)+.1
    #b=10e-2
    #b=yMax
    #b=xtest1[0]-(xtest1[0]-.001)
    b=0
    #xtest2[0] = x[0]
    ytest1 = (xtest1+b)*multFactorY1
    #ytest1 = ytest1-ytest1[0]+y[0]+b
    ytest2 = (np.square(xtest2)+b)*multFactorY2
    #ytest2 = ytest2-ytest2[0]+y[0]+b
    ytest3 = (np.power(xtest3,0.5))*multFactorY3
    
    iSt = 0;
    plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
    plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
    plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)
        
    plt.xscale('log')
    plt.yscale('log')
    
    plt.xlabel(r'$\tau$ (s)',fontsize=afs)
    plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)
    
    plt.xticks(fontsize=tfs)
    plt.yticks(fontsize=tfs)
    plt.xlim(x[0]-.05,x[-1]+1)
    plt.ylim(.001,.2)
    
    
    
    if doSave:
        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'.pdf',bbox_inches="tight")
    #...
    
    plt.legend(loc=4,frameon=False,fontsize=lfs)

    if doSave:
        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'_withLegend.pdf',bbox_inches="tight")
        
    
#    plt.clf()
#%%

#
#
#
#
#%% dltd delta floa delta flot
infoLists  =  [
        [[92],['DltD Ctrl','∆floA','∆floT']],
        ]

doSave = 0
saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'

# legend font size
lfs = 10
# axis font size
afs = 24
# tick font size
tfs = 18
# line width
lw = 2

# choose plot colors
col3 = '#00fffb'
col2 = '#e600ff'
col1 = '#000000'

# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict

# pick opacity for shade plot
alphaVal = 0.25

# pick error format
errorFormat = 'ci'
#errorFormat = 'sem'
#errorFormat = 'std'

counter = 0
for infoList in infoLists:
    
    counter += 1
    
    print(infoList)
    print('info list:')
    print(counter)
    print(' ')
    # get indices 
    indsMain = infoList[0]
    
    # get legend entires
    legEntries = infoList[1]
    
    # define filename as string of legend entries
    fname = legEntries[0].replace(" ","")+'_'\
    +legEntries[1].replace(" ","")+'_'\
    +legEntries[2].replace(" ","")


#    returnList = [
#                    dict_vals[subKey0],
#                    dict_cis[subKey0],
#                    dict_cis2[subKey0],
#                    dict_stds[subKey0],
#                    dict_nSamples[subKey0],
#                    subKey0,
#                    ]

    # get lists with plotting info
    returnList_1 = returnListForPlotting(dictToUse,keysToUse,indsMain[0])
    


    # do plotting
#    ['DltD Ctrl','∆floA','∆floT']
    yExt1,x1,y1 = doPlotCondition(returnList_1[47-5:47],col1,legEntries[0],lw,alphaVal,errorFormat)
            
    yExt2,x2,y2 = doPlotCondition(returnList_1[35-5:35],col2,legEntries[1],lw,alphaVal,errorFormat)
        
    yExt3,x3,y3 = doPlotCondition(returnList_1[41-5:41],col3,legEntries[2],lw,alphaVal,errorFormat)
    

    # get min and max y
    yExt = np.concatenate((yExt1,yExt2),axis=0)
    yMin = np.min(yExt)
    yMax = np.max(yExt)
    
    # get array for x values
    x = x1
    
    # for reference curves
    xMin = 0
    xMax = -1
    #multFactor = yMax
    
    # adjust start
    multFactorX = 1
    multFactorY1 = .07
    multFactorY2 = .3
    multFactorY3 = .032
    #multFactorY1 = .1
    #multFactorY2 = 1
    xDelta = .1
    
    # adjusts end of ref line
    xMaxOffset1 = 0
    xMaxOffset2 = .8
    xMaxOffset3 = 0
    
    xMinScale = .6
    
    xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
    xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
    xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX
    
    #b = np.min(x)+.1
    #b=10e-2
    #b=yMax
    #b=xtest1[0]-(xtest1[0]-.001)
    b=0
    #xtest2[0] = x[0]
    ytest1 = (xtest1+b)*multFactorY1
    #ytest1 = ytest1-ytest1[0]+y[0]+b
    ytest2 = (np.square(xtest2)+b)*multFactorY2
    #ytest2 = ytest2-ytest2[0]+y[0]+b
    ytest3 = (np.power(xtest3,0.5))*multFactorY3
    
    iSt = 0;
    plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
    plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
    plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)
    
    
    plt.xscale('log')
    plt.yscale('log')
    
    plt.xlabel(r'$\tau$ (s)',fontsize=afs)
    plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)
    
    plt.xticks(fontsize=tfs)
    plt.yticks(fontsize=tfs)
    plt.xlim(x[0]-.05,x[-1]+1)
    plt.ylim(.001,.2)
    
    
    
    if doSave:
        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'.pdf',bbox_inches="tight")
    #...
    
    plt.legend(loc=4,frameon=False,fontsize=lfs)

    
    if doSave:
        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'_withLegend.pdf',bbox_inches="tight")
    #...
#    plt.clf()
    
  
    
#%% dltd antibiotic plots
    
infoLists  =  [
        [[92],['DltD Ctrl','Amp','Fos','Tun','Van']],
        ]

doSave = 0
saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'

# legend font size
lfs = 10
# axis font size
afs = 24
# tick font size
tfs = 18
# line width
lw = 2

# choose plot colors
col4 = '#0004ff'
col5 = '#ff0000'
col3 = '#00fffb'
col2 = '#e600ff'
col1 = '#000000'

# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict

# pick opacity for shade plot
alphaVal = 0.25

# pick error format
errorFormat = 'ci'
#errorFormat = 'sem'
#errorFormat = 'std'

counter = 0
for infoList in infoLists:
    
    counter += 1
    
    print(infoList)
    print('info list:')
    print(counter)
    print(' ')
    # get indices 
    indsMain = infoList[0]
    
    # get legend entires
    legEntries = infoList[1]
    
    # define filename as string of legend entries
    fname = legEntries[0].replace(" ","")+'_'\
    +legEntries[1].replace(" ","")+'_'\
    +legEntries[2].replace(" ","")+'_'\
    +legEntries[3].replace(" ","")+'_'\
    +legEntries[4].replace(" ","")
    
#    returnList = [
#                    dict_vals[subKey0],
#                    dict_cis[subKey0],
#                    dict_cis2[subKey0],
#                    dict_stds[subKey0],
#                    dict_nSamples[subKey0],
#                    subKey0,
#                    ]

    # get lists with plotting info
    returnList_1 = returnListForPlotting(dictToUse,keysToUse,indsMain[0])
    


    # do plotting
#    [[92],['DltD Ctrl','Amp','Fos','Tun','Van']
    yExt1,x1,y1 = doPlotCondition(returnList_1[47-5:47],col1,legEntries[0],lw,alphaVal,errorFormat)
            
    yExt2,x2,y2 = doPlotCondition(returnList_1[6:11],col2,legEntries[1],lw,alphaVal,errorFormat)
        
    yExt3,x3,y3 = doPlotCondition(returnList_1[24:29],col3,legEntries[2],lw,alphaVal,errorFormat)
    
    yExt4,x4,y4 = doPlotCondition(returnList_1[12:17],col4,legEntries[3],lw,alphaVal,errorFormat)
    
    yExt5,x5,y5 = doPlotCondition(returnList_1[0:5],col5,legEntries[4],lw,alphaVal,errorFormat)

    # get min and max y
    yExt = np.concatenate((yExt1,yExt2),axis=0)
    yMin = np.min(yExt)
    yMax = np.max(yExt)
    
    # get array for x values
    x = x1
    
    # for reference curves
    xMin = 0
    xMax = -1
    #multFactor = yMax
    
    # adjust start
    multFactorX = 1
    multFactorY1 = .07
    multFactorY2 = .3
    multFactorY3 = .032
    #multFactorY1 = .1
    #multFactorY2 = 1
    xDelta = .1
    
    # adjusts end of ref line
    xMaxOffset1 = 0
    xMaxOffset2 = .8
    xMaxOffset3 = 0
    
    xMinScale = .6
    
    xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
    xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
    xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX
    
    #b = np.min(x)+.1
    #b=10e-2
    #b=yMax
    #b=xtest1[0]-(xtest1[0]-.001)
    b=0
    #xtest2[0] = x[0]
    ytest1 = (xtest1+b)*multFactorY1
    #ytest1 = ytest1-ytest1[0]+y[0]+b
    ytest2 = (np.square(xtest2)+b)*multFactorY2
    #ytest2 = ytest2-ytest2[0]+y[0]+b
    ytest3 = (np.power(xtest3,0.5))*multFactorY3
    
    iSt = 0;
    plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
    plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
    plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)
    
    
    plt.xscale('log')
    plt.yscale('log')
    
    plt.xlabel(r'$\tau$ (s)',fontsize=afs)
    plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)
    
    plt.xticks(fontsize=tfs)
    plt.yticks(fontsize=tfs)
    plt.xlim(x[0]-.05,x[-1]+1)
    plt.ylim(.001,.2)
    
    
    
    if doSave:
        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'.pdf',bbox_inches="tight")
    #...
    
    plt.legend(loc=4,frameon=False,fontsize=lfs)
    
    if doSave:
        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'_withLegend.pdf',bbox_inches="tight")
    #...

#    plt.clf()
    
    
#%% dltd cells protoplasts
infoLists  =  [
        [[92],['DltD cells','DltD protoplasts']],
        ]

doSave = 0
saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'

# legend font size
lfs = 10
# axis font size
afs = 24
# tick font size
tfs = 18
# line width
lw = 2

# choose plot colors
col2 = '#0004ff'
col1 = '#ff0000'


# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict

# pick opacity for shade plot
alphaVal = 0.25

# pick error format
errorFormat = 'ci'
#errorFormat = 'sem'
#errorFormat = 'std'

counter = 0
for infoList in infoLists:
    
    counter += 1
    
    print(infoList)
    print('info list:')
    print(counter)
    print(' ')
    # get indices 
    indsMain = infoList[0]
    
    # get legend entires
    legEntries = infoList[1]
    
    # define filename as string of legend entries
    fname = legEntries[0].replace(" ","")+'_'\
    +legEntries[1].replace(" ","")
    
#    returnList = [
#                    dict_vals[subKey0],
#                    dict_cis[subKey0],
#                    dict_cis2[subKey0],
#                    dict_stds[subKey0],
#                    dict_nSamples[subKey0],
#                    subKey0,
#                    ]

    # get lists with plotting info
    returnList_1 = returnListForPlotting(dictToUse,keysToUse,indsMain[0])
    


    # do plotting
    
    yExt1,x1,y1 = doPlotCondition(returnList_1[48:53],col1,legEntries[0],lw,alphaVal,errorFormat)
            
    yExt2,x2,y2 = doPlotCondition(returnList_1[18:23],col2,legEntries[1],lw,alphaVal,errorFormat)

    # get min and max y
    yExt = np.concatenate((yExt1,yExt2),axis=0)
    yMin = np.min(yExt)
    yMax = np.max(yExt)
    
    # get array for x values
    x = x1
    
    # for reference curves
    xMin = 0
    xMax = -1
    #multFactor = yMax
    
    # adjust start
    multFactorX = 1
    multFactorY1 = .07
    multFactorY2 = .3
    multFactorY3 = .032
    #multFactorY1 = .1
    #multFactorY2 = 1
    xDelta = .1
    
    # adjusts end of ref line
    xMaxOffset1 = 0
    xMaxOffset2 = .8
    xMaxOffset3 = 0
    
    xMinScale = .6
    
    xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
    xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
    xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX
    
    #b = np.min(x)+.1
    #b=10e-2
    #b=yMax
    #b=xtest1[0]-(xtest1[0]-.001)
    b=0
    #xtest2[0] = x[0]
    ytest1 = (xtest1+b)*multFactorY1
    #ytest1 = ytest1-ytest1[0]+y[0]+b
    ytest2 = (np.square(xtest2)+b)*multFactorY2
    #ytest2 = ytest2-ytest2[0]+y[0]+b
    ytest3 = (np.power(xtest3,0.5))*multFactorY3
    
    iSt = 0;
    plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
    plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
    plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)
    
    
    plt.xscale('log')
    plt.yscale('log')
    
    plt.xlabel(r'$\tau$ (s)',fontsize=afs)
    plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)
    
    plt.xticks(fontsize=tfs)
    plt.yticks(fontsize=tfs)
    plt.xlim(x[0]-.05,x[-1]+1)
    plt.ylim(.001,.2)
    
    
    
    if doSave:
        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'.pdf',bbox_inches="tight")
    #...
    
    plt.legend(loc=4,frameon=False,fontsize=lfs)
    
    if doSave:
        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'_withLegend.pdf',bbox_inches="tight")
    #...

#    plt.clf()
    
    
#%% FloT plots
infoLists  =  [
        [[48,46,72,74,70],['FloT Ctrl','∆dlt','∆dlt + Mg','∆dlt','pH8.1']],
        ]

doSave = 0
saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'

# legend font size
lfs = 10
# axis font size
afs = 24
# tick font size
tfs = 18
# line width
lw = 2

# choose plot colors
col1 = '#000000'
col2 = '#0004ff'
col3 = '#ff0000'
col4 = '#00fffb'
col5 = '#e600ff'

# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict

# pick opacity for shade plot
alphaVal = 0.25

# pick error format
errorFormat = 'ci'
#errorFormat = 'sem'
#errorFormat = 'std'

counter = 0
for infoList in infoLists:
    
    counter += 1
    
    print(infoList)
    print('info list:')
    print(counter)
    print(' ')
    # get indices 
    indsMain = infoList[0]
    
    # get legend entires
    legEntries = infoList[1]
    
    # define filename as string of legend entries
    fname = legEntries[0].replace(" ","")+'_'\
    +legEntries[1].replace(" ","")+'_'\
    +legEntries[2].replace(" ","")+'_'\
    +legEntries[3].replace(" ","")+'_'\
    +legEntries[4].replace(" ","")

#    returnList = [
#                    dict_vals[subKey0],
#                    dict_cis[subKey0],
#                    dict_cis2[subKey0],
#                    dict_stds[subKey0],
#                    dict_nSamples[subKey0],
#                    subKey0,
#                    ]

    # get lists with plotting info
    returnList_1 = returnListForPlotting(dictToUse,keysToUse,indsMain[0])
    
    # get lists with plotting info
    returnList_2 = returnListForPlotting(dictToUse,keysToUse,indsMain[1])
    
    # get lists with plotting info
    returnList_3 = returnListForPlotting(dictToUse,keysToUse,indsMain[2])
    
    # get lists with plotting info
    returnList_4 = returnListForPlotting(dictToUse,keysToUse,indsMain[3])
    
    # get lists with plotting info
    returnList_5 = returnListForPlotting(dictToUse,keysToUse,indsMain[4])


    # do plotting
    
    yExt1,x1,y1 = doPlotCondition(returnList_1[6:11],col1,legEntries[0],lw,alphaVal,errorFormat)
            
    yExt2,x2,y2 = doPlotCondition(returnList_2[0:5],col2,legEntries[1],lw,alphaVal,errorFormat)
        
    yExt3,x3,y3 = doPlotCondition(returnList_3[0:5],col3,legEntries[2],lw,alphaVal,errorFormat)
    
    yExt4,x4,y4 = doPlotCondition(returnList_4[0:5],col4,legEntries[3],lw,alphaVal,errorFormat)

    yExt5,x5,y5 = doPlotCondition(returnList_5[0:5],col5,legEntries[4],lw,alphaVal,errorFormat)

    # get min and max y
    yExt = np.concatenate((yExt1,yExt2),axis=0)
    yMin = np.min(yExt)
    yMax = np.max(yExt)
    
    # get array for x values
    x = x1
    
    # for reference curves
    xMin = 0
    xMax = -1
    #multFactor = yMax
    
    # adjust start
    multFactorX = 1
    multFactorY1 = .07
    multFactorY2 = .3
    multFactorY3 = .032
    #multFactorY1 = .1
    #multFactorY2 = 1
    xDelta = .1
    
    # adjusts end of ref line
    xMaxOffset1 = 0
    xMaxOffset2 = .8
    xMaxOffset3 = 0
    
    xMinScale = .6
    
    xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
    xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
    xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX
    
    #b = np.min(x)+.1
    #b=10e-2
    #b=yMax
    #b=xtest1[0]-(xtest1[0]-.001)
    b=0
    #xtest2[0] = x[0]
    ytest1 = (xtest1+b)*multFactorY1
    #ytest1 = ytest1-ytest1[0]+y[0]+b
    ytest2 = (np.square(xtest2)+b)*multFactorY2
    #ytest2 = ytest2-ytest2[0]+y[0]+b
    ytest3 = (np.power(xtest3,0.5))*multFactorY3
    
    iSt = 0;
    plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
    plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
    plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)
        
    plt.xscale('log')
    plt.yscale('log')
    
    plt.xlabel(r'$\tau$ (s)',fontsize=afs)
    plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)
    
    plt.xticks(fontsize=tfs)
    plt.yticks(fontsize=tfs)
    plt.xlim(x[0]-.05,x[-1]+1)
    plt.ylim(.001,.2)
    
    
    
    if doSave:
        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'.pdf',bbox_inches="tight")
    #...
    
    plt.legend(loc=4,frameon=False,fontsize=lfs)

    if doSave:
        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'_withLegends.pdf',bbox_inches="tight")
    #...
#    plt.clf()
    
#%%
#
#
#
#
    
#%% MreB# and mtirf35
infoLists  =  [
        [[96,30],['MreB Protoplasts','MreB Cells','FloA','FloT']],
        ]

doSave = 1
saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'

# legend font size
lfs = 10
# axis font size
afs = 24
# tick font size
tfs = 18
# line width
lw = 2

# choose plot colors
col4 = '#00fffb'
col3 = '#e600ff'
col1 = '#0004ff'
col2 = '#ff0000'

# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict

# pick opacity for shade plot
alphaVal = 0.25

# pick error format
errorFormat = 'ci'
#errorFormat = 'sem'
#errorFormat = 'std'

counter = 0
for infoList in infoLists:
    
    counter += 1
    
    print(infoList)
    print('info list:')
    print(counter)
    print(' ')
    # get indices 
    indsMain = infoList[0]
    
    # get legend entires
    legEntries = infoList[1]
    
    # define filename as string of legend entries
    fname = legEntries[0].replace(" ","")+'_'\
    +legEntries[1].replace(" ","")+'_'\
    +legEntries[2].replace(" ","")+'_'\
    +legEntries[3].replace(" ","")


#    returnList = [
#                    dict_vals[subKey0],
#                    dict_cis[subKey0],
#                    dict_cis2[subKey0],
#                    dict_stds[subKey0],
#                    dict_nSamples[subKey0],
#                    subKey0,
#                    ]

    # get lists with plotting info
    returnList_1 = returnListForPlotting(dictToUse,keysToUse,indsMain[0])
    
    # get lists with plotting info
    returnList_2 = returnListForPlotting(dictToUse,keysToUse,indsMain[1])

    # do plotting
    
    yExt1,x1,y1 = doPlotCondition(returnList_1[0:5],col1,legEntries[0],lw,alphaVal,errorFormat)
            
    yExt2,x2,y2 = doPlotCondition(returnList_1[6:11],col2,legEntries[1],lw,alphaVal,errorFormat)
        
    yExt3,x3,y3 = doPlotCondition(returnList_2[0:5],col3,legEntries[2],lw,alphaVal,errorFormat)
    
    yExt4,x4,y4 = doPlotCondition(returnList_2[6:11],col4,legEntries[3],lw,alphaVal,errorFormat)

    # get min and max y
    yExt = np.concatenate((yExt1,yExt2),axis=0)
    yMin = np.min(yExt)
    yMax = np.max(yExt)
    
    # get array for x values
    x = x1
    
    # for reference curves
    xMin = 0
    xMax = -1
    #multFactor = yMax
    
    # adjust start
    multFactorX = 1
    multFactorY1 = .07
    multFactorY2 = .3
    multFactorY3 = .032
    #multFactorY1 = .1
    #multFactorY2 = 1
    xDelta = .1
    
    # adjusts end of ref line
    xMaxOffset1 = 0
    xMaxOffset2 = .8
    xMaxOffset3 = 0
    
    xMinScale = .6
    
    xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
    xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
    xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX
    
    #b = np.min(x)+.1
    #b=10e-2
    #b=yMax
    #b=xtest1[0]-(xtest1[0]-.001)
    b=0
    #xtest2[0] = x[0]
    ytest1 = (xtest1+b)*multFactorY1
    #ytest1 = ytest1-ytest1[0]+y[0]+b
    ytest2 = (np.square(xtest2)+b)*multFactorY2
    #ytest2 = ytest2-ytest2[0]+y[0]+b
    ytest3 = (np.power(xtest3,0.5))*multFactorY3
    
    iSt = 0;
    plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
    plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
    plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)
    
    
    plt.xscale('log')
    plt.yscale('log')
    
    plt.xlabel(r'$\tau$ (s)',fontsize=afs)
    plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)
    
    plt.xticks(fontsize=tfs)
    plt.yticks(fontsize=tfs)
#    plt.xlim(x[0]-.05,x[-1]+1)
    plt.xlim(x[0]-.05,13)
    plt.ylim(.001,.2)
    
    
    
    if doSave:
        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'.pdf',bbox_inches="tight")
    #...
    
    plt.legend(loc=4,frameon=False,fontsize=lfs)

    if doSave:
        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'_withLegend.pdf',bbox_inches="tight")
    #...
    
#    plt.clf()


#
#
#
#
    
    
#%% MreB# and mtirf35
infoLists  =  [
        [[96,30],['MreB Cells','FloA','FloT']],
        ]

doSave = 0
saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'

# legend font size
lfs = 10
# axis font size
afs = 24
# tick font size
tfs = 18
# line width
lw = 2

# choose plot colors
col4 = '#00fffb'
col3 = '#e600ff'
col1 = '#0004ff'
col2 = '#ff0000'

# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict

# pick opacity for shade plot
alphaVal = 0.25

# pick error format
errorFormat = 'ci'
#errorFormat = 'sem'
#errorFormat = 'std'

counter = 0
for infoList in infoLists:
    
    counter += 1
    
    print(infoList)
    print('info list:')
    print(counter)
    print(' ')
    # get indices 
    indsMain = infoList[0]
    
    # get legend entires
    legEntries = infoList[1]
    
    # define filename as string of legend entries
    fname = legEntries[0].replace(" ","")+'_'\
    +legEntries[1].replace(" ","")+'_'\
    +legEntries[2].replace(" ","")


#    returnList = [
#                    dict_vals[subKey0],
#                    dict_cis[subKey0],
#                    dict_cis2[subKey0],
#                    dict_stds[subKey0],
#                    dict_nSamples[subKey0],
#                    subKey0,
#                    ]

    # get lists with plotting info
    returnList_1 = returnListForPlotting(dictToUse,keysToUse,indsMain[0])
    
    # get lists with plotting info
    returnList_2 = returnListForPlotting(dictToUse,keysToUse,indsMain[1])

    # do plotting
                
    yExt1,x1,y1 = doPlotCondition(returnList_1[6:11],colBlack,legEntries[0],lw,alphaVal,errorFormat)
        
    yExt2,x2,y2 = doPlotCondition(returnList_2[0:5],colMagenta,legEntries[1],lw,alphaVal,errorFormat)
    
    yExt3,x3,y3 = doPlotCondition(returnList_2[6:11],colCyan,legEntries[2],lw,alphaVal,errorFormat)

    # get min and max y
    yExt = np.concatenate((yExt1,yExt2),axis=0)
    yMin = np.min(yExt)
    yMax = np.max(yExt)
    
    # get array for x values
    x = x1
    
    # for reference curves
    xMin = 0
    xMax = -1
    #multFactor = yMax
    
    # adjust start
    multFactorX = 1
    multFactorY1 = .07
    multFactorY2 = .3
    multFactorY3 = .032
    #multFactorY1 = .1
    #multFactorY2 = 1
    xDelta = .1
    
    # adjusts end of ref line
    xMaxOffset1 = 0
    xMaxOffset2 = .8
    xMaxOffset3 = 0
    
    xMinScale = .6
    
    xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
    xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
    xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX
    
    #b = np.min(x)+.1
    #b=10e-2
    #b=yMax
    #b=xtest1[0]-(xtest1[0]-.001)
    b=0
    #xtest2[0] = x[0]
    ytest1 = (xtest1+b)*multFactorY1
    #ytest1 = ytest1-ytest1[0]+y[0]+b
    ytest2 = (np.square(xtest2)+b)*multFactorY2
    #ytest2 = ytest2-ytest2[0]+y[0]+b
    ytest3 = (np.power(xtest3,0.5))*multFactorY3
    
    iSt = 0;
    plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
    plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
    plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)
    
    
    plt.xscale('log')
    plt.yscale('log')
    
    plt.xlabel(r'$\tau$ (s)',fontsize=afs)
    plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)
    
    plt.xticks(fontsize=tfs)
    plt.yticks(fontsize=tfs)
#    plt.xlim(x[0]-.05,x[-1]+1)
    plt.xlim(x[0]-.05,13)
    plt.ylim(.001,.2)
    
    
    
    if doSave:
        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'.pdf',bbox_inches="tight")
    #...
    
    plt.legend(loc=4,frameon=False,fontsize=lfs)

    if doSave:
        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'_withLegend.pdf',bbox_inches="tight")
    #...
    
#    plt.clf()


#
#
#
#
#%% mtirf100
infoLists  =  [
        [[32],['FloA','FloT']],
        ]

doSave = 0
saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'

# legend font size
lfs = 10
# axis font size
afs = 24
# tick font size
tfs = 18
# line width
lw = 2

# choose plot colors
col1 = colMagenta
col2 = colCyan

# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict

# pick opacity for shade plot
alphaVal = 0.25

# pick error format
errorFormat = 'ci'
#errorFormat = 'sem'
#errorFormat = 'std'

counter = 0
for infoList in infoLists:
    
    counter += 1
    
    print(infoList)
    print('info list:')
    print(counter)
    print(' ')
    # get indices 
    indsMain = infoList[0]
    
    # get legend entires
    legEntries = infoList[1]
    
    # define filename as string of legend entries
    fname = legEntries[0].replace(" ","")+'_'\
    +legEntries[1].replace(" ","")

#    returnList = [
#                    dict_vals[subKey0],
#                    dict_cis[subKey0],
#                    dict_cis2[subKey0],
#                    dict_stds[subKey0],
#                    dict_nSamples[subKey0],
#                    subKey0,
#                    ]

    # get lists with plotting info
    returnList_1 = returnListForPlotting(dictToUse,keysToUse,indsMain[0])
    
#    # get lists with plotting info
#    returnList_2 = returnListForPlotting(dictToUse,keysToUse,indsMain[1])

    # do plotting
    
    yExt1,x1,y1 = doPlotCondition(returnList_1[0:5],col1,legEntries[0],lw,alphaVal,errorFormat)
            
    yExt2,x2,y2 = doPlotCondition(returnList_1[6:11],col2,legEntries[1],lw,alphaVal,errorFormat)
        
#    yExt3,x3,y3 = doPlotCondition(returnList_2[0:5],col3,legEntries[2],lw,alphaVal,errorFormat)
#    
#    yExt4,x4,y4 = doPlotCondition(returnList_2[6:11],col4,legEntries[3],lw,alphaVal,errorFormat)

    # get min and max y
    yExt = np.concatenate((yExt1,yExt2),axis=0)
    yMin = np.min(yExt)
    yMax = np.max(yExt)
    
    # get array for x values
    x = x1
    
    # for reference curves
    xMin = 0
    xMax = -1
    #multFactor = yMax
    
    # adjust start
    multFactorX = 1
    multFactorY1 = .07
    multFactorY2 = .3
    multFactorY3 = .032
    #multFactorY1 = .1
    #multFactorY2 = 1
    xDelta = .1
    
    # adjusts end of ref line
    xMaxOffset1 = 0
    xMaxOffset2 = .8
    xMaxOffset3 = 0
    
    xMinScale = .6
    
    xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
    xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
    xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX
    
    #b = np.min(x)+.1
    #b=10e-2
    #b=yMax
    #b=xtest1[0]-(xtest1[0]-.001)
    b=0
    #xtest2[0] = x[0]
    ytest1 = (xtest1+b)*multFactorY1
    #ytest1 = ytest1-ytest1[0]+y[0]+b
    ytest2 = (np.square(xtest2)+b)*multFactorY2
    #ytest2 = ytest2-ytest2[0]+y[0]+b
    ytest3 = (np.power(xtest3,0.5))*multFactorY3
    
    iSt = 0;
    plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
    plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
    plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)
    
    
    plt.xscale('log')
    plt.yscale('log')
    
    plt.xlabel(r'$\tau$ (s)',fontsize=afs)
    plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)
    
    plt.xticks(fontsize=tfs)
    plt.yticks(fontsize=tfs)
    plt.xlim(x[0]-.05,x[-1]+1)
    plt.ylim(.001,.2)
    
    
    
    if doSave:
        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'_tirf100.pdf',bbox_inches="tight")
    #...
    
    plt.legend(loc=4,frameon=False,fontsize=lfs)

    if doSave:
        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'_tirf100_withLegend.pdf',bbox_inches="tight")
    #...
    
#    plt.clf()
    
    
#%% all below DEPRECATED
#
#
#
#
#
#
#
#%% plotting DEPRECATED
    
# inds for m epi mad
m_epi_mad_ind = 0

# get the d
m_epi_mad_all = msd_dict[keys_msd_dict[m_epi_mad_ind]]

#%% plotting, control

doSave = 1
saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'
fname = 'control_epi_mad_FloA_FloT'

# legend font size
lfs = 14

# axis font size
afs = 24

# tick font size
tfs = 18

# line width
lw = 2

#mainDictList = dictToUse[keysToUse[indMain]]
# set the index of the relevant experiment
# for m epi mad control
indMain = 8
#indMain = 0 # for debugging

# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict

returnList_epiMadCtrl = returnListForPlotting(dictToUse,keysToUse,indMain)
    
x=returnList_epiMadCtrl[0]

# choose plot colors
col1 = '#00fffb'
col2 = '#e600ff'

# pick opacity for shade plot
alphaVal = 0.25

# pick error format
errorFormat = 'ci'
#errorFormat = 'sem'
#errorFormat = 'ci2'
#errorFormat = 'std'

# legend names
leg1 = 'FloT'
leg2 = 'FloA'

# do plotting

#yExt1,x1,y1 = doPlotCondition(returnList_epiMadCtrl[0:4],col1,leg1,lw,alphaVal,errorFormat)
yExt1,x1,y1 = doPlotCondition(returnList_epiMadCtrl[0:5],col1,leg1,lw,alphaVal,errorFormat)

#yExt2,x2,y2 = doPlotCondition(returnList_epiMadCtrl[5:9],col2,leg2,lw,alphaVal,errorFormat)
yExt2,x2,y2 = doPlotCondition(returnList_epiMadCtrl[6:11],col2,leg2,lw,alphaVal,errorFormat)

# get min and max y
yExt = np.concatenate((yExt1,yExt2),axis=0)
yMin = np.min(yExt)
yMax = np.max(yExt)

# get array for x values
x = x1


# for reference curves
xMin = 0
xMax = -1
#multFactor = yMax

# adjust start
multFactorX = 1
multFactorY1 = .07
multFactorY2 = .3
multFactorY3 = .032
#multFactorY1 = .1
#multFactorY2 = 1
xDelta = .1

# adjusts end of ref line
xMaxOffset1 = 0
xMaxOffset2 = .8
xMaxOffset3 = 0

xMinScale = .6

xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX

#b = np.min(x)+.1
#b=10e-2
#b=yMax
#b=xtest1[0]-(xtest1[0]-.001)
b=0
#xtest2[0] = x[0]
ytest1 = (xtest1+b)*multFactorY1
#ytest1 = ytest1-ytest1[0]+y[0]+b
ytest2 = (np.square(xtest2)+b)*multFactorY2
#ytest2 = ytest2-ytest2[0]+y[0]+b
ytest3 = (np.power(xtest3,0.5))*multFactorY3

iSt = 0;
plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)

plt.legend(loc=4,frameon=False,fontsize=lfs)

plt.xscale('log')
plt.yscale('log')

plt.xlabel(r'$\tau$ (s)',fontsize=afs)
plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)

plt.xticks(fontsize=tfs)
plt.yticks(fontsize=tfs)
plt.xlim(x[0]-.05,x[-1]+1)
plt.ylim(yMin-1,yMax+.065)

if doSave:
    plt.savefig(saveDirPlot+fname+'_'+errorFormat+'.pdf',bbox_inches="tight")

#...


#%% for benzyl alcohol expts
    
doSave = 1
saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'
fname = 'bnz_vs_control_epi_mad_FloA_FloT'

# legend font size
lfs = 10
# axis font size
afs = 24
# tick font size
tfs = 18
# line width
lw = 2


# set the index of the relevant experiment
# for m epi mad control
indMain = 8
# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict
returnList_epiMadCtrl = returnListForPlotting(dictToUse,keysToUse,indMain)

# for benzyl alcohol
indMain = 6
dictToUse = msd_dict
keysToUse = keys_msd_dict
returnList_bnz = returnListForPlotting(dictToUse,keysToUse,indMain)

# choose plot colors
col1 = '#00fffb'
col2 = '#e600ff'
col3 = '#0004ff'
col4 = '#ff0000'

# pick opacity for shade plot
alphaVal = 0.25

# pick error format
errorFormat = 'ci'
#errorFormat = 'sem'
#errorFormat = 'std'

# legend names
leg1 = 'FloT control'
leg2 = 'FloA control'
leg3 = 'FloT BNZ'
leg4 = 'FloA BNZ'

# do plotting

yExt1,x1,y1 = doPlotCondition(returnList_epiMadCtrl[0:5],col1,leg1,lw,alphaVal,errorFormat)
    
yExt2,x2,y2 = doPlotCondition(returnList_epiMadCtrl[6:11],col2,leg2,lw,alphaVal,errorFormat)

yExt3,x3,y3 = doPlotCondition(returnList_bnz[0:5],col3,leg3,lw,alphaVal,errorFormat)
    
yExt4,x4,y4 = doPlotCondition(returnList_bnz[6:11],col4,leg4,lw,alphaVal,errorFormat)



# get min and max y
yExt = np.concatenate((yExt1,yExt2,yExt3,yExt4),axis=0)
yMin = np.min(yExt)
yMax = np.max(yExt)

# get array for x values
x = x1



# for reference curves
xMin = 0
xMax = -1
#multFactor = yMax

# adjust start
multFactorX = 1
multFactorY1 = .07
multFactorY2 = .3
multFactorY3 = .032
#multFactorY1 = .1
#multFactorY2 = 1
xDelta = .1

# adjusts end of ref line
xMaxOffset1 = 0
xMaxOffset2 = .8
xMaxOffset3 = 0

xMinScale = .6

xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX

#b = np.min(x)+.1
#b=10e-2
#b=yMax
#b=xtest1[0]-(xtest1[0]-.001)
b=0
#xtest2[0] = x[0]
ytest1 = (xtest1+b)*multFactorY1
#ytest1 = ytest1-ytest1[0]+y[0]+b
ytest2 = (np.square(xtest2)+b)*multFactorY2
#ytest2 = ytest2-ytest2[0]+y[0]+b
ytest3 = (np.power(xtest3,0.5))*multFactorY3

iSt = 0;
plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)

plt.legend(loc=4,frameon=False,fontsize=lfs)

plt.xscale('log')
plt.yscale('log')

plt.xlabel(r'$\tau$ (s)',fontsize=afs)
plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)

plt.xticks(fontsize=tfs)
plt.yticks(fontsize=tfs)
plt.xlim(x[0]-.05,x[-1]+1)
plt.ylim(yMin-1,yMax+.065)



if doSave:
    plt.savefig(saveDirPlot+fname+'_'+errorFormat+'.pdf',bbox_inches="tight")

#...
    
    
#%% benzyl alcohol version 2

doSave = 1
saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'
fname = 'bnz_vs_mCtrlBNZ_vs_control_epi_mad_FloA_FloT'

# legend font size
lfs = 10
# axis font size
afs = 24
# tick font size
tfs = 18
# line width
lw = 2


# set the index of the relevant experiment

# for m epi mad control
indMain = 10
# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict
returnList_epiMadCtrl = returnListForPlotting(dictToUse,keysToUse,indMain)

# for mctrl bnz
indMain = 8
dictToUse = msd_dict
keysToUse = keys_msd_dict
returnList_CtrlBnz = returnListForPlotting(dictToUse,keysToUse,indMain)


# for benzyl alcohol
indMain = 6
dictToUse = msd_dict
keysToUse = keys_msd_dict
returnList_bnz = returnListForPlotting(dictToUse,keysToUse,indMain)

# choose plot colors
col1 = '#00fffb'
col2 = '#e600ff'
col3 = '#0004ff'
col4 = '#ff0000'
col5 = '#00ab03'
col6 = '#d68800'


# pick opacity for shade plot
alphaVal = 0.25

# pick error format
errorFormat = 'ci'
#errorFormat = 'sem'
#errorFormat = 'std'

# legend names
leg1 = 'FloT global control'
leg2 = 'FloA global control'
leg3 = 'FloT BNZ'
leg4 = 'FloA BNZ'
leg5 = 'FloA BNZ control'
leg6 = 'FloT BNZ control'

# do plotting

yExt1,x1,y1 = doPlotCondition(returnList_epiMadCtrl[0:5],col1,leg1,lw,alphaVal,errorFormat)
    
yExt2,x2,y2 = doPlotCondition(returnList_epiMadCtrl[6:11],col2,leg2,lw,alphaVal,errorFormat)

yExt5,x5,y5 = doPlotCondition(returnList_CtrlBnz[0:5],col5,leg5,lw,alphaVal,errorFormat)
    
yExt6,x6,y6 = doPlotCondition(returnList_CtrlBnz[6:11],col6,leg6,lw,alphaVal,errorFormat)

yExt3,x3,y3 = doPlotCondition(returnList_bnz[0:5],col3,leg3,lw,alphaVal,errorFormat)
    
yExt4,x4,y4 = doPlotCondition(returnList_bnz[6:11],col4,leg4,lw,alphaVal,errorFormat)



# get min and max y
yExt = np.concatenate((yExt1,yExt2,yExt3,yExt4),axis=0)
yMin = np.min(yExt)
yMax = np.max(yExt)

# get array for x values
x = x1



# for reference curves
xMin = 0
xMax = -1
#multFactor = yMax

# adjust start
multFactorX = 1
multFactorY1 = .07
multFactorY2 = .3
multFactorY3 = .032
#multFactorY1 = .1
#multFactorY2 = 1
xDelta = .1

# adjusts end of ref line
xMaxOffset1 = 0
xMaxOffset2 = .8
xMaxOffset3 = 0

xMinScale = .6

xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX

#b = np.min(x)+.1
#b=10e-2
#b=yMax
#b=xtest1[0]-(xtest1[0]-.001)
b=0
#xtest2[0] = x[0]
ytest1 = (xtest1+b)*multFactorY1
#ytest1 = ytest1-ytest1[0]+y[0]+b
ytest2 = (np.square(xtest2)+b)*multFactorY2
#ytest2 = ytest2-ytest2[0]+y[0]+b
ytest3 = (np.power(xtest3,0.5))*multFactorY3

iSt = 0;
plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)

plt.legend(loc=4,frameon=False,fontsize=lfs)

plt.xscale('log')
plt.yscale('log')

plt.xlabel(r'$\tau$ (s)',fontsize=afs)
plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)

plt.xticks(fontsize=tfs)
plt.yticks(fontsize=tfs)
plt.xlim(x[0]-.05,x[-1]+1)
plt.ylim(yMin-1,yMax+.065)



if doSave:
    plt.savefig(saveDirPlot+fname+'_'+errorFormat+'.pdf',bbox_inches="tight")

#...
    
#%% delta pbpc version 2

doSave = 1
saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'
fname = 'delta_pbpc_vs_mCtrlPbpc_control_epi_mad_FloA_FloT'


# legend font size
lfs = 10
# axis font size
afs = 24
# tick font size
tfs = 18
# line width
lw = 2


# set the index of the relevant experiment

# for m epi mad control
indMain = 10
# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict
returnList_epiMadCtrl = returnListForPlotting(dictToUse,keysToUse,indMain)

# for mctrl delta pbpc
indMain = 40
dictToUse = msd_dict
keysToUse = keys_msd_dict
returnList_CtrlDeltaPbpc = returnListForPlotting(dictToUse,keysToUse,indMain)


# for delta pbpc
indMain = 38
dictToUse = msd_dict
keysToUse = keys_msd_dict
returnList_delta_pbpc = returnListForPlotting(dictToUse,keysToUse,indMain)

# choose plot colors
col1 = '#00fffb'
col2 = '#e600ff'
col3 = '#0004ff'
col4 = '#ff0000'
col5 = '#00ab03'
col6 = '#d68800'


# pick opacity for shade plot
alphaVal = 0.25

# pick error format
errorFormat = 'ci'
#errorFormat = 'sem'
#errorFormat = 'std'

# legend names
leg1 = 'FloT global control'
leg2 = 'FloA global control'
leg3 = 'FloA ∆pbpC'
leg4 = 'FloT ∆pbpC'
leg5 = 'FloT ∆pbpC control'
leg6 = 'FloA ∆pbpC control'

# do plotting

yExt1,x1,y1 = doPlotCondition(returnList_epiMadCtrl[0:5],col1,leg1,lw,alphaVal,errorFormat)
    
yExt2,x2,y2 = doPlotCondition(returnList_epiMadCtrl[6:11],col2,leg2,lw,alphaVal,errorFormat)

yExt5,x5,y5 = doPlotCondition(returnList_CtrlDeltaPbpc[0:5],col5,leg5,lw,alphaVal,errorFormat)
    
yExt6,x6,y6 = doPlotCondition(returnList_CtrlDeltaPbpc[6:11],col6,leg6,lw,alphaVal,errorFormat)

yExt3,x3,y3 = doPlotCondition(returnList_delta_pbpc[0:5],col3,leg3,lw,alphaVal,errorFormat)
    
yExt4,x4,y4 = doPlotCondition(returnList_delta_pbpc[6:11],col4,leg4,lw,alphaVal,errorFormat)



# get min and max y
yExt = np.concatenate((yExt1,yExt2,yExt3,yExt4),axis=0)
yMin = np.min(yExt)
yMax = np.max(yExt)

# get array for x values
x = x1



# for reference curves
xMin = 0
xMax = -1
#multFactor = yMax

# adjust start
multFactorX = 1
multFactorY1 = .07
multFactorY2 = .3
multFactorY3 = .032
#multFactorY1 = .1
#multFactorY2 = 1
xDelta = .1

# adjusts end of ref line
xMaxOffset1 = 0
xMaxOffset2 = .8
xMaxOffset3 = 0

xMinScale = .6

xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX

#b = np.min(x)+.1
#b=10e-2
#b=yMax
#b=xtest1[0]-(xtest1[0]-.001)
b=0
#xtest2[0] = x[0]
ytest1 = (xtest1+b)*multFactorY1
#ytest1 = ytest1-ytest1[0]+y[0]+b
ytest2 = (np.square(xtest2)+b)*multFactorY2
#ytest2 = ytest2-ytest2[0]+y[0]+b
ytest3 = (np.power(xtest3,0.5))*multFactorY3

iSt = 0;
plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)

plt.legend(loc=4,frameon=False,fontsize=lfs)

plt.xscale('log')
plt.yscale('log')

plt.xlabel(r'$\tau$ (s)',fontsize=afs)
plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)

plt.xticks(fontsize=tfs)
plt.yticks(fontsize=tfs)
plt.xlim(x[0]-.05,x[-1]+1)
plt.ylim(yMin-1,yMax+.065)



if doSave:
    plt.savefig(saveDirPlot+fname+'_'+errorFormat+'.pdf',bbox_inches="tight")

#...
    
#%% protoplasts version 2

doSave = 1
saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'
fname = 'protoplasts_vs_mCtrlProtoplasts_control_epi_mad_FloA_FloT'


# legend font size
lfs = 10
# axis font size
afs = 24
# tick font size
tfs = 18
# line width
lw = 2


# set the index of the relevant experiment

# for m epi mad control
indMain = 10
# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict
returnList_epiMadCtrl = returnListForPlotting(dictToUse,keysToUse,indMain)

# for mctrl delta pbpc
indMain = 20
dictToUse = msd_dict
keysToUse = keys_msd_dict
returnList_CtrlProtoplasts = returnListForPlotting(dictToUse,keysToUse,indMain)


# for protoplasts
indMain = 18
dictToUse = msd_dict
keysToUse = keys_msd_dict
returnList_protoplasts = returnListForPlotting(dictToUse,keysToUse,indMain)

# choose plot colors
col1 = '#00fffb'
col2 = '#e600ff'
col3 = '#0004ff'
col4 = '#ff0000'
col5 = '#00ab03'
col6 = '#d68800'


# pick opacity for shade plot
alphaVal = 0.25

# pick error format
errorFormat = 'ci'
#errorFormat = 'sem'
#errorFormat = 'std'

# legend names
leg1 = 'FloT global control'
leg2 = 'FloA global control'
leg3 = 'FloT protoplasts'
leg4 = 'FloA protoplasts'
leg5 = 'FloA protoplasts control'
leg6 = 'FloT protoplasts control'

# do plotting

yExt1,x1,y1 = doPlotCondition(returnList_epiMadCtrl[0:5],col1,leg1,lw,alphaVal,errorFormat)
    
yExt2,x2,y2 = doPlotCondition(returnList_epiMadCtrl[6:11],col2,leg2,lw,alphaVal,errorFormat)

yExt5,x5,y5 = doPlotCondition(returnList_CtrlProtoplasts[0:5],col5,leg5,lw,alphaVal,errorFormat)
    
yExt6,x6,y6 = doPlotCondition(returnList_CtrlProtoplasts[6:11],col6,leg6,lw,alphaVal,errorFormat)

yExt3,x3,y3 = doPlotCondition(returnList_protoplasts[0:5],col3,leg3,lw,alphaVal,errorFormat)
    
yExt4,x4,y4 = doPlotCondition(returnList_protoplasts[6:11],col4,leg4,lw,alphaVal,errorFormat)



# get min and max y
yExt = np.concatenate((yExt1,yExt2,yExt3,yExt4),axis=0)
yMin = np.min(yExt)
yMax = np.max(yExt)

# get array for x values
x = x1



# for reference curves
xMin = 0
xMax = -1
#multFactor = yMax

# adjust start
multFactorX = 1
multFactorY1 = .07
multFactorY2 = .3
multFactorY3 = .032
#multFactorY1 = .1
#multFactorY2 = 1
xDelta = .1

# adjusts end of ref line
xMaxOffset1 = 0
xMaxOffset2 = .8
xMaxOffset3 = 0

xMinScale = .6

xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX

#b = np.min(x)+.1
#b=10e-2
#b=yMax
#b=xtest1[0]-(xtest1[0]-.001)
b=0
#xtest2[0] = x[0]
ytest1 = (xtest1+b)*multFactorY1
#ytest1 = ytest1-ytest1[0]+y[0]+b
ytest2 = (np.square(xtest2)+b)*multFactorY2
#ytest2 = ytest2-ytest2[0]+y[0]+b
ytest3 = (np.power(xtest3,0.5))*multFactorY3

iSt = 0;
plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)

plt.legend(loc=4,frameon=False,fontsize=lfs)

plt.xscale('log')
plt.yscale('log')

plt.xlabel(r'$\tau$ (s)',fontsize=afs)
plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)

plt.xticks(fontsize=tfs)
plt.yticks(fontsize=tfs)
plt.xlim(x[0]-.05,x[-1]+1)
plt.ylim(yMin-1,yMax+.065)



if doSave:
    plt.savefig(saveDirPlot+fname+'_'+errorFormat+'.pdf',bbox_inches="tight")

#...
    
#%% for delta mreb, bh, mbl version 2
    
doSave = 1
saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'
fname = 'delta_mreb_bh_mbl_vs_control_s_delta_rsgI'

# legend font size
lfs = 10
# axis font size
afs = 24
# tick font size
tfs = 18
# line width
lw = 2

# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict

# set the index of the relevant experiment

# for m epi mad control
indMain = 10

# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict
returnList_epiMadCtrl = returnListForPlotting(dictToUse,keysToUse,indMain)

# for m epi mad control
indMain = 54
returnList_s_delta_rsgI = returnListForPlotting(dictToUse,keysToUse,indMain)

indMain = 36
returnList_delta_mreb_bh_mbl = returnListForPlotting(dictToUse,keysToUse,indMain)


# choose plot colors
col1 = '#00fffb'
col2 = '#e600ff'
col3 = '#0004ff'
col4 = '#ff0000'
col5 = '#00ab03'
col6 = '#d68800'


# pick opacity for shade plot
alphaVal = 0.25

# pick error format
errorFormat = 'ci'
#errorFormat = 'sem'
#errorFormat = 'std'

# legend names
leg5 = 'FloT ∆rsgI'
leg6 = 'FloA ∆rsgI'
leg3 = 'FloT ∆mreB ∆mreBH ∆mbl'
leg4 = 'FloA ∆mreB ∆mreBH ∆mbl'
leg1 = 'FloT global control'
leg2 = 'FloA global control'

# do plotting
yExt1,x1,y1 = doPlotCondition(returnList_epiMadCtrl[0:5],col1,leg1,lw,alphaVal,errorFormat)
    
yExt2,x2,y2 = doPlotCondition(returnList_epiMadCtrl[6:11],col2,leg2,lw,alphaVal,errorFormat)

yExt5,x5,y5 = doPlotCondition(returnList_s_delta_rsgI[0:5],col5,leg5,lw,alphaVal,errorFormat)
    
yExt6,x6,y6 = doPlotCondition(returnList_s_delta_rsgI[6:11],col6,leg6,lw,alphaVal,errorFormat)

yExt3,x3,y3 = doPlotCondition(returnList_delta_mreb_bh_mbl[0:5],col3,leg3,lw,alphaVal,errorFormat)
    
yExt4,x4,y4 = doPlotCondition(returnList_delta_mreb_bh_mbl[6:11],col4,leg4,lw,alphaVal,errorFormat)

#%

# get min and max y
yExt = np.concatenate((yExt1,yExt2,yExt3,yExt4),axis=0)
yMin = np.min(yExt)
yMax = np.max(yExt)

# get array for x values
x = x1




# for reference curves
xMin = 0
xMax = -1
#multFactor = yMax

# adjust start
multFactorX = 1
multFactorY1 = .07
multFactorY2 = .3
multFactorY3 = .032
#multFactorY1 = .1
#multFactorY2 = 1
xDelta = .1

# adjusts end of ref line
xMaxOffset1 = 0
xMaxOffset2 = .8
xMaxOffset3 = 0

xMinScale = .6

xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX

#b = np.min(x)+.1
#b=10e-2
#b=yMax
#b=xtest1[0]-(xtest1[0]-.001)
b=0
#xtest2[0] = x[0]
ytest1 = (xtest1+b)*multFactorY1
#ytest1 = ytest1-ytest1[0]+y[0]+b
ytest2 = (np.square(xtest2)+b)*multFactorY2
#ytest2 = ytest2-ytest2[0]+y[0]+b
ytest3 = (np.power(xtest3,0.5))*multFactorY3

iSt = 0;
plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)

plt.legend(loc=4,frameon=False,fontsize=lfs)

plt.xscale('log')
plt.yscale('log')

plt.xlabel(r'$\tau$ (s)',fontsize=afs)
plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)

plt.xticks(fontsize=tfs)
plt.yticks(fontsize=tfs)
plt.xlim(x[0]-.05,x[-1]+1)
plt.ylim(yMin-1,yMax+.065)



if doSave:
    plt.savefig(saveDirPlot+fname+'_'+errorFormat+'.pdf',bbox_inches="tight")

#...
    

#%% for lbo
    
doSave = 1
saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'
fname = 'm_plus_lbo_vs_m_minus_lbo_vs_ctrl_epi_MAD'

# legend font size
lfs = 10
# axis font size
afs = 24
# tick font size
tfs = 18
# line width
lw = 2

# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict

# set the index of the relevant experiment
# for m epi mad control
indMain = 8
returnList_epiMadCtrl = returnListForPlotting(dictToUse,keysToUse,indMain)

indMain = 0
returnList_mPlusLbo = returnListForPlotting(dictToUse,keysToUse,indMain)

indMain = 2
returnList_mMinusLbo = returnListForPlotting(dictToUse,keysToUse,indMain)


# choose plot colors
col1 = '#00fffb'
col2 = '#e600ff'
col3 = '#0004ff'
col4 = '#ff0000'
col5 = '#00ab03'
col6 = '#d68800'

# pick opacity for shade plot
alphaVal = 0.25

# pick error format
errorFormat = 'ci'
#errorFormat = 'sem'
#errorFormat = 'std'

# legend names
leg1 = 'FloT control'
leg2 = 'FloA control'
leg3 = 'FloA +LBO2'
leg4 = 'FloT +LBO2'
leg5 = 'FloA -LBO2'
leg6 = 'FloT -LBO2'

# do plotting

yExt1,x1,y1 = doPlotCondition(returnList_epiMadCtrl[0:5],col1,leg1,lw,alphaVal,errorFormat)
    
yExt2,x2,y2 = doPlotCondition(returnList_epiMadCtrl[6:11],col2,leg2,lw,alphaVal,errorFormat)

yExt3,x3,y3 = doPlotCondition(returnList_mPlusLbo[0:5],col3,leg3,lw,alphaVal,errorFormat)
    
yExt4,x4,y4 = doPlotCondition(returnList_mPlusLbo[6:11],col4,leg4,lw,alphaVal,errorFormat)

yExt5,x5,y5 = doPlotCondition(returnList_mMinusLbo[0:5],col5,leg5,lw,alphaVal,errorFormat)
    
yExt6,x6,y6 = doPlotCondition(returnList_mMinusLbo[6:11],col6,leg6,lw,alphaVal,errorFormat)

#%

# get min and max y
yExt = np.concatenate((yExt1,yExt2,yExt3,yExt4),axis=0)
yMin = np.min(yExt)
yMax = np.max(yExt)

# get array for x values
x = x1




# for reference curves
xMin = 0
xMax = -1
#multFactor = yMax

# adjust start
multFactorX = 1
multFactorY1 = .07
multFactorY2 = .3
multFactorY3 = .032
#multFactorY1 = .1
#multFactorY2 = 1
xDelta = .1

# adjusts end of ref line
xMaxOffset1 = 0
xMaxOffset2 = .8
xMaxOffset3 = 0

xMinScale = .6

xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX

#b = np.min(x)+.1
#b=10e-2
#b=yMax
#b=xtest1[0]-(xtest1[0]-.001)
b=0
#xtest2[0] = x[0]
ytest1 = (xtest1+b)*multFactorY1
#ytest1 = ytest1-ytest1[0]+y[0]+b
ytest2 = (np.square(xtest2)+b)*multFactorY2
#ytest2 = ytest2-ytest2[0]+y[0]+b
ytest3 = (np.power(xtest3,0.5))*multFactorY3

iSt = 0;
plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)

plt.legend(loc=4,frameon=False,fontsize=lfs)

plt.xscale('log')
plt.yscale('log')

plt.xlabel(r'$\tau$ (s)',fontsize=afs)
plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)

plt.xticks(fontsize=tfs)
plt.yticks(fontsize=tfs)
plt.xlim(x[0]-.05,x[-1]+1)
plt.ylim(yMin-1,yMax+.065)




if doSave:
    print('SAVING')
    plt.savefig(saveDirPlot+fname+'_'+errorFormat+'.pdf',bbox_inches="tight")

#...
    
    

#%% for delta pbpc
    
doSave = 1
saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'
fname = 'delta_pbpc_vs_control_epi_mad_FloA_FloT'

# legend font size
lfs = 10
# axis font size
afs = 24
# tick font size
tfs = 18
# line width
lw = 2

# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict

# set the index of the relevant experiment
# for m epi mad control
indMain = 8
returnList_epiMadCtrl = returnListForPlotting(dictToUse,keysToUse,indMain)

returnList_ctrl

indMain = 38
returnList_delta_pbpC = returnListForPlotting(dictToUse,keysToUse,indMain)


# choose plot colors
col1 = '#00fffb'
col2 = '#e600ff'
col3 = '#0004ff'
col4 = '#ff0000'

# pick opacity for shade plot
alphaVal = 0.25

# pick error format
errorFormat = 'ci'
#errorFormat = 'sem'
#errorFormat = 'std'

# legend names
leg1 = 'FloT control'
leg2 = 'FloA control'
leg3 = 'FloA, ∆pbpC'
leg4 = 'FloT, ∆pbpC'

# do plotting

yExt1,x1,y1 = doPlotCondition(returnList_epiMadCtrl[0:5],col1,leg1,lw,alphaVal,errorFormat)
    
yExt2,x2,y2 = doPlotCondition(returnList_epiMadCtrl[6:11],col2,leg2,lw,alphaVal,errorFormat)

yExt3,x3,y3 = doPlotCondition(returnList_delta_pbpC[0:5],col3,leg3,lw,alphaVal,errorFormat)
    
yExt4,x4,y4 = doPlotCondition(returnList_delta_pbpC[6:11],col4,leg4,lw,alphaVal,errorFormat)

#%

# get min and max y
yExt = np.concatenate((yExt1,yExt2,yExt3,yExt4),axis=0)
yMin = np.min(yExt)
yMax = np.max(yExt)

# get array for x values
x = x1




# for reference curves
xMin = 0
xMax = -1
#multFactor = yMax

# adjust start
multFactorX = 1
multFactorY1 = .07
multFactorY2 = .3
multFactorY3 = .032
#multFactorY1 = .1
#multFactorY2 = 1
xDelta = .1

# adjusts end of ref line
xMaxOffset1 = 0
xMaxOffset2 = .8
xMaxOffset3 = 0

xMinScale = .6

xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX

#b = np.min(x)+.1
#b=10e-2
#b=yMax
#b=xtest1[0]-(xtest1[0]-.001)
b=0
#xtest2[0] = x[0]
ytest1 = (xtest1+b)*multFactorY1
#ytest1 = ytest1-ytest1[0]+y[0]+b
ytest2 = (np.square(xtest2)+b)*multFactorY2
#ytest2 = ytest2-ytest2[0]+y[0]+b
ytest3 = (np.power(xtest3,0.5))*multFactorY3

iSt = 0;
plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)

plt.legend(loc=4,frameon=False,fontsize=lfs)

plt.xscale('log')
plt.yscale('log')

plt.xlabel(r'$\tau$ (s)',fontsize=afs)
plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)

plt.xticks(fontsize=tfs)
plt.yticks(fontsize=tfs)
plt.xlim(x[0]-.05,x[-1]+1)
plt.ylim(yMin-1,yMax+.065)




if doSave:
    print('SAVING')
    plt.savefig(saveDirPlot+fname+'_'+errorFormat+'.pdf',bbox_inches="tight")

#...


#%% for chimeric constructs
    
doSave = 1
saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'
fname = 'chimeras_vs_control_epi_mad_FloA_FloT'

# legend font size
lfs = 10
# axis font size
afs = 24
# tick font size
tfs = 18
# line width
lw = 2

# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict

# set the index of the relevant experiment
# for m epi mad control
indMain = 8
returnList_epiMadCtrl = returnListForPlotting(dictToUse,keysToUse,indMain)

# for FloAntTct
indMain = 88
returnList_FloAntTct = returnListForPlotting(dictToUse,keysToUse,indMain)

# for FloTntAct
indMain = 90
returnList_FloTntAct = returnListForPlotting(dictToUse,keysToUse,indMain)



# choose plot colors
col1 = '#00fffb'
col2 = '#e600ff'
col3 = '#0004ff'
col4 = '#ff0000'

# pick opacity for shade plot
alphaVal = 0.25

# pick error format
errorFormat = 'ci'
#errorFormat = 'sem'
#errorFormat = 'std'

# legend names
leg1 = 'FloT control'
leg2 = 'FloA control'
leg3 = 'FloAntTct'
leg4 = 'FloTntAct'

# do plotting

yExt1,x1,y1 = doPlotCondition(returnList_epiMadCtrl[0:5],col1,leg1,lw,alphaVal,errorFormat)
    
yExt2,x2,y2 = doPlotCondition(returnList_epiMadCtrl[6:11],col2,leg2,lw,alphaVal,errorFormat)

yExt3,x3,y3 = doPlotCondition(returnList_FloAntTct[0:5],col3,leg3,lw,alphaVal,errorFormat)
    
yExt4,x4,y4 = doPlotCondition(returnList_FloTntAct[0:5],col4,leg4,lw,alphaVal,errorFormat)

#%

# get min and max y
yExt = np.concatenate((yExt1,yExt2,yExt3,yExt4),axis=0)
yMin = np.min(yExt)
yMax = np.max(yExt)

# get array for x values
x = x1




# for reference curves
xMin = 0
xMax = -1
#multFactor = yMax

# adjust start
multFactorX = 1
multFactorY1 = .07
multFactorY2 = .3
multFactorY3 = .032
#multFactorY1 = .1
#multFactorY2 = 1
xDelta = .1

# adjusts end of ref line
xMaxOffset1 = 0
xMaxOffset2 = .8
xMaxOffset3 = 0

xMinScale = .6

xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX

#b = np.min(x)+.1
#b=10e-2
#b=yMax
#b=xtest1[0]-(xtest1[0]-.001)
b=0
#xtest2[0] = x[0]
ytest1 = (xtest1+b)*multFactorY1
#ytest1 = ytest1-ytest1[0]+y[0]+b
ytest2 = (np.square(xtest2)+b)*multFactorY2
#ytest2 = ytest2-ytest2[0]+y[0]+b
ytest3 = (np.power(xtest3,0.5))*multFactorY3

iSt = 0;
plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)

plt.legend(loc=4,frameon=False,fontsize=lfs)

plt.xscale('log')
plt.yscale('log')

plt.xlabel(r'$\tau$ (s)',fontsize=afs)
plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)

plt.xticks(fontsize=tfs)
plt.yticks(fontsize=tfs)
plt.xlim(x[0]-.05,x[-1]+1)
plt.ylim(yMin-1,yMax+.065)




if doSave:
    plt.savefig(saveDirPlot+fname+'_'+errorFormat+'.pdf',bbox_inches="tight")

#...


#%% for delta dlt
    
doSave = 1
saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'
fname = 'delta_dlt_vs_control_epi_mad_FloA_FloT'

# legend font size
lfs = 10
# axis font size
afs = 24
# tick font size
tfs = 18
# line width
lw = 2

# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict

# set the index of the relevant experiment
# for m epi mad control
indMain = 8
returnList_epiMadCtrl = returnListForPlotting(dictToUse,keysToUse,indMain)

indMain = 28
returnList_delta_dlt = returnListForPlotting(dictToUse,keysToUse,indMain)


# choose plot colors
col1 = '#00fffb'
col2 = '#e600ff'
col3 = '#0004ff'
col4 = '#ff0000'

# pick opacity for shade plot
alphaVal = 0.25

# pick error format
errorFormat = 'ci'
#errorFormat = 'sem'
#errorFormat = 'std'

# legend names
leg1 = 'FloT control'
leg2 = 'FloA control'
leg3 = 'FloT, ∆dltA-E'
leg4 = 'FloA, ∆dltA-E'

# do plotting

yExt1,x1,y1 = doPlotCondition(returnList_epiMadCtrl[0:5],col1,leg1,lw,alphaVal,errorFormat)
    
yExt2,x2,y2 = doPlotCondition(returnList_epiMadCtrl[6:11],col2,leg2,lw,alphaVal,errorFormat)

yExt3,x3,y3 = doPlotCondition(returnList_delta_dlt[0:5],col3,leg3,lw,alphaVal,errorFormat)
    
yExt4,x4,y4 = doPlotCondition(returnList_delta_dlt[6:11],col4,leg4,lw,alphaVal,errorFormat)

#%

# get min and max y
yExt = np.concatenate((yExt1,yExt2,yExt3,yExt4),axis=0)
yMin = np.min(yExt)
yMax = np.max(yExt)

# get array for x values
x = x1




# for reference curves
xMin = 0
xMax = -1
#multFactor = yMax

# adjust start
multFactorX = 1
multFactorY1 = .07
multFactorY2 = .3
multFactorY3 = .032
#multFactorY1 = .1
#multFactorY2 = 1
xDelta = .1

# adjusts end of ref line
xMaxOffset1 = 0
xMaxOffset2 = .8
xMaxOffset3 = 0

xMinScale = .6

xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX

#b = np.min(x)+.1
#b=10e-2
#b=yMax
#b=xtest1[0]-(xtest1[0]-.001)
b=0
#xtest2[0] = x[0]
ytest1 = (xtest1+b)*multFactorY1
#ytest1 = ytest1-ytest1[0]+y[0]+b
ytest2 = (np.square(xtest2)+b)*multFactorY2
#ytest2 = ytest2-ytest2[0]+y[0]+b
ytest3 = (np.power(xtest3,0.5))*multFactorY3

iSt = 0;
plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)

plt.legend(loc=4,frameon=False,fontsize=lfs)

plt.xscale('log')
plt.yscale('log')

plt.xlabel(r'$\tau$ (s)',fontsize=afs)
plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)

plt.xticks(fontsize=tfs)
plt.yticks(fontsize=tfs)
plt.xlim(x[0]-.05,x[-1]+1)
plt.ylim(yMin-1,yMax+.065)




if doSave:
    print('SAVING')
    plt.savefig(saveDirPlot+fname+'_'+errorFormat+'.pdf',bbox_inches="tight")

#...
    
#%% for protoplasts
    
doSave = 1
saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'
fname = 'protoplasts_vs_control_epi_mad_FloA_FloT'

# legend font size
lfs = 10
# axis font size
afs = 24
# tick font size
tfs = 18
# line width
lw = 2

# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict

# set the index of the relevant experiment
# for m epi mad control
indMain = 8
returnList_epiMadCtrl = returnListForPlotting(dictToUse,keysToUse,indMain)

indMain = 16
returnList_protoplasts = returnListForPlotting(dictToUse,keysToUse,indMain)


# choose plot colors
col1 = '#00fffb'
col2 = '#e600ff'
col3 = '#0004ff'
col4 = '#ff0000'

# pick opacity for shade plot
alphaVal = 0.25

# pick error format
errorFormat = 'ci'
#errorFormat = 'sem'
#errorFormat = 'std'

# legend names
leg1 = 'FloT control'
leg2 = 'FloA control'
leg3 = 'FloT, protoplasts'
leg4 = 'FloA, protoplasts'

# do plotting

yExt1,x1,y1 = doPlotCondition(returnList_epiMadCtrl[0:5],col1,leg1,lw,alphaVal,errorFormat)
    
yExt2,x2,y2 = doPlotCondition(returnList_epiMadCtrl[6:11],col2,leg2,lw,alphaVal,errorFormat)

yExt3,x3,y3 = doPlotCondition(returnList_protoplasts[0:5],col3,leg3,lw,alphaVal,errorFormat)
    
yExt4,x4,y4 = doPlotCondition(returnList_protoplasts[6:11],col4,leg4,lw,alphaVal,errorFormat)

#%

# get min and max y
yExt = np.concatenate((yExt1,yExt2,yExt3,yExt4),axis=0)
yMin = np.min(yExt)
yMax = np.max(yExt)

# get array for x values
x = x1




# for reference curves
xMin = 0
xMax = -1
#multFactor = yMax

# adjust start
multFactorX = 1
multFactorY1 = .07
multFactorY2 = .3
multFactorY3 = .032
#multFactorY1 = .1
#multFactorY2 = 1
xDelta = .1

# adjusts end of ref line
xMaxOffset1 = 0
xMaxOffset2 = .8
xMaxOffset3 = 0

xMinScale = .6

xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX

#b = np.min(x)+.1
#b=10e-2
#b=yMax
#b=xtest1[0]-(xtest1[0]-.001)
b=0
#xtest2[0] = x[0]
ytest1 = (xtest1+b)*multFactorY1
#ytest1 = ytest1-ytest1[0]+y[0]+b
ytest2 = (np.square(xtest2)+b)*multFactorY2
#ytest2 = ytest2-ytest2[0]+y[0]+b
ytest3 = (np.power(xtest3,0.5))*multFactorY3

iSt = 0;
plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)

plt.legend(loc=4,frameon=False,fontsize=lfs)

plt.xscale('log')
plt.yscale('log')

plt.xlabel(r'$\tau$ (s)',fontsize=afs)
plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)

plt.xticks(fontsize=tfs)
plt.yticks(fontsize=tfs)
plt.xlim(x[0]-.05,x[-1]+1)
plt.ylim(yMin-1,yMax+.065)



if doSave:
    print('SAVING')
    plt.savefig(saveDirPlot+fname+'_'+errorFormat+'.pdf',bbox_inches="tight")

#...
    
#%% for delta mreb, bh, mbl
    
doSave = 1
saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'
fname = 'delta_mreb_bh_mbl_vs_control_epi_mad_FloA_FloT'

# legend font size
lfs = 10
# axis font size
afs = 24
# tick font size
tfs = 18
# line width
lw = 2

# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict

# set the index of the relevant experiment
# for m epi mad control
indMain = 8
returnList_epiMadCtrl = returnListForPlotting(dictToUse,keysToUse,indMain)

indMain = 32
returnList_delta_mreb_bh_mbl = returnListForPlotting(dictToUse,keysToUse,indMain)


# choose plot colors
col1 = '#00fffb'
col2 = '#e600ff'
col3 = '#0004ff'
col4 = '#ff0000'

# pick opacity for shade plot
alphaVal = 0.25

# pick error format
errorFormat = 'ci'
#errorFormat = 'sem'
#errorFormat = 'std'

# legend names
leg1 = 'FloT control'
leg2 = 'FloA control'
leg3 = 'FloT, ∆mreB, ∆mreBH, ∆mbl'
leg4 = 'FloA, ∆mreB, ∆mreBH, ∆mbl'

# do plotting

yExt1,x1,y1 = doPlotCondition(returnList_epiMadCtrl[0:5],col1,leg1,lw,alphaVal,errorFormat)
    
yExt2,x2,y2 = doPlotCondition(returnList_epiMadCtrl[6:11],col2,leg2,lw,alphaVal,errorFormat)

yExt3,x3,y3 = doPlotCondition(returnList_delta_mreb_bh_mbl[0:5],col3,leg3,lw,alphaVal,errorFormat)
    
yExt4,x4,y4 = doPlotCondition(returnList_delta_mreb_bh_mbl[6:11],col4,leg4,lw,alphaVal,errorFormat)

#%

# get min and max y
yExt = np.concatenate((yExt1,yExt2,yExt3,yExt4),axis=0)
yMin = np.min(yExt)
yMax = np.max(yExt)

# get array for x values
x = x1




# for reference curves
xMin = 0
xMax = -1
#multFactor = yMax

# adjust start
multFactorX = 1
multFactorY1 = .07
multFactorY2 = .3
multFactorY3 = .032
#multFactorY1 = .1
#multFactorY2 = 1
xDelta = .1

# adjusts end of ref line
xMaxOffset1 = 0
xMaxOffset2 = .8
xMaxOffset3 = 0

xMinScale = .6

xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX

#b = np.min(x)+.1
#b=10e-2
#b=yMax
#b=xtest1[0]-(xtest1[0]-.001)
b=0
#xtest2[0] = x[0]
ytest1 = (xtest1+b)*multFactorY1
#ytest1 = ytest1-ytest1[0]+y[0]+b
ytest2 = (np.square(xtest2)+b)*multFactorY2
#ytest2 = ytest2-ytest2[0]+y[0]+b
ytest3 = (np.power(xtest3,0.5))*multFactorY3

iSt = 0;
plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)

plt.legend(loc=4,frameon=False,fontsize=lfs)

plt.xscale('log')
plt.yscale('log')

plt.xlabel(r'$\tau$ (s)',fontsize=afs)
plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)

plt.xticks(fontsize=tfs)
plt.yticks(fontsize=tfs)
plt.xlim(x[0]-.05,x[-1]+1)
plt.ylim(yMin-1,yMax+.065)



if doSave:
    plt.savefig(saveDirPlot+fname+'_'+errorFormat+'.pdf',bbox_inches="tight")

#...
    
#%% plotting, control TIRF

doSave = 1
saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'
fname = 'control_tirf_FloA_FloT'

# legend font size
lfs = 14
# axis font size
afs = 24
# tick font size
tfs = 18
# line width
lw = 2


# set the index of the relevant experiment
indMain = 18

# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict

returnList_tirf100Ctrl = returnListForPlotting(dictToUse,keysToUse,indMain)
    
x=returnList_tirf100Ctrl[0]

# choose plot colors
col1 = '#00fffb'
col2 = '#e600ff'

# pick opacity for shade plot
alphaVal = 0.25

# pick error format
errorFormat = 'ci'
#errorFormat = 'sem'
#errorFormat = 'std'

# legend names
leg1 = 'FloT'
leg2 = 'FloA'

# do plotting

yExt1,x1,y1 = doPlotCondition(returnList_tirf100Ctrl[0:5],col1,leg1,lw,alphaVal,errorFormat)
    
yExt2,x2,y2 = doPlotCondition(returnList_tirf100Ctrl[6:11],col2,leg2,lw,alphaVal,errorFormat)

# get min and max y
yExt = np.concatenate((yExt1,yExt2),axis=0)
yMin = np.min(yExt)
yMax = np.max(yExt)

# get array for x values
x = x1



# for reference curves
xMin = 0
xMax = -1
#multFactor = yMax

# adjust start
multFactorX = 1
multFactorY1 = .07
multFactorY2 = .3
multFactorY3 = .032
#multFactorY1 = .1
#multFactorY2 = 1
xDelta = .1

# adjusts end of ref line
xMaxOffset1 = 0
xMaxOffset2 = .8
xMaxOffset3 = 0

xMinScale = .6

xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX

#b = np.min(x)+.1
#b=10e-2
#b=yMax
#b=xtest1[0]-(xtest1[0]-.001)
b=0
#xtest2[0] = x[0]
ytest1 = (xtest1+b)*multFactorY1
#ytest1 = ytest1-ytest1[0]+y[0]+b
ytest2 = (np.square(xtest2)+b)*multFactorY2
#ytest2 = ytest2-ytest2[0]+y[0]+b
ytest3 = (np.power(xtest3,0.5))*multFactorY3

iSt = 0;
plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)

plt.legend(loc=4,frameon=False,fontsize=lfs)

plt.xscale('log')
plt.yscale('log')

plt.xlabel(r'$\tau$ (s)',fontsize=afs)
plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)

plt.xticks(fontsize=tfs)
plt.yticks(fontsize=tfs)
plt.xlim(x[0]-.05,x[-1]+1)
plt.ylim(yMin-1,yMax+.065)

if doSave:
    print('SAVING')
    plt.savefig(saveDirPlot+fname+'_'+errorFormat+'.pdf',bbox_inches="tight")

#...


#%% for fos
    
doSave = 1
saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'
fname = 'fos_vs_control_epi_mad_FloA_FloT'

# legend font size
lfs = 10
# axis font size
afs = 24
# tick font size
tfs = 18
# line width
lw = 2

# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict

# set the index of the relevant experiment
# for m epi mad control
indMain = 8
returnList_epiMadCtrl = returnListForPlotting(dictToUse,keysToUse,indMain)

indMain = 12
returnList_fos = returnListForPlotting(dictToUse,keysToUse,indMain)


# choose plot colors
col1 = '#00fffb'
col2 = '#e600ff'
col3 = '#0004ff'
col4 = '#ff0000'

# pick opacity for shade plot
alphaVal = 0.25

# pick error format
errorFormat = 'ci'
#errorFormat = 'sem'
#errorFormat = 'std'

# legend names
leg1 = 'FloT control'
leg2 = 'FloA control'
leg3 = 'FloA, fos'
leg4 = 'FloT, fos'

# do plotting

yExt1,x1,y1 = doPlotCondition(returnList_epiMadCtrl[0:5],col1,leg1,lw,alphaVal,errorFormat)
    
yExt2,x2,y2 = doPlotCondition(returnList_epiMadCtrl[6:11],col2,leg2,lw,alphaVal,errorFormat)

yExt3,x3,y3 = doPlotCondition(returnList_fos[0:5],col3,leg3,lw,alphaVal,errorFormat)
    
yExt4,x4,y4 = doPlotCondition(returnList_fos[6:11],col4,leg4,lw,alphaVal,errorFormat)

#%

# get min and max y
yExt = np.concatenate((yExt1,yExt2,yExt3,yExt4),axis=0)
yMin = np.min(yExt)
yMax = np.max(yExt)

# get array for x values
x = x1




# for reference curves
xMin = 0
xMax = -1
#multFactor = yMax

# adjust start
multFactorX = 1
multFactorY1 = .07
multFactorY2 = .3
multFactorY3 = .032
#multFactorY1 = .1
#multFactorY2 = 1
xDelta = .1

# adjusts end of ref line
xMaxOffset1 = 0
xMaxOffset2 = .8
xMaxOffset3 = 0

xMinScale = .6

xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX

#b = np.min(x)+.1
#b=10e-2
#b=yMax
#b=xtest1[0]-(xtest1[0]-.001)
b=0
#xtest2[0] = x[0]
ytest1 = (xtest1+b)*multFactorY1
#ytest1 = ytest1-ytest1[0]+y[0]+b
ytest2 = (np.square(xtest2)+b)*multFactorY2
#ytest2 = ytest2-ytest2[0]+y[0]+b
ytest3 = (np.power(xtest3,0.5))*multFactorY3

iSt = 0;
plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)

plt.legend(loc=4,frameon=False,fontsize=lfs)

plt.xscale('log')
plt.yscale('log')

plt.xlabel(r'$\tau$ (s)',fontsize=afs)
plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)

plt.xticks(fontsize=tfs)
plt.yticks(fontsize=tfs)
plt.xlim(x[0]-.05,x[-1]+1)
plt.ylim(yMin-1,yMax+.065)



if doSave:
    print('SAVING')
    plt.savefig(saveDirPlot+fname+'_'+errorFormat+'.pdf',bbox_inches="tight")

#...
    
    
#%% for tun
    
doSave = 1
saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'
fname = 'tun_vs_control_epi_mad_FloA_FloT'

# legend font size
lfs = 10
# axis font size
afs = 24
# tick font size
tfs = 18
# line width
lw = 2

# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict

# set the index of the relevant experiment
# for m epi mad control
indMain = 8
returnList_epiMadCtrl = returnListForPlotting(dictToUse,keysToUse,indMain)

indMain = 22
returnList_tun = returnListForPlotting(dictToUse,keysToUse,indMain)


# choose plot colors
col1 = '#00fffb'
col2 = '#e600ff'
col3 = '#0004ff'
col4 = '#ff0000'

# pick opacity for shade plot
alphaVal = 0.25

# pick error format
errorFormat = 'ci'
#errorFormat = 'sem'
#errorFormat = 'std'

# legend names
leg1 = 'FloT control'
leg2 = 'FloA control'
leg3 = 'FloA, tun'
leg4 = 'FloT, tun'

# do plotting

yExt1,x1,y1 = doPlotCondition(returnList_epiMadCtrl[0:5],col1,leg1,lw,alphaVal,errorFormat)
    
yExt2,x2,y2 = doPlotCondition(returnList_epiMadCtrl[6:11],col2,leg2,lw,alphaVal,errorFormat)

yExt3,x3,y3 = doPlotCondition(returnList_tun[0:5],col3,leg3,lw,alphaVal,errorFormat)
    
yExt4,x4,y4 = doPlotCondition(returnList_tun[6:11],col4,leg4,lw,alphaVal,errorFormat)

#%

# get min and max y
yExt = np.concatenate((yExt1,yExt2,yExt3,yExt4),axis=0)
yMin = np.min(yExt)
yMax = np.max(yExt)

# get array for x values
x = x1




# for reference curves
xMin = 0
xMax = -1
#multFactor = yMax

# adjust start
multFactorX = 1
multFactorY1 = .07
multFactorY2 = .3
multFactorY3 = .032
#multFactorY1 = .1
#multFactorY2 = 1
xDelta = .1

# adjusts end of ref line
xMaxOffset1 = 0
xMaxOffset2 = .8
xMaxOffset3 = 0

xMinScale = .6

xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX

#b = np.min(x)+.1
#b=10e-2
#b=yMax
#b=xtest1[0]-(xtest1[0]-.001)
b=0
#xtest2[0] = x[0]
ytest1 = (xtest1+b)*multFactorY1
#ytest1 = ytest1-ytest1[0]+y[0]+b
ytest2 = (np.square(xtest2)+b)*multFactorY2
#ytest2 = ytest2-ytest2[0]+y[0]+b
ytest3 = (np.power(xtest3,0.5))*multFactorY3

iSt = 0;
plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)

plt.legend(loc=4,frameon=False,fontsize=lfs)

plt.xscale('log')
plt.yscale('log')

plt.xlabel(r'$\tau$ (s)',fontsize=afs)
plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)

plt.xticks(fontsize=tfs)
plt.yticks(fontsize=tfs)
plt.xlim(x[0]-.05,x[-1]+1)
plt.ylim(yMin-1,yMax+.065)



if doSave:
    print('SAVING')
    plt.savefig(saveDirPlot+fname+'_'+errorFormat+'.pdf',bbox_inches="tight")

#...


#%% for amp
    
doSave = 1
saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'
fname = 'amp_vs_control_epi_mad_FloA_FloT'

# legend font size
lfs = 10
# axis font size
afs = 24
# tick font size
tfs = 18
# line width
lw = 2

# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict

# set the index of the relevant experiment
# for m epi mad control
indMain = 8
returnList_epiMadCtrl = returnListForPlotting(dictToUse,keysToUse,indMain)

indMain = 4
returnList_amp = returnListForPlotting(dictToUse,keysToUse,indMain)


# choose plot colors
col1 = '#00fffb'
col2 = '#e600ff'
col3 = '#0004ff'
col4 = '#ff0000'

# pick opacity for shade plot
alphaVal = 0.25

# pick error format
errorFormat = 'ci'
#errorFormat = 'sem'
#errorFormat = 'std'

# legend names
leg1 = 'FloT control'
leg2 = 'FloA control'
leg3 = 'FloA, amp'
leg4 = 'FloT, amp'

# do plotting

yExt1,x1,y1 = doPlotCondition(returnList_epiMadCtrl[0:5],col1,leg1,lw,alphaVal,errorFormat)
    
yExt2,x2,y2 = doPlotCondition(returnList_epiMadCtrl[6:11],col2,leg2,lw,alphaVal,errorFormat)

yExt3,x3,y3 = doPlotCondition(returnList_amp[0:5],col3,leg3,lw,alphaVal,errorFormat)
    
yExt4,x4,y4 = doPlotCondition(returnList_amp[6:11],col4,leg4,lw,alphaVal,errorFormat)

#%

# get min and max y
yExt = np.concatenate((yExt1,yExt2,yExt3,yExt4),axis=0)
yMin = np.min(yExt)
yMax = np.max(yExt)

# get array for x values
x = x1




# for reference curves
xMin = 0
xMax = -1
#multFactor = yMax

# adjust start
multFactorX = 1
multFactorY1 = .07
multFactorY2 = .3
multFactorY3 = .032
#multFactorY1 = .1
#multFactorY2 = 1
xDelta = .1

# adjusts end of ref line
xMaxOffset1 = 0
xMaxOffset2 = .8
xMaxOffset3 = 0

xMinScale = .6

xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX

#b = np.min(x)+.1
#b=10e-2
#b=yMax
#b=xtest1[0]-(xtest1[0]-.001)
b=0
#xtest2[0] = x[0]
ytest1 = (xtest1+b)*multFactorY1
#ytest1 = ytest1-ytest1[0]+y[0]+b
ytest2 = (np.square(xtest2)+b)*multFactorY2
#ytest2 = ytest2-ytest2[0]+y[0]+b
ytest3 = (np.power(xtest3,0.5))*multFactorY3

iSt = 0;
plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)

plt.legend(loc=4,frameon=False,fontsize=lfs)

plt.xscale('log')
plt.yscale('log')

plt.xlabel(r'$\tau$ (s)',fontsize=afs)
plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)

plt.xticks(fontsize=tfs)
plt.yticks(fontsize=tfs)
plt.xlim(x[0]-.05,x[-1]+1)
plt.ylim(yMin-1,yMax+.065)



if doSave:
    print('SAVING')
    plt.savefig(saveDirPlot+fname+'_'+errorFormat+'.pdf',bbox_inches="tight")

#...
    
    
#%% for amp
    
doSave = 1
saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'
fname = 'van_vs_control_epi_mad_FloA_FloT'

# legend font size
lfs = 10
# axis font size
afs = 24
# tick font size
tfs = 18
# line width
lw = 2

# choose the dictionary and associated list of keys(msd or vac)
dictToUse = msd_dict
keysToUse = keys_msd_dict

# set the index of the relevant experiment
# for m epi mad control
indMain = 8
returnList_epiMadCtrl = returnListForPlotting(dictToUse,keysToUse,indMain)

indMain = 26
returnList_van = returnListForPlotting(dictToUse,keysToUse,indMain)


# choose plot colors
col1 = '#00fffb'
col2 = '#e600ff'
col3 = '#0004ff'
col4 = '#ff0000'

# pick opacity for shade plot
alphaVal = 0.25

# pick error format
errorFormat = 'ci'
#errorFormat = 'sem'
#errorFormat = 'std'

# legend names
leg1 = 'FloT control'
leg2 = 'FloA control'
leg3 = 'FloA, van'
leg4 = 'FloT, van'

# do plotting

yExt1,x1,y1 = doPlotCondition(returnList_epiMadCtrl[0:5],col1,leg1,lw,alphaVal,errorFormat)
    
yExt2,x2,y2 = doPlotCondition(returnList_epiMadCtrl[6:11],col2,leg2,lw,alphaVal,errorFormat)

yExt3,x3,y3 = doPlotCondition(returnList_van[0:5],col3,leg3,lw,alphaVal,errorFormat)
    
yExt4,x4,y4 = doPlotCondition(returnList_van[6:11],col4,leg4,lw,alphaVal,errorFormat)

#%

# get min and max y
yExt = np.concatenate((yExt1,yExt2,yExt3,yExt4),axis=0)
yMin = np.min(yExt)
yMax = np.max(yExt)

# get array for x values
x = x1




# for reference curves
xMin = 0
xMax = -1
#multFactor = yMax

# adjust start
multFactorX = 1
multFactorY1 = .07
multFactorY2 = .3
multFactorY3 = .032
#multFactorY1 = .1
#multFactorY2 = 1
xDelta = .1

# adjusts end of ref line
xMaxOffset1 = 0
xMaxOffset2 = .8
xMaxOffset3 = 0

xMinScale = .6

xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX

#b = np.min(x)+.1
#b=10e-2
#b=yMax
#b=xtest1[0]-(xtest1[0]-.001)
b=0
#xtest2[0] = x[0]
ytest1 = (xtest1+b)*multFactorY1
#ytest1 = ytest1-ytest1[0]+y[0]+b
ytest2 = (np.square(xtest2)+b)*multFactorY2
#ytest2 = ytest2-ytest2[0]+y[0]+b
ytest3 = (np.power(xtest3,0.5))*multFactorY3

iSt = 0;
plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)

plt.legend(loc=4,frameon=False,fontsize=lfs)

plt.xscale('log')
plt.yscale('log')

plt.xlabel(r'$\tau$ (s)',fontsize=afs)
plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)

plt.xticks(fontsize=tfs)
plt.yticks(fontsize=tfs)
plt.xlim(x[0]-.05,x[-1]+1)
plt.ylim(yMin-1,yMax+.065)



if doSave:
    print('SAVING')
    plt.savefig(saveDirPlot+fname+'_'+errorFormat+'.pdf',bbox_inches="tight")

#...
    
#%% for other data/supplemental figures
    
#indsAlreadyUsed = [8, 6, 54, 56, 34, 28, 16, 32, 18, 12, 22, 4, 26]


doSave = 1
saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'
for i in range(0,len(keysToUse),2):
#    if i in indsAlreadyUsed:
#        continue
    #...
    print(i)
    
#    doSave = 1
#    saveDirPlot = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/figures/'
    domainsTracked = keysToUse[i+1]
    if len(domainsTracked) == 2:
        if 'FloA' in domainsTracked[0]:
            leg3 = 'FloA '+keysToUse[i][1:]
            
        elif 'FloT' in domainsTracked[0]:
            leg3 = 'FloT '+keysToUse[i][1:]
            
        else:
            leg3 = domainsTracked[0].split('_')[0]+' '+keysToUse[i]
            
        if 'FloA' in domainsTracked[1]:
            leg4 = 'FloA '+keysToUse[i][1:]
            
        elif 'FloT' in domainsTracked[1]:
            leg4 = 'FloT '+keysToUse[i][1:]
            
        else:
            leg4 = domainsTracked[1].split('_')[0]+' '+keysToUse[i]
            
        #...
        
    elif len(domainsTracked) == 1:
        if 'FloA' in domainsTracked[0]:
            leg3 = 'FloA '+keysToUse[i][1:]
            
        elif 'FloT' in domainsTracked[0]:
            leg3 = 'FloA '+keysToUse[i][1:]
            
        else:
            leg3 = domainsTracked[0].split('_')[0]+' '+keysToUse[i]
            
        #...
        
    #...
    
    fname = keysToUse[i]+'_vs_control_epi_mad_FloA_FloT'
    
    # legend font size
    lfs = 10
    # axis font size
    afs = 24
    # tick font size
    tfs = 18
    # line width
    lw = 2
    
    
    # set the index of the relevant experiment
    # for m epi mad control
    indMain = 10
    # choose the dictionary and associated list of keys(msd or vac)
    dictToUse = msd_dict
    keysToUse = keys_msd_dict
    returnList_epiMadCtrl = returnListForPlotting(dictToUse,keysToUse,indMain)
    
    # for experimental condition
    dictToUse = msd_dict
    keysToUse = keys_msd_dict
    returnList_expt = returnListForPlotting(dictToUse,keysToUse,i)
    
    # choose plot colors
    col1 = '#00fffb'
    col2 = '#e600ff'
    col3 = '#0004ff'
    col4 = '#ff0000'
    
    # pick opacity for shade plot
    alphaVal = 0.25
    
    # pick error format
    errorFormat = 'ci'
    #errorFormat = 'sem'
    #errorFormat = 'std'
    
    # legend names
    leg1 = 'FloT control'
    leg2 = 'FloA control'

    
    # do plotting
    plt.figure()

    yExt1,x1,y1 = doPlotCondition(returnList_epiMadCtrl[0:5],col1,leg1,lw,alphaVal,errorFormat)
        
    yExt2,x2,y2 = doPlotCondition(returnList_epiMadCtrl[6:11],col2,leg2,lw,alphaVal,errorFormat)
    
    yExt3,x3,y3 = doPlotCondition(returnList_expt[0:5],col3,leg3,lw,alphaVal,errorFormat)
    
    if len(domainsTracked) == 2:
        yExt4,x4,y4 = doPlotCondition(returnList_expt[6:11],col4,leg4,lw,alphaVal,errorFormat)
    #...
    
    
    # get min and max y
    if len(domainsTracked) == 2:
        yExt = np.concatenate((yExt1,yExt2,yExt3,yExt4),axis=0)
        
    else:
        yExt = np.concatenate((yExt1,yExt2,yExt3),axis=0)
        
    #...
    
    yMin = np.min(yExt)
    yMax = np.max(yExt)
    
    # get array for x values
    x = x1
    
    
    
    # for reference curves
    xMin = 0
    xMax = -1
    #multFactor = yMax
    
    # adjust start
    multFactorX = 1
    multFactorY1 = .07
    multFactorY2 = .3
    multFactorY3 = .032
    #multFactorY1 = .1
    #multFactorY2 = 1
    xDelta = .1
    
    # adjusts end of ref line
    xMaxOffset1 = 0
    xMaxOffset2 = .8
    xMaxOffset3 = 0
    
    xMinScale = .6
    
    xtest1 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset1,xDelta))*multFactorX
    xtest2 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset2,xDelta))*multFactorX
    xtest3 = (np.arange(x[xMin]*xMinScale,x[xMax]-xMaxOffset3,xDelta))*multFactorX
    
    #b = np.min(x)+.1
    #b=10e-2
    #b=yMax
    #b=xtest1[0]-(xtest1[0]-.001)
    b=0
    #xtest2[0] = x[0]
    ytest1 = (xtest1+b)*multFactorY1
    #ytest1 = ytest1-ytest1[0]+y[0]+b
    ytest2 = (np.square(xtest2)+b)*multFactorY2
    #ytest2 = ytest2-ytest2[0]+y[0]+b
    ytest3 = (np.power(xtest3,0.5))*multFactorY3
    
    iSt = 0;
    plt.plot(xtest3[iSt:],ytest3[iSt:],'--',label='Slope = 0.5',linewidth=lw)
    plt.plot(xtest1[iSt:],ytest1[iSt:],'--',label='Slope = 1',linewidth=lw)
    plt.plot(xtest2[iSt:],ytest2[iSt:],'--',label='Slope = 2',linewidth=lw)
    
    plt.legend(loc=4,frameon=False,fontsize=lfs)
    
    plt.xscale('log')
    plt.yscale('log')
    
    plt.xlabel(r'$\tau$ (s)',fontsize=afs)
    plt.ylabel(r'MSD ($\mu$m$^2$)',fontsize=afs)
    
    plt.xticks(fontsize=tfs)
    plt.yticks(fontsize=tfs)
    plt.xlim(x[0]-.05,x[-1]+1)
    plt.ylim(yMin-1,yMax+.065)

    if doSave:
        plt.savefig(saveDirPlot+fname+'_'+errorFormat+'.pdf',bbox_inches="tight")
    #...
#    plt.clf()
#...