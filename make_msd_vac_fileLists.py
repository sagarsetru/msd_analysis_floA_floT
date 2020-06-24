#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 12:25:34 2020

@author: sagarsetru
"""

import os

startDir = '/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis'

os.chdir(startDir)

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
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/MreB'
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

#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/DltD/Amp',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/DltD/Cells',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/DltD/Ctrl',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/DltD/Fos',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/DltD/Protoplasts',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/DltD/Tun',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/DltD/Van',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/DltD/∆floA',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/DltD/∆flot',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/PBP3-2/Amp',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/PBP3-2/Cells',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/PBP3-2/Ctrl',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/PBP3-2/Fos',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/PBP3-2/Protoplasts',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/PBP3-2/Tun',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/PBP3-2/Van',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/PBP3-2/∆floA',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/PBP3-2/∆flot',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/MreB',
#'/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/PBP3',
#]


#%%
#find ~/Documents/Princeton/wingreen/for\ NED/analysis/Flotillins/m∆pbpC/FloT/* -type f -name "*.csv" | grep "vac_track" > m_delta_pbpC_FloT_vacTracks.txt

# define strings
str1 = 'find '
str2 = '/* -type f -name "*.csv" | grep '
str_vacFind = '"vac_track" > '
str_msdFind = '"msd_track" > '
str_FloT_dir = '/FloT'
str_FloA_dir = '/FloA'
str_FloA_fName = '_FloA'
str_FloT_fName = '_FloT'
str_vac_fName = '_vacTracks'
str_msd_fName = '_msdTracks'
str_ext = '.txt'

# loop over main directories
for mainDir in dirs_all:
    
    # get the experimental condition for this data
    exptCond = '/'+os.path.split(mainDir)[1]
    
    # go to appropriate directory
    #os.chdir(mainDir)
    
    # add slash for spaces in string
    mainDirRun = mainDir
    mainDirRun = mainDirRun.replace(" ", "\\ ")

    print(mainDir)
    
    # define conditions
    cond_FloA = os.path.exists(mainDir+str_FloA_dir)
    cond_FloT = os.path.exists(mainDir+str_FloT_dir)
    cond_Cells = os.path.exists(mainDir+'/Cells')
    cond_Protoplasts = os.path.exists(mainDir+'/Protoplasts')
    cond_Amp = os.path.exists(mainDir+'/Amp')
    cond_Ctrl = os.path.exists(mainDir+'/Ctrl')
    cond_Fos = os.path.exists(mainDir+'/Fos')
    cond_Tun = os.path.exists(mainDir+'/Tun')
    cond_Van = os.path.exists(mainDir+'/Van')
    cond_delta_floA = os.path.exists(mainDir+'/∆floA')
    cond_delta_floT = os.path.exists(mainDir+'/∆floT')
    
    cond_main = cond_FloA or cond_FloT or cond_Cells or cond_Protoplasts or cond_Amp or cond_Ctrl or cond_Fos or cond_Tun or cond_Van or cond_delta_floA or cond_delta_floT
    
    
    if cond_main:
    
        # if there is a subdirectory specific to "FloA"
        if cond_FloA:
            print('FloA')
            
            # define the command for vac and run
            command = str1+mainDirRun+str_FloA_dir+str2+str_vacFind+mainDirRun+exptCond+str_FloA_fName+str_vac_fName+str_ext
            os.system(command)
            
            # define the command for msd and run
            command = str1+mainDirRun+str_FloA_dir+str2+str_msdFind+mainDirRun+exptCond+str_FloA_fName+str_msd_fName+str_ext
            os.system(command)
        
        #...
            
        # if there is a subdirectory specific for "FloT"
        if cond_FloT:
            print('FloT')
            
            # define the command for vac and run
            command = str1+mainDirRun+str_FloT_dir+str2+str_vacFind+mainDirRun+exptCond+str_FloT_fName+str_vac_fName+str_ext
            os.system(command)
            
            # define the command for msd and run
            command = str1+mainDirRun+str_FloT_dir+str2+str_msdFind+mainDirRun+exptCond+str_FloT_fName+str_msd_fName+str_ext
            os.system(command)

        # for MreB, dltd, pbp3-2 folder
        if cond_Cells:
            print('Cells')
            
            # command for vac and run
            command = str1+mainDirRun+'/Cells'+str2+str_vacFind+mainDirRun+'/Cells'+str_vac_fName+str_ext
            os.system(command)
            
            # command for msd and run
            command = str1+mainDirRun+'/Cells'+str2+str_msdFind+mainDirRun+'/Cells'+str_msd_fName+str_ext
            os.system(command)
            
        #...
        
        # for MreB, dltd, pbp3-2 folder
        if cond_Protoplasts:
            print('Protoplasts')
            
            # command for vac and run
            command = str1+mainDirRun+'/Protoplasts'+str2+str_vacFind+mainDirRun+'/Protoplasts'+str_vac_fName+str_ext
            os.system(command)
            
            # command for msd and run
            command = str1+mainDirRun+'/Protoplasts'+str2+str_msdFind+mainDirRun+'/Protoplasts'+str_msd_fName+str_ext
            os.system(command)
        
        #...
        
        # for MreB, dltd, pbp3-2 folder
        if cond_Amp:
            print('Amp')
            
            # command for vac and run
            command = str1+mainDirRun+'/Amp'+str2+str_vacFind+mainDirRun+'/Amp'+str_vac_fName+str_ext
            os.system(command)
            
            # command for msd and run
            command = str1+mainDirRun+'/Amp'+str2+str_msdFind+mainDirRun+'/Amp'+str_msd_fName+str_ext
            os.system(command)
        
        #...
        
        # for MreB, dltd, pbp3-2 folder
        if cond_Ctrl:
            print('Ctrl')
            
            # command for vac and run
            command = str1+mainDirRun+'/Ctrl'+str2+str_vacFind+mainDirRun+'/Ctrl'+str_vac_fName+str_ext
            os.system(command)
            
            # command for msd and run
            command = str1+mainDirRun+'/Ctrl'+str2+str_msdFind+mainDirRun+'/Ctrl'+str_msd_fName+str_ext
            os.system(command)
        
        #...
        
        # for MreB, dltd, pbp3-2 folder
        if cond_Fos:
            print('Fos')
            
            # command for vac and run
            command = str1+mainDirRun+'/Fos'+str2+str_vacFind+mainDirRun+'/Fos'+str_vac_fName+str_ext
            os.system(command)
            
            # command for msd and run
            command = str1+mainDirRun+'/Fos'+str2+str_msdFind+mainDirRun+'/Fos'+str_msd_fName+str_ext
            os.system(command)
        
        #...
        
        # for MreB, dltd, pbp3-2 folder
        if cond_Tun:
            print('Tun')
            
            # command for vac and run
            command = str1+mainDirRun+'/Tun'+str2+str_vacFind+mainDirRun+'/Tun'+str_vac_fName+str_ext
            os.system(command)
            
            # command for msd and run
            command = str1+mainDirRun+'/Tun'+str2+str_msdFind+mainDirRun+'/Tun'+str_msd_fName+str_ext
            os.system(command)
        
        #...
        
        # for MreB, dltd, pbp3-2 folder
        if cond_Van:
            print('Van')
            
            # command for vac and run
            command = str1+mainDirRun+'/Van'+str2+str_vacFind+mainDirRun+'/Van'+str_vac_fName+str_ext
            os.system(command)
            
            # command for msd and run
            command = str1+mainDirRun+'/Van'+str2+str_msdFind+mainDirRun+'/Van'+str_msd_fName+str_ext
            os.system(command)
        
        #...
        
        # for MreB, dltd, pbp3-2 folder
        if cond_delta_floA:
            print('∆floA')
            
            # command for vac and run
            command = str1+mainDirRun+'/∆floA'+str2+str_vacFind+mainDirRun+'/∆floA'+str_vac_fName+str_ext
            os.system(command)
            
            # command for msd and run
            command = str1+mainDirRun+'/∆floA'+str2+str_msdFind+mainDirRun+'/∆floA'+str_msd_fName+str_ext
            os.system(command)
        
        #...
        
        # for MreB, dltd, pbp3-2 folder
        if cond_delta_floT:
            print('∆floT')
            
            # command for vac and run
            command = str1+mainDirRun+'/∆floT'+str2+str_vacFind+mainDirRun+'/∆floT'+str_vac_fName+str_ext
            os.system(command)
            
            # command for msd and run
            command = str1+mainDirRun+'/∆floT'+str2+str_msdFind+mainDirRun+'/∆floT'+str_msd_fName+str_ext
            os.system(command)
        
        #...
    
    # for fusion constructs
    else:
        print('no sub category')
        
        # define the command for vac and run
        command = str1+mainDirRun+str2+str_vacFind+mainDirRun+exptCond+str_vac_fName+str_ext
        os.system(command)
        
        # define the command for msd and run
        command = str1+mainDirRun+str2+str_msdFind+mainDirRun+exptCond+str_msd_fName+str_ext
        os.system(command)
        
    #...
    print(' ')
    
#...