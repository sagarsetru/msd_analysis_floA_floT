# Sagar Setru
# 2016 07 14

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import collections
import itertools
import bootstrapped.bootstrap as bs
import bootstrapped.stats_functions as bs_stats


#def loadDataFile( baseDir, subBaseDir, fileName):
#    "helper function to load bed files"
#    return pd.read_csv(baseDir+subBaseDir+'/'+fileName, sep='\t', header=None)
##...

def compute_MSD(path):
   totalsize=len(path)
   msd=[]
   for i in range(totalsize-1):
       j=i+1
       msd.append(np.sum((path[0:-j]-path[j::])**2)/float(totalsize-j))
    #...
   msd=np.array(msd)
   return msd
#...
#http://stackoverflow.com/questions/26472653/computing-the-mean-square-displacement-of-a-2d-random-walk-in-python

def compute_vel(path,t):
    totalsize=len(path)
    velx=[]
    vely=[]
    for i in range(totalsize-1):
        j=i+1
        delta_t = t[j]-t[i]
        velx.append( (path[j][0]-path[i][0]) / delta_t)
        vely.append( (path[j][1]-path[i][1]) / delta_t)
    #...
    velx=np.array(velx)
    vely=np.array(vely)
    return velx, vely, delta_t

def compute_VAC(velx,vely):
    totalsize = len(velx)
    vac=[]
    for tau in range(totalsize):
        tau=tau+1
        deltaCoords = np.multiply(velx[0+tau:totalsize],velx[0:totalsize-tau]) + np.multiply(vely[0+tau:totalsize],vely[0:totalsize-tau])
        val = np.sum(deltaCoords) # dx^2+dy^2+dz^2
        vac.append(np.mean(val)) # average
       # print velx[1+tau:totalsize]
       # print ' '
       # print deltaCoords
       # print ' '
        #print tau
    #...
    vac = np.array(vac)
    return vac

def insertIntoDict(key,val,aDict):
    if not key in aDict:
        aDict[key] = [val]
    else:
        aDict[key].append(val)
# http://stackoverflow.com/questions/18817789/how-to-add-values-to-existing-dictionary-key-python


#
#
#
#
#


# directory to save analysis
saveDir = "/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/"


# =1 if you want to display ID, directory
verbose = 0
# =1 if plotting
doPlot = 0

##base directory for data
#baseDir = "/Users/sagarsetru/Documents/Princeton/wingreen/for NED/"
#
#subBaseDir_a1 = "FloA/1/"
#subBaseDir_a2 = "FloA/2/"
#subBaseDir_a3 = "FloA/3/"
#subBaseDir_t1 = "FloT/1/"
#subBaseDir_t2 = "FloT/2/"
#subBaseDir_t3 = "FloT/3/"
#subBaseDirs = ["FloA/1/", "FloA/2/", "FloA/3/", "FloT/1/", "FloT/2/", "FloT/3/"]
##subBaseDirs = ["FloA/1/", "FloA/2/",]
#
#fileNames_tracks = ["01FloA-tracks.xls",
#"02FloA-tracks.xls",
#"03FloA-tracks.xls",
#"01FloT-tracks.xls",
#"02FloT-tracks.xls",
#"03FloT-tracks.xls"]
#    
#fileNames_spots = ["01FloA-spots.xls",
#"02FloA-spots.xls",
#"03FloA-spots.xls",
#"01FloT-spots.xls",
#"02FloT-spots.xls",
#"03FloT-spots.xls"]

#fnameSpots = '/Volumes/Extreme SSD/lipidDomainData/files_spots.txt'
#fnameTracks = '/Volumes/Extreme SSD/lipidDomainData/files_tracks.txt'
fnameSpots = '/Volumes/Extreme SSD/lipidDomainData/Sagar MSD/files_spots.txt'
fnameTracks = '/Volumes/Extreme SSD/lipidDomainData/Sagar MSD/files_tracks.txt'
#fname = '/Users/sagarsetru/Documents/Princeton/woods hole physio 2019/data analysis/data processing/files.txt'

#%% load files

# from https://stackoverflow.com/questions/3277503/how-to-read-a-file-line-by-line-into-a-list
with open(fnameSpots) as f:
    txtfilesSpots = f.readlines()
txtfilesSpots = [x.strip() for x in txtfilesSpots]
txtfilesSpots.sort()

with open(fnameTracks) as f:
    txtfilesTracks = f.readlines()
txtfilesTracks = [x.strip() for x in txtfilesTracks]
txtfilesTracks.sort()

#%% make dictionaries for tracks

#dict for corrupted tracks IDs, use subBaseDir as key
corruptedTracks = {}
#dict for track lengths, use (subBaseDir, trackID) as key
trackLengths = {}


#for spotFile, trackFile in zip(txtfilesSpots,txtfilesTracks):
#    print(spotFile)
#    print(trackFile)   

counter1 = -1
counterCorrupted = 0

# loop to fill dictionaries
#for subBaseDir in subBaseDirs:
for spotFile, trackFile in zip(txtfilesSpots,txtfilesTracks):

    counter1 += 1
    
#    if counter1 == 255:
#        break
##        continue
#    #...
#    
#    if counter1 == 256:
#        continue
#    #...
#    
#    if counter1 == 257:
#        continue
#    #...
#    
#    if counter1 == 258:
#        continue
#    #...
#    
#    if counter1 == 259:
#        continue
#    #...
    
    #get current file names
#    fileName_track = fileNames_tracks[counter1]
#    fileName_spot = fileNames_spots[counter1]
    
    # load track file
#    flo_tracks = loadDataFile( baseDir, subBaseDir, fileName_track )
    flo_tracks = pd.read_csv(trackFile, header=None)
    
    # load spot file
#    flo_spots = loadDataFile( baseDir, subBaseDir, fileName_spot )
    flo_spots = pd.read_csv(spotFile,header=None)
    
    # get paths of each file, make sure they are the same
    fPathTrack = os.path.split(trackFile)[0]
    fPathSpot = os.path.split(spotFile)[0]
    
    # get number of the movie within this replicate
    trackMovieNum = os.path.split(trackFile)[1][0]
    spotMovieNum = os.path.split(spotFile)[1][0]
    
    # if they are equal, then move on, otherwise there is an error
    if fPathTrack == fPathSpot and trackMovieNum == spotMovieNum:
        
        fPath = fPathTrack
        
        movieNum = trackMovieNum
        
        keyForDict = fPath+movieNum

    else:
        
        if fPathTrack != fPathSpot:
            
            print('Paths of spots and tracks files are not identical! Fix.')
            
            print('Path for track file:')
            print(fPathTrack)
            
            print('Path for spot file:')
            print(fPathSpot)
            
        #...
        
        if trackMovieNum != spotMovieNum:
            
            print('Movie numbers within this replicate are not the same! Fix.')
            
            print('Track movie number:')
            print(trackMovieNum)
            print('Path for track file:')
            print(fPathTrack)
            
            print('Spot movie number:')
            print(spotMovieNum)
            print('Path for spot file:')
            print(fPathSpot)
            
            break
        #...
        
    #continue
    #break
    #...
    
    # get track ID from spot file
#    ID_spots = flo_spots.loc[1:,3]
    ID_spots = flo_spots.loc[1:,2]

    # get ndarray of unique IDs
    ID_unique = ID_spots.unique()
    # get number of unique IDs
    n_ID = ID_unique.size
    
    # get x position from spot file
#    x_spots = flo_spots.loc[1:,5]
    x_spots = flo_spots.loc[1:,4]
    
    # get y position from spot file
#    y_spots = flo_spots.loc[1:,6]
    y_spots = flo_spots.loc[1:,5]
    
    # get time stamp from spot file
#    t_spots = flo_spots.loc[1:,8]
    t_spots = flo_spots.loc[1:,7]
    
    print('Working on: ')
    print(keyForDict)
    print('file '+str(counter1)+' out of '+str(len(txtfilesSpots)))
    print(' ')

    counter_IDs = -1
    # loop through unique IDs
    for ID in ID_unique:
        #print(counter_IDs)
        counter_IDs += 1
        if verbose == 1:
#            print(subBaseDir)
            print('Path:')
            print(fPath)
            print()
            print('Movie number:')
            print(movieNum)
            print()
            print('ID:')
            print(ID)
        #...
        
        ID_indices = ID_spots[ID_spots == ID].index.tolist()
        # get x, y, t values
        x_df = x_spots[ID_indices]
        y_df = y_spots[ID_indices]
        t_df = t_spots[ID_indices]
        # get x y positions as a list of tuples [ [x1,y1],[x2,y2],...]
        xy = []
        t=[]
        for IDind in ID_indices:
            x_val = float(x_df[IDind])
            y_val = float(y_df[IDind])
            t_val = float(t_df[IDind])
            xy.append([x_val,y_val])
            t.append([t_val])
        #...
        # convert to numpy array
        xy=np.array(xy)
        t=np.array(t)

        # check if there are repeat times
        repeatInds = np.setdiff1d(np.arange(len(t)), np.unique(t, return_index=True)[1])
        if not repeatInds.size:
            corrupted = 0
        else:
            corrupted = 1
            counterCorrupted += 1
            #add to dictionary of corrupted tracks
#            insertIntoDict(subBaseDir,ID,corruptedTracks)
            insertIntoDict(keyForDict,ID,corruptedTracks)
#            if not keyForDict in corruptedTracks:
#                corruptedTracks[keyForDict] = [len(t)]
#            else:
#                corruptedTracks[keyForDict].append(len(t))
            #print 'corrupted'
            #print ID
            continue
        
        # add length of track to dictionary
#        insertIntoDict((subBaseDir,ID),len(t),trackLengths)
#        if ('/Volumes/Extreme SSD/lipidDomainData/Sagar MSD/Flotillins/mAmp/FloA/A1', '1') == keyForDict:
#            break
#        #...
        
        insertIntoDict((keyForDict,ID),len(t),trackLengths)
#        print(len(t))
        if len(trackLengths[(keyForDict,ID)]) > 1:
            break
        #...
        #print('inserting')
#        if not keyforDict in trackLengths:
#            trackLengths[(keyForDict,ID)] = [len(t)]
#        else:
#            trackLengths[(keyForDict,ID)].append(len(t))
#        #...

    #...
#...
    
#%%

#get shortest and longest track lengths
trackMinMax = [(min(a), max(a)) for a in zip(*trackLengths.values())]

# get total number of tracks
nTracks_all = len(trackLengths.keys())

# initialize matrices for all msd, vac values
msd_all = np.zeros([nTracks_all,trackMinMax[0][1]-1])
vac_all = np.zeros([nTracks_all,trackMinMax[0][1]-2])


# initialize dictionary to count number of tracks per video
nTracksPerVid = {}
# initialize dictionary for max track length per video
trackMaxPerVid = {}

#%% fill dictionaries:
# -to count number of tracks per video
# -to get longest track length
#for subBaseDir in subBaseDirs:
counter2 = -1
for spotFile, trackFile in zip(txtfilesSpots,txtfilesTracks):
    
    counter2 += 1

    counter = 0
    
    # get paths of each file, make sure they are the same
    fPathTrack = os.path.split(trackFile)[0]
    fPathSpot = os.path.split(spotFile)[0]
    
    # get number of the movie within this replicate
    trackMovieNum = os.path.split(trackFile)[1][0]
    spotMovieNum = os.path.split(spotFile)[1][0]
    
    # if they are equal, then move on, otherwise there is an error
    if fPathTrack == fPathSpot and trackMovieNum == spotMovieNum:
        
        fPath = fPathTrack
        
        movieNum = trackMovieNum
        
        keyForDict = fPath+movieNum

    else:
        
        if fPathTrack != fPathSpot:
            
            print('Paths of spots and tracks files are not identical! Fix.')
            
            print('Path for track file:')
            print(fPathTrack)
            
            print('Path for spot file:')
            print(fPathSpot)
            
        #...
        
        if trackMovieNum != spotMovieNum:
            
            print('Movie numbers within this replicate are not the same! Fix.')
            
            print('Track movie number:')
            print(trackMovieNum)
            
            print('Spot movie number:')
            print(spotMovieNum)
            
        #...
        
        continue
    #...
    
#    if counter2 == 33:
#        break
#    #...
    
    print('Working on: ')
    print(keyForDict)
    print('file '+str(counter2)+' out of '+str(len(txtfilesSpots)))
    print(' ')
    
    maxTrack = [0]
    for key in trackLengths.keys():
        key1 = key[0]
        if keyForDict == key1:
            counter += 1
            keyUse = key
#            print(keyUse)
            track_L = trackLengths[key]
#            print(track_L)
#            if len(track_L) > 1:
#                break
#            #...
#            print(track_L)
            if track_L > maxTrack:
#                print(track_L)
#                print(maxTrack)
#                print(key)
#                print(track_L)
                maxTrack = [max(track_L)]
#                maxTrack = track_L
            #...
        #...
    #...
    
#    print(maxTrack)
#    if len(maxTrack) > 1:
#        break
    
#    insertIntoDict(subBaseDir,maxTrack,trackMaxPerVid)
#    insertIntoDict(subBaseDir,counter,nTracksPerVid)
    insertIntoDict(keyForDict,maxTrack,trackMaxPerVid)
    insertIntoDict(keyForDict,counter,nTracksPerVid)
    
#...



#%% loop to fill msd, vac matrices, and save
    
doSave = 1
counter1 = -1
counterCorrupted = 0

# loop to fill dictionaries
#for subBaseDir in subBaseDirs:
for spotFile, trackFile in zip(txtfilesSpots,txtfilesTracks):

    counter1 += 1
    #get current file names
    fileName_track = txtfilesTracks[counter1]
    fileName_spot = txtfilesSpots[counter1]

    # get number of the movie within this replicate
    trackMovieNum = os.path.split(trackFile)[1][0]
    spotMovieNum = os.path.split(spotFile)[1][0]
    
    # if they are equal, then move on, otherwise there is an error
    if fPathTrack == fPathSpot and trackMovieNum == spotMovieNum:
        
        fPath = fPathTrack
        
        subBaseDir = fPath[len(baseDirString):]
        
        movieNum = trackMovieNum
        
        keyForDict = fPath+movieNum

    else:
        
        if fPathTrack != fPathSpot:
            
            print('Paths of spots and tracks files are not identical! Fix.')
            
            print('Path for track file:')
            print(fPathTrack)
            
            print('Path for spot file:')
            print(fPathSpot)
            
        #...
        
        if trackMovieNum != spotMovieNum:
            
            print('Movie numbers within this replicate are not the same! Fix.')
            
            print('Track movie number:')
            print(trackMovieNum)
            
            print('Spot movie number:')
            print(spotMovieNum)
            
            #break
        
        #...
    #...

    # get longest track length for this vid
#    maxTrackLength = trackMaxPerVid[subBaseDir]
    maxTrackLength = trackMaxPerVid[keyForDict]
    
    # get the total number of tracsk for this vid
#    nTracks = nTracksPerVid[subBaseDir]
    nTracks = nTracksPerVid[keyForDict]

    # initialize msd and vac matrices for all tracks from this vid
#    msd_vid = np.zeros([nTracks[0],maxTrackLength[0][0]])
#    vac_vid = np.zeros([nTracks[0],maxTrackLength[0][0]])
#    msd_vid = np.zeros([nTracks[0],maxTrackLength[0][0]-1])
#    vac_vid = np.zeros([nTracks[0],maxTrackLength[0][0]-2])
    msd_vid = np.zeros([nTracks[0],205])
    vac_vid = np.zeros([nTracks[0],205])

    # initialize list for indices where msd and vac are zero
    msd_zero = []
    vac_zero = []

    # arrays for mean and sem msd and vac values per tau
    msd_mean = np.zeros([1,maxTrackLength[0][0]-1])
    msd_sem = np.zeros([1,maxTrackLength[0][0]-1])
    msd_std = np.zeros([1,maxTrackLength[0][0]-1])
    msd_btstrp = np.zeros([2,maxTrackLength[0][0]-1])
    
    vac_mean = np.zeros([1,maxTrackLength[0][0]-2])
    vac_sem = np.zeros([1,maxTrackLength[0][0]-2])  
    vac_std = np.zeros([1,maxTrackLength[0][0]-2])
    vac_btstrp = np.zeros([2,maxTrackLength[0][0]-2])
    
    # load track file
#    flo_tracks = loadDataFile( baseDir, subBaseDir, fileName_track )
    flo_tracks = pd.read_csv(trackFile, header=None)
    
    # load spot file
#    flo_spots = loadDataFile( baseDir, subBaseDir, fileName_spot )
    flo_spots = pd.read_csv(spotFile,header=None)
    
    
    # get track ID from spot file
#    ID_spots = flo_spots.loc[1:,3]
    ID_spots = flo_spots.loc[1:,2]

    # get ndarray of unique IDs
    ID_unique = ID_spots.unique()
    # get number of unique IDs
    n_ID = ID_unique.size
    
    # get x position from spot file
#    x_spots = flo_spots.loc[1:,5]
    x_spots = flo_spots.loc[1:,4]
    
    # get y position from spot file
#    y_spots = flo_spots.loc[1:,6]
    y_spots = flo_spots.loc[1:,5]
    
    # get time stamp from spot file
#    t_spots = flo_spots.loc[1:,8]
    t_spots = flo_spots.loc[1:,7]

    # get list of corrupted tracks for this video
#    cTracks = corruptedTracks[subBaseDir]
    cTracks = corruptedTracks[keyForDict]

    counter_IDs = -1
    counter_IDs_keep = -1
    # loop through unique IDs
    for ID in ID_unique:
        counter_IDs += 1
        if verbose == 1:
#            print(subBaseDir)
#            print(ID)
            print('Path:')
            print(fPath)
            print()
            print('Movie number:')
            print(movieNum)
            print()
            print('ID:')
            print(ID)
        #...
        
        if ID in cTracks:
            if verbose == 1:
                print('Skipping this track (corrupted)')
            #...
            
            continue
        
        #...
        
        counter_IDs_keep += 1
        
        ID_indices = ID_spots[ID_spots == ID].index.tolist()
        
        # get x, y, t values
        x_df = x_spots[ID_indices]
        y_df = y_spots[ID_indices]
        t_df = t_spots[ID_indices]
        
        # get x y positions as a list of tuples [ [x1,y1],[x2,y2],...]
        xy = []
        t = []
        
        for IDind in ID_indices:
            x_val = float(x_df[IDind])
            y_val = float(y_df[IDind])
            t_val = float(t_df[IDind])
            xy.append([x_val,y_val])
            t.append([t_val])
        #...
        # convert to numpy array
        xy=np.array(xy)
        t=np.array(t)

        # calculate the MSD
        msd = compute_MSD(xy)

        # get velocities, delta_t
        velx, vely, delta_t = compute_vel(xy,t)
        #print 'delta_t: ',delta_t

        # get velocity velocity autocorrelation
        vac = compute_VAC(velx,vely)
        #print 'msd: ',msd
        #print 'vac: ',vac
        #print 'vac[0:len(vac)-1]: ',vac[0:len(vac)-1]

        if 0 in msd:
            #msd_zero.append([counter_IDs_keep,int(np.where(msd==0)[0])])
            msd_zero.append(int(np.where(msd==0)[0]))
            #print '0 in msd'
            #print 'msd: ',msd
            #print subBaseDir, ID
        #...
        if 0 in vac[0:len(vac)-1]:
            #vac_zero.append([counter_IDs_keep,int(np.where(vac[0:len(vac)-1]==0)[0])])
            vac_zero.append(int(np.where(vac[0:len(vac)-1]==0)[0]))
            #print '0 in vac'
            #print 'vac: ',vac
            #print 'vac[0:len(vac)-1]: ',vac[0:len(vac)-1]
            #print subBaseDir, ID
        #..
        
        # load msd, vac into msd, vac matrix for this vid
        msd_vid[counter_IDs_keep,0:len(msd)] = msd
        vac_vid[counter_IDs_keep,0:len(vac)] = vac
#        vac_vid[counter_IDs_keep,0:len(vac)-1] = vac[0:len(vac)-1]
    #...
    
    #print 'msd_zero: ', msd_zero
    #print 'vac_zero: ',vac_zero
    #print subBaseDir
    #print 'delta_t: ',delta_t
    
    # loop through each tau, collect mean and sem for MSD
    for tau in range(maxTrackLength[0][0]-1):
        # get all msd for this tau
        msd_tau = msd_vid[:,tau]
        # get non zero indices
        ind_nonzero = np.nonzero(msd_tau)
        msd_tau_keep = msd_tau[ind_nonzero]
        #add zero value if an msd value in this tau range was zero
        if tau in msd_zero:
            dict_zero = collections.Counter(msd_zero)
            nZero = dict_zero[tau]
            #print 'msd, nZero: ',nZero,', tau: ',tau
            msd_tau_keep = np.append(msd_tau_keep,np.zeros([1,nZero]))
        #...
        # calculate and store mean, std, sem, bootstrap values for msd
        msd_mean[0,tau] = np.mean(msd_tau_keep)
        msd_std[0,tau] = np.std(msd_tau_keep)
        msd_sem[0,tau] = np.divide(np.std(msd_tau_keep),np.sqrt(len(msd_tau_keep)))
        
        btstrp = bs.bootstrap(msd_tau_keep,stat_func=bs_stats.mean)
        msd_btstrp[0,tau] = btstrp.lower_bound
        msd_btstrp[1,tau] = btstrp.upper_bound
    #...
        
    # loop through each tau, collect mean and sem for VAC        
    for tau in range(maxTrackLength[0][0]-2):
        # get all vac foro this tau
        vac_tau = vac_vid[:,tau]
        # get non zero indices
        ind_nonzero = np.nonzero(vac_tau)
        vac_tau_keep = vac_tau[ind_nonzero]
        #add zero values if any msd values in this tau range was zero
        if tau in vac_zero:
            dict_zero = collections.Counter(vac_zero)
            nZero = dict_zero[tau]
            #print 'vac, nZero: ',nZero,', tau: ',tau
            vac_tau_keep = np.append(vac_tau_keep,np.zeros([1,nZero]))
        #...
        # calculate and store mean and sem values for vac
        vac_mean[0,tau] = np.mean(vac_tau_keep)
        vac_std[0,tau] = np.std(msd_tau_keep)
        vac_sem[0,tau] = np.divide(np.std(vac_tau_keep),np.sqrt(len(vac_tau_keep)))
        
        btstrp = bs.bootstrap(vac_tau_keep,stat_func=bs_stats.mean)
        vac_btstrp[0,tau] = btstrp.lower_bound
        vac_btstrp[1,tau] = btstrp.upper_bound
    #...
    
    if doSave == 1:
        #move to analysis directory for saving data
        os.chdir(saveDir)

        #make directory to save files
        if not os.path.exists(subBaseDir):
            os.makedirs(subBaseDir)
        #...
        os.chdir(subBaseDir)
        #save MSD to csv file
        if printSave == 1:
            print('Saving msd_track'+ID+'.csv')
        #...
        np.savetxt('msd_track_mean'+ID+'.csv',msd_mean,delimiter=',')

        #save VAC to csv file
        if printSave == 1:
            print('Saving vac_track'+ID+'.csv')
        #...
        np.savetxt('vac_track_mean'+ID+'.csv',vac_mean,delimiter=',')
            
    # do plotting
    if doPlot:
        # set number of frames for short plots
        nFrames_short = 9;
        os.chdir(saveDir)

        #make directory to save files
        if not os.path.exists(subBaseDir):
            os.makedirs(subBaseDir)
        #...
        os.chdir(subBaseDir)

        # set x arrays for plotting
        short_x = np.multiply(np.array(range(nFrames_short))+1,delta_t)
        long_x = np.multiply(np.array(range(len(msd_mean[0])))+1,delta_t)

        #set vac arrays for plotting
        short_vac = vac_mean[0][0:nFrames_short]
        long_vac = vac_mean[0]
        short_vac_sem = vac_sem[0][0:nFrames_short]
        long_vac_sem = vac_sem[0]

        #set msd arrays for plotting
        short_msd = msd_mean[0][0:nFrames_short]
        long_msd = msd_mean[0]
        short_msd_sem = msd_sem[0][0:nFrames_short]
        long_msd_sem = msd_sem[0]
        
        #plot mean VAC LONG
        plt.figure(0)
        plt.errorbar(long_x[0:len(long_x)-1],long_vac,yerr=long_vac_sem)
        plt.xlabel(r'$\tau$',fontsize=18)
        #plt.ylabel(r'$\langle$VAC$\rangle\pm$SEM',fontsize=18)
        plt.ylabel('Avg. VAC +/- SEM',fontsize=18)
        plt.tick_params(axis='both', which='major', labelsize=15)
        plt.xlim([0,max(long_x)+delta_t])
        if 'FloA' in subBaseDir:
            plt.ylim([-0.2,0.2])
        else:
            plt.ylim([-0.06,0.06])
        plt.figure(0)
        plt.savefig('vac_mean_'+fileName_spot[0:6]+'_long.pdf')
        plt.clf()

        #plot mean VAC SHORT
        plt.figure(1)
        plt.errorbar(short_x,short_vac,yerr=short_vac_sem)
        plt.xlabel(r'$\tau$',fontsize=18)
        #plt.ylabel(r'$\langle$VAC$\rangle\pm$SEM',fontsize=18)
        plt.ylabel('Avg. VAC +/- SEM',fontsize=18)
        plt.tick_params(axis='both', which='major', labelsize=15)
        plt.xlim([0,max(short_x)+delta_t])
        #if 'FloA' in subBaseDir:
        plt.ylim([-0.06,0.06])
        #else:
        #plt.ylim([-0.06,0.06])
        plt.figure(1)
        plt.savefig('vac_mean_'+fileName_spot[0:6]+'_short.pdf')
        plt.clf()

        #plot mean MSD LONG
        plt.figure(2)
        plt.errorbar(long_x,long_msd,yerr=long_msd_sem)
        plt.xlabel(r'$\tau$',fontsize=18)
        plt.ylabel('Avg. MSD +/- SEM',fontsize=18)
        plt.tick_params(axis='both', which='major', labelsize=15)
        plt.figure(2)
        plt.xlim([0,max(long_x)+delta_t])
        #plt.ylim([0,.1])
        plt.savefig('msd_mean_'+fileName_spot[0:6]+'_long.pdf')
        plt.clf()
        
        #plot mean MSD SHORT
        plt.figure(3)
        plt.errorbar(short_x,short_msd,yerr=short_msd_sem)
        plt.xlabel(r'$\tau$',fontsize=18)
        plt.ylabel('Avg. MSD +/- SEM',fontsize=18)
        plt.tick_params(axis='both', which='major', labelsize=15)
        plt.figure(3)
        plt.ylim([0,.1])
        plt.xlim([0,max(short_x)+delta_t])
        plt.savefig('msd_mean_'+fileName_spot[0:6]+'_short.pdf')
        plt.clf()

        #plot mean MSD on log log with ref curves LONG
        xtest1 = np.multiply(range(len(msd_mean[0])),delta_t)
        xtest2 = xtest1
        xtest2[0]=0.1
        ytest1 = xtest1
        ytest1 = ytest1-ytest1[0]+short_msd[0]
        ytest2 = np.square(xtest2)
        ytest2 = ytest2-ytest2[0]+short_msd[0]
        plt.figure(4)
        #plt.errorbar(np.multiply(range(len(msd_mean[0])),delta_t),msd_mean[0],yerr=msd_sem[0],label='Avg. MSD')
        plt.errorbar(long_x,long_msd,yerr=long_msd_sem,label='Avg. MSD')
        plt.plot(xtest2,ytest2,label='Slope=2')
        plt.plot(xtest1,ytest1,label='Slope=1')
        plt.xlabel(r'$\tau$',fontsize=18)
        plt.ylabel('Avg. MSD +/- SEM',fontsize=18)
        plt.xscale('log')
        plt.yscale('log')
        plt.tick_params(axis='both', which='major', labelsize=15)
        plt.legend(loc=2)
        plt.xlim([min(long_x),max(long_x)])
        plt.ylim([0,1000])
        plt.figure(4)
        plt.savefig('msd_mean_'+fileName_spot[0:6]+'_loglog_long.pdf')
        plt.clf()

        #plot mean MSD on log log with ref curves SHORT
        plt.figure(5)
        #plt.errorbar(np.multiply(range(len(msd_mean[0])),delta_t),msd_mean[0],yerr=msd_sem[0],label='Avg. MSD')
        plt.errorbar(short_x,short_msd,yerr=short_msd_sem,label='Avg. MSD')
        plt.plot(xtest2,ytest2,label='Slope=2')
        plt.plot(xtest1,ytest1,label='Slope=1')
        plt.xlabel(r'$\tau$',fontsize=18)
        plt.ylabel('Avg. MSD +/- SEM',fontsize=18)
        plt.xscale('log')
        plt.yscale('log')
        plt.tick_params(axis='both', which='major', labelsize=15)
        plt.legend(loc=2)
        plt.ylim([0,1000])
        plt.xlim([min(short_x),max(short_x)])
        #plt.ylim([0,.1])
        plt.figure(5)
        plt.savefig('msd_mean_'+fileName_spot[0:6]+'_loglog_short.pdf')
        plt.clf()

# make figure of histogram of all track lengths

# convert values in dict to list of lists, then to list
all_trackLengths = list(trackLengths.values())
all_trackLengths = list(itertools.chain.from_iterable(all_trackLengths))

plt.figure(2)
#plt.hist(zip(*all_trackLengths),bins=range(trackMinMax[0][1]))
#plt.hist(list(all_trackLengths),bins=range(trackMinMax[0][1]))
plt.hist(all_trackLengths,bins=range(trackMinMax[0][1]))
#plt.bar(list(all_trackLengths),bins=range(trackMinMax[0][1]))
plt.ylabel('Frequency',fontsize=18)
plt.xlabel('No. frames per track',fontsize=18)
#plt.tick_params(axis='both', which='major', labelsize=15)
os.chdir(saveDir)
plt.savefig('hist_allTrackLengths.pdf')
plt.clf()
