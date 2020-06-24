# Sagar Setru
# 2016 07 14

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import itertools


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
#%%
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
       # print tau
    #...
    vac = np.array(vac)
    return vac

# directory to save analysis
saveDir = "/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/"

# =1 if saving msd, vac data
doSave = 1
# =1 to display files that get saved
printFilesSaved = 0
# =1 if you want to display ID, directory
verbose = 0
# =1 if plotting
doPlot = 0

#=1 if setting a maximum tau
setMax = 0
#set maximum number of time steps to plot, will only use if setMax == 1
maxTau = 10

#base directory for data
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


#fnameSpots = '/Volumes/Extreme SSD/lipidDomainData/Sagar MSD/files_spots.txt'
#fnameTracks = '/Volumes/Extreme SSD/lipidDomainData/Sagar MSD/files_tracks.txt'
#
#fnameSpots = '/Volumes/Extreme SSD/lipidDomainData/Sagar MSD/Flotillins/files_fig13e_spots.txt'
#fnameTracks = '/Volumes/Extreme SSD/lipidDomainData/Sagar MSD/Flotillins/files_fig13e_tracks.txt'

fnameTracks = '/Volumes/Extreme SSD/lipidDomainData/Sagar MSD/Flotillins/files_fig13e_ctrl_tracks.txt'
fnameSpots = '/Volumes/Extreme SSD/lipidDomainData/Sagar MSD/Flotillins/files_fig13e_ctrl_spots.txt'

#fnameSpots = '/Volumes/Extreme SSD/lipidDomainData/Sagar MSD/files_spots_ctrlVal.txt'
#fnameTracks = '/Volumes/Extreme SSD/lipidDomainData/Sagar MSD/files_tracks_ctrlVal.txt'

#fnameSpots = '/Volumes/Extreme SSD/lipidDomainData/Sagar MSD/files_spots_ugt.txt'
#fnameTracks = '/Volumes/Extreme SSD/lipidDomainData/Sagar MSD/files_tracks_ugt.txt'

#fnameSpots = '/Volumes/Extreme SSD/lipidDomainData/Sagar MSD/files2_spots.txt'
#fnameTracks = '/Volumes/Extreme SSD/lipidDomainData/Sagar MSD/files2_tracks.txt'
#fnameSpots = '/Volumes/Extreme SSD/lipidDomainData/Sagar MSD/files3_spots.txt'
#fnameTracks = '/Volumes/Extreme SSD/lipidDomainData/Sagar MSD/files3_tracks.txt'
baseDirString = '/Volumes/Extreme SSD/lipidDomainData/Sagar MSD/'

# from https://stackoverflow.com/questions/3277503/how-to-read-a-file-line-by-line-into-a-list
with open(fnameSpots) as f:
    txtfilesSpots = f.readlines()
txtfilesSpots = [x.strip() for x in txtfilesSpots]
txtfilesSpots.sort()

with open(fnameTracks) as f:
    txtfilesTracks = f.readlines()
txtfilesTracks = [x.strip() for x in txtfilesTracks]
txtfilesTracks.sort()

#%%
counter1 = 0
counterCorrupted = 0

# loop to fill dictionaries
#for subBaseDir in subBaseDirs:
for spotFile, trackFile in zip(txtfilesSpots,txtfilesTracks):

    counter1 += 1
    
    
#    if counter1 == 255:
#        continue
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
#    print(fPathTrack)
#    print(fPathSpot)
    
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
            
        #...
        
        continue
    #...
    
    print('Working on: ')
    print(keyForDict)
    print('file '+str(counter1)+' out of '+str(len(txtfilesSpots)))
    print(' ')
    
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

    counter_IDs = -1
    # loop through unique IDs
    for ID in ID_unique:
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

        # skip if only going up to a certain length
        if setMax == 1:
            if len(t) > maxTau:
#                print('skipping')
                continue
            #...
        #...
        
        # skip if the track is just 1 time point long
        if len(xy) <= 1:
            continue
        #...

        # check if there are repeat times
        repeatInds = np.setdiff1d(np.arange(len(t)), np.unique(t, return_index=True)[1])
        if not repeatInds.size:
            corrupted = 0
        else:
            corrupted = 1
            counterCorrupted += 1
            #print 'corrupted'
            #print ID
#            print('skipping')
            continue
            
        # calculate the MSD
        msd = compute_MSD(xy)

        # get velocities, delta_t
        velx, vely, delta_t = compute_vel(xy,t)
        
#        if counter1 == 1:
#            break
        #...

        # get velocity velocity autocorrelation
        vac = compute_VAC(velx,vely)
        
        # get tau values for msd
        tau_msd = np.multiply(np.array(range(len(msd)))+1,delta_t)
        
        # get tau values for vac
        tau_vac = np.multiply(np.array(range(len(vac)))+1,delta_t)
        
        msd_and_tau = np.vstack((tau_msd,msd)).T
        
        vac_and_tau = np.vstack((tau_vac,vac)).T

        #plt.plot(msd)
        #plt.plot(vac)
        #print vac
        if doSave == 1:
            #move to analysis directory for saving data
            os.chdir(saveDir)

            #make directory to save files
            if not os.path.exists(subBaseDir):
                os.makedirs(subBaseDir)
            #...
            os.chdir(subBaseDir)
            #save MSD to csv file
            if printFilesSaved == 1:
                print('Saving msd_track'+ID+'.csv')
            #...
            np.savetxt('msd_track'+ID+'.csv',msd_and_tau,delimiter=',')

            #save VAC to csv file
            if printFilesSaved == 1:
                print('Saving vac_track'+ID+'.csv')
            #...
            np.savetxt('vac_track'+ID+'.csv',vac_and_tau,delimiter=',')
        #...

        # do plotting
        if doPlot == 1:
            #move to analysis directory for saving plots
            os.chdir(saveDir)

            #make directory to save files
            if not os.path.exists(subBaseDir):
                os.makedirs(subBaseDir)
            #...
            os.chdir(subBaseDir)
            
            plt.figure(0)
            plt.plot(np.multiply(np.array(range(len(vac)))+1,delta_t),vac)
            plt.xlabel(r'$\tau$',fontsize=18)
            plt.ylabel('VAC',fontsize=18)
            plt.tick_params(axis='both', which='major', labelsize=15)


            plt.figure(1)
            plt.plot(np.multiply(np.array(range(len(msd)))+1,delta_t),msd)
            plt.xlabel(r'$\tau$',fontsize=18)
            plt.ylabel('MSD',fontsize=18)
            plt.tick_params(axis='both', which='major', labelsize=15)
            
    if doPlot == 1:
        if setMax == 1:
            plt.figure(0)
            plt.savefig('vac_'+str(maxTau)+'tau'+fileName_spot[0:6]+'.pdf')
            plt.clf()

            plt.figure(1)
            plt.savefig('msd_'+str(maxTau)+'tau'+fileName_spot[0:6]+'.pdf')
            plt.clf()
        else:
            plt.figure(0)
            plt.savefig('vac_all'+fileName_spot[0:6]+'.pdf')
            plt.clf()

            plt.figure(1)
            plt.savefig('msd_all'+fileName_spot[0:6]+'.pdf')
            plt.clf()

            #with PdfPages('vac_curves.pdf') as pdf:
            #    plt.figure(0)
            #    pdf.savefig()
            #    plt.close()
                
            #with PdfPages('msd_curves.pdf') as pdf:
            #    plt.figure(1)
            #    pdf.savefig()
            #    plt.close()       
                #plt.style.use('ggplot')

