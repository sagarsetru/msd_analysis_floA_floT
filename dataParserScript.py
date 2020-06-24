# Sagar Setru
# 2016 07 14

<<<<<<< HEAD

=======
>>>>>>> 9685d8720f079278ea5aad7fad293fe80babb4bb
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import itertools


def loadDataFile( baseDir, subBaseDir, fileName):
    "helper function to load bed files"
    return pd.read_csv(baseDir+subBaseDir+'/'+fileName, sep='\t', header=None)
#...

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

# directory to save analysis
saveDir = "/Users/sagarsetru/Documents/Princeton/wingreen/for NED/analysis/"

# =1 if saving msd, vac data
doSave = 0
# =1 if printing files that get saved
printSave = 0
# =1 if you want to display ID, directory
verbose = 0
# =1 if plotting
doPlot = 1

#=1 if setting a maximum tau
setMax = 0
#set maximum number of time steps to plot
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

fnameSpots = '/Volumes/Extreme SSD/lipidDomainData/files_spots.txt'
fnameTracks = '/Volumes/Extreme SSD/lipidDomainData/files_tracks.txt'

# from https://stackoverflow.com/questions/3277503/how-to-read-a-file-line-by-line-into-a-list
with open(fnameSpots) as f:
    txtfilesSpots = f.readlines()
txtfilesSpots = [x.strip() for x in txtfilesSpots]

with open(fnameTracks) as f:
    txtfilesTracks = f.readlines()
txtfilesTracks = [x.strip() for x in txtfilesTracks]

counter1 = -1
counterCorrupted = 0
for subBaseDir in subBaseDirs:
    counter1 += 1
    #get current file names
    fileName_track = fileNames_tracks[counter1]
    fileName_spot = fileNames_spots[counter1]
    
    # load track file
    flo_tracks = loadDataFile( baseDir, subBaseDir, fileName_track )
    # load spot file
    flo_spots = loadDataFile( baseDir, subBaseDir, fileName_spot )
    
    # get track ID from spot file
    ID_spots = flo_spots.loc[1:,3]
    # get ndarray of unique IDs
    ID_unique = ID_spots.unique()
    # get number of unique IDs
    n_ID = ID_unique.size
    # get x position from spot file
    x_spots = flo_spots.loc[1:,5]
    # get y position from spot file
    y_spots = flo_spots.loc[1:,6]
    # get time stamp from spot file
    t_spots = flo_spots.loc[1:,8]

    counter_IDs = -1
    # loop through unique IDs
    for ID in ID_unique:
        counter_IDs += 1
        if verbose == 1:
            print(subBaseDir)
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

        if setMax == 1:
            if len(t) > maxTau:
                continue
            #...
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
            continue
            
        # calculate the MSD
        msd = compute_MSD(xy)

        # get velocities, delta_t
        velx, vely, delta_t = compute_vel(xy,t)

        # get velocity velocity autocorrelation
        vac = compute_VAC(velx,vely)

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
            if printSave == 1:
                print('Saving msd_track'+ID+'.csv')
            #...
            np.savetxt('msd_track'+ID+'.csv',msd,delimiter=',')

            #save VAC to csv file
            if printSave == 1:
                print('Saving vac_track'+ID+'.csv')
            #...
            np.savetxt('vac_track'+ID+'.csv',vac,delimiter=',')
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

