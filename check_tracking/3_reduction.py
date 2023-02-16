#! /usr/bin/python

import os, sys, numpy as np
from astropy.io import fits

# 1. read path
sourceDir 	= input('\nEnter directory path to your images:\t')
color		= input('Enter the filter you used for the images:\t')
darkPath	= input('Enter path to your the masterdark:\t\t')
flatPath	= input('Enter path to your the masterflat:\t\t')

if os.path.exists('%s'%(flatPath))==False:
	print("Sorry, path to masterflat not found. Exiting \n")
	sys.exit(0)
	
if os.path.exists('%s'%(darkPath))==False:
	print("Sorry, path to masterdark not found. Exiting \n")
	sys.exit(0)
	
if os.path.exists('%s/'%(sourceDir))==False:
	print("Sorry, path to directory not found. Exiting \n")
	sys.exit(0)

else:
	fileList = os.listdir(sourceDir)
	fileList.sort()
	
	# remove not fits entries
	tempList = []
	for i in range(0, len(fileList)):
		if(fileList[i].find(".FIT")!=-1 or fileList[i].find(".fit")!=-1 or fileList[i].find(".FITS")!=-1 or fileList[i].find(".fits")!=-1):
			tempList.append("%s"%(fileList[i]))
	fileList=tempList			
	#DEBUG: print fileList
	numberOfFiles=int(len(fileList))
	

# 2. create directory for the reduced pictures
if os.path.exists('reduced')==False:
	os.system('mkdir reduced')
if os.path.exists('reduced/%s'%color)==True:
	os.system('rm -R reduced/%s'%color)
	os.system('mkdir reduced/%s'%color)
else:
	os.system('mkdir reduced/%s'%color)


# 3. open the master darkframe and read the first hdu (header, data)	
darkFile = fits.open('%s'%(darkPath))
darkData = darkFile[0].data


# 4. open the master darkframe and read the first hdu (header, data)	
flatFile = fits.open('%s'%(flatPath))
flatData = flatFile[0].data
flatMedian = np.median(flatData)
flatData = flatData/flatMedian


# 5. open the image files, substract dark frame, divide by normalized flatfield and writ files to reduced directory
for i in range(0,len(fileList)):
	temp = fits.open('%s/%s'%(sourceDir,fileList[i]))
	tempHeader=temp[0].header
	tempData=temp[0].data - darkData
	tempData=tempData / flatData
	hdu=fits.PrimaryHDU(tempData,tempHeader)
	# scaling of image, has to be done for convert to work
	# hdu.scale('int16', '', bzero=32768)
	hduList=fits.HDUList([hdu])
	hduList.writeto('reduced/%s/%s'%(color,fileList[i]))


print('done: wrote files to directory reduced/%s\n'%(color))
