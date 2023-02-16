#! /usr/bin/python

import os, sys, numpy as np
from astropy.io import fits

# this script will stack all presented fits files on top of each other
# for the purpose of investigating the tracking error we encounter with the OST :-/

targetDir = 'stacked'

# 1. get fileList
sourceDir 	= input('\nEnter directory path to your images:\t')

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

# 2. create output dir

if os.path.exists(targetDir)==False:
	os.system('mkdir %s'%(targetDir))
else:
	os.system('rm -r %s'%(targetDir))
	os.system('mkdir %s'%(targetDir))

# 3. copy first file
os.system('cp %s/%s %s'%(sourceDir,fileList[0],targetDir))

# 4. open first file
main	= fits.open('%s/%s'%(targetDir,fileList[0]))
mainData= main[0].data

# 4.1. get image dimensions
xDim=mainData.shape[1]
yDim=mainData.shape[0]

# 5. loop: for all i stack all images from 1 to i
for i in range(1,len(fileList)):
	
	temp 		= fits.open('%s/%s'%(sourceDir,fileList[i]))
	tempData	= temp[0].data
	tempHeader	= temp[0].header
	
	# 5.1 add files, update # of snapshots
	tempHeader['SNAPSHOT']='%s'%(i+1)
	mainData = mainData + tempData  
	tempData = mainData
	for x in range(0,xDim-1):
		for y in range(0,yDim-1):
			if tempData[y,x]>65535.:
				tempData[y,x]=65535.
	
	# 5.2 create new header (+ scale the image)
	hdu=fits.PrimaryHDU(tempData,tempHeader)
	
	# 5.3 insert hdu-unit into hdu-list (images only contain 1 hdu -> list only contains 1 element)
	hduList=fits.HDUList([hdu])
	
	# 5.4 write hdu-list to file 
	
	hduList.writeto('%s/%s'%(targetDir,fileList[i]))
	
	# 5.5 fix scaling problem :-/
	temp2		= fits.open('%s/%s'%(targetDir,fileList[i]))
	temp2Data	= temp2[0].data
	temp2Header	= temp2[0].header
	hdu2 = fits.PrimaryHDU(temp2Data,temp2Header)
	#
	# something is broken here:
	#
	# hdu2.scale('int16', '', bzero=32768)
	#hdu2.scale('int16', 'minmax')
	hduList2=fits.HDUList([hdu2])
	os.system('rm %s/%s'%(targetDir,fileList[i]))
	hduList2.writeto('%s/%s'%(targetDir,fileList[i]))
	#
	print(('done stacking file # %s'%(i)))
	
# 6. done
