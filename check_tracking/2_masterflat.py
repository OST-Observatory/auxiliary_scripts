#! /usr/bin/python
 
import os, sys, numpy as np
from astropy.io import fits

# 1. read path
sourceDir 	= input('\nEnter directory path to your flatfields:\t')
color		= input('Enter the filter you used for the flatfields:\t')
darkPath	= input('Enter path to your the masterdark:\t\t')
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
	
# 2. open the master darkframe and read the first hdu (header, data)	
darkFile = fits.open('%s'%(darkPath))
darkData = darkFile[0].data	

# 3. open the first flatfield file and read the first hdu (header, data) and subtract the master darkframe
mainFile  = fits.open('%s/%s'%(sourceDir,fileList[1]))
flatData  = mainFile[0].data - darkData
flatHeader= mainFile[0].header

# 4. open the remaining flatfield files subtract the master darkframe and add pixel-values of the image arrays to the first one
for i in range(1,len(fileList)):
	temp = fits.open('%s/%s'%(sourceDir,fileList[i]))
	tempData=temp[0].data - darkData
	flatData=flatData + tempData
	
# 5. get exposure time from the fits-header, round it & convert it to string
exposureTime=flatHeader['EXPTIME']
exposureTime=str(np.around(exposureTime, decimals=1))

# 6. update value for number of added files in the fits-header
flatHeader['SNAPSHOT']='%s'%numberOfFiles

# 7. dividing added data by number of files -> mean value
flatData=flatData/numberOfFiles

# 8. round data to integers
flatData=np.around(flatData, decimals=0)

# 9. insert header and data into hdu-unit and scaling the image
hdu=fits.PrimaryHDU(flatData,flatHeader)
# hdu.scale('int16', '', bzero=32768)

# 10. insert hdu-unit into hdu-list (images only contain 1 hdu -> list only contains 1 element)
hduList=fits.HDUList([hdu])

# 11. write hdu-list to file (remove it in case it already exists)
if os.path.exists('masterflat_%s.fits'%(color))==True:
	os.system('rm masterflat_%s.fits'%(color))

hduList.writeto('masterflat_%s.fits'%(color))
print('done: wrote masterflat_%s.fits\n'%(color))
