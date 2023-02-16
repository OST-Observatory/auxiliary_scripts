#! /usr/bin/python
 
import os, sys, numpy as np
from astropy.io import fits

# 1. read path
sourceDir 	= input('\nEnter directory path to your dark frames:\t')
if os.path.exists('%s/'%(sourceDir))==False:
	print("Sorry, path not found. Exiting \n")
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
	
# 2. open the first darkframe file and read the first hdu (header, data)
mainFile=fits.open('%s/%s'%(sourceDir,fileList[1]))
darkData=mainFile[0].data
darkHeader=mainFile[0].header

# 3. open the remaining darkframe files and add pixel-values of the image arrays to the first one
for i in range(1,len(fileList)):
	temp=fits.open('%s/%s'%(sourceDir,fileList[i]))
	tempData=temp[0].data
	darkData=darkData+tempData

# 4. get exposure time from the fits-header, round it & convert it to string
exposureTime=darkHeader['EXPTIME']
exposureTime=str(np.around(exposureTime, decimals=1))

# 5. update value for number of added files in the fits-header
darkHeader['SNAPSHOT']='%s'%numberOfFiles

# 6. dividing added data by number of files -> mean value
darkData=darkData/numberOfFiles

# 7. round data to integers
darkData=np.around(darkData, decimals=0)

# 8. insert header and data into hdu-unit and scaling the image
hdu=fits.PrimaryHDU(darkData,darkHeader)
# hdu.scale('int16', '', bzero=32768)

# 9. insert hdu-unit into hdu-list (images only contain 1 hdu -> list only contains 1 element)
hduList=fits.HDUList([hdu])

# 10. write hdu-list to file (remove it in case it already exists)
if os.path.exists('masterdark_%ss.fits'%exposureTime)==True:
	os.system('rm masterdark_%ss.fits'%exposureTime)

hduList.writeto('masterdark_%ss.fits'%exposureTime)
print('done: wrote masterdark_%ss.fits\n'%exposureTime)


#######     some useful commands     ########
####### maybe write them to another file ####

#### opening a fits-file
# f = fits.open('fitsfile.fits')

#### accessing the image
# data=f[0].data

#### accessing header data
# header=f[0].header

#### access an attribute from header
# exposureTime = f[0].header['exptime']

#### printing value of of image pixel
# print data[1,1]

#### x-dimension of image
# data_poins_x=data.shape[1]
