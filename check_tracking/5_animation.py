#! /usr/bin/python

import os

sourceDir 	= input('\nEnter directory path to images to animate:\t')

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

for i in range(0,numberOfFiles):
	os.system('convert %s/%s -normalize -resize 765x510 %s/%s.gif'%(sourceDir,fileList[i],sourceDir,fileList[i]))
	print(('done converting file # %s'%(i)))

os.system('mkdir %s/animation'%sourceDir)
os.system('mv %s/*.gif %s/animation'%(sourceDir,sourceDir))
os.system('convert %s/animation/*.gif -delay 100 -loop 0 %s/animation.gif'%(sourceDir, sourceDir))
os.system('rm -r %s/animation'%sourceDir)