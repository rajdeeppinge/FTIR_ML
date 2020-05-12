
# Author - Rajdeep Pinge
# Date - 25th May, 2017

'''
Code to extract IR spectrum data from spectrum images taken from Spectral DB for Organic Compounds

Link to DB: http://sdbs.db.aist.go.jp/sdbs/cgi-bin/direct_frame_top.cgi
'''

from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
import math
import os
from os import listdir, path
from os.path import isfile, join
from pandas import DataFrame


# directory where images generated from extracted data are stored
imgSavePath = "spectrumGraph/"


# source directory where all the data images have been stored
dataPath = "DataImages/"

# Extract data from all files in the directory. Modification done on 2nd June, 2017
onlyfiles = [f for f in listdir(dataPath) if isfile(join(dataPath, f))]

for imgFileName in onlyfiles:

	fileNameComp = imgFileName.split(".")

	if fileNameComp[-1] != "gif":		# is not an image file
		continue

	imgFileExt = "." + fileNameComp[-1]
	fileNameComp.remove(fileNameComp[-1])
	imgFileName = ".".join(fileNameComp)

	imgFile = Image.open(dataPath + imgFileName + imgFileExt, 'r')

	imgFile = imgFile.convert('L')

	cols, rows = imgFile.size
	#print rows, cols

	cropImg = imgFile.crop((30, 97, cols-1, 417))		#(30, 97, col, 417)

	cropPath = "CropImages/"
	cropFileName = "crop_" + imgFileName + ".png"

	# check if a directory exists to store cropped images
	if not os.path.exists(cropPath):
		os.makedirs(cropPath)		# if not, create one

	cropImg.save(cropPath + cropFileName)

	cropCols, cropRows = cropImg.size

	pix_val = list(cropImg.getdata())
	pix_val = np.array(pix_val).reshape((cropRows, cropCols))

	#print pix_val.shape

	reqPts = []

	for rowno in range(cropRows):
		for colno in range(cropCols):
			if pix_val[rowno, colno] == 0:
				reqPts.append((colno, ((cropRows-rowno)/float(cropRows)*100.0) ))

	reqPts = sorted(reqPts)

	pts = zip(*reqPts)		# splits x and y values from the tuple and makes their individual lists so as to make it easier to plot

	plt.figure(imgFileName + "_1")
	plt.plot(pts[0], pts[1], linewidth=1)
	#plt.axvline(x=max(pts[0])*2/5.2)
	plt.xlabel('Wave Number ($cm^{-1}$)')
	plt.ylabel('Transmittance (%T)')
	plt.title("Original Extracted FT-IR spectrum " + str(len(pts[0])) + " points")
	plt.grid(True)


	if not os.path.exists(imgSavePath):
				os.makedirs(imgSavePath)

	plt.savefig(imgSavePath + imgFileName + '_orig.png')


	#########################################################################################################
	############################## Truncating the extracted data ############################################
	'''
	Modification on 1st June, 2017

	'''

	x_axis = pts[0]
	y_axis = pts[1]

	count = 0
	sumTransmittance = 0
	currX = None
	prevX = None

	data = []

	for index in range(len(x_axis)):
		currX = x_axis[index]
		
		if (prevX is not None) and (currX != prevX):
			data.append( (prevX, sumTransmittance/float(count)) )
			sumTransmittance = 0
			count = 0

		sumTransmittance += y_axis[index]
		count += 1

		prevX = currX

	data = sorted(data)

	pts = zip(*data)

	splitPt = int(math.ceil(max(pts[0])*2/5.2))
	#print splitPt


	plt.figure(imgFileName + "_2")
	plt.plot(pts[0], pts[1], linewidth=1)
	#plt.axvline(x=splitPt
	plt.xlabel('Wave Number ($cm^{-1}$)')
	plt.ylabel('Transmittance (%T)')
	plt.title("Truncated FT-IR spectrum " + str(len(pts[0])) + " points")
	plt.grid(True)

	plt.savefig(imgSavePath + imgFileName + '_trunc.png')

	############################## Truncating the extracted data Ends #######################################


	#########################################################################################################
	############################### Normalization and interpolation #########################################
	'''
	Modification on 1st June, 2017

	'''

	x1_range = pts[0][:splitPt]
	x2_range = pts[0][splitPt:]

	y_range = pts[1]
	#y1_range = pts[1][:splitPt]
	#y2_range = pts[1][splitPt:]

	x1_range_diff = 4000 - 2000
	x2_range_diff = 2000 - 400

	x1_step = x1_range_diff / float(x1_range[-1] + 1)
	x2_step = x2_range_diff / float(x2_range[-1] - x2_range[0] + 1)

	val = 4000.0
	epsilon = 0.1
	x1_new_range = [val]
	y1_new_range = [y_range[0]]
	yiter = 1

	while val >= 2000.0 + x1_step + epsilon and yiter < len(x1_range):
		val -= x1_step
		correspYval = y_range[yiter]

		interpX = (val + x1_new_range[-1]) / 2.0
		interpY = (correspYval + y1_new_range[-1]) / 2.0

		x1_new_range.append(interpX)
		y1_new_range.append(interpY)

		x1_new_range.append(val)
		y1_new_range.append(correspYval)

		yiter += 1

#	print len(x1_new_range), len(y1_new_range)

	val = 2000.0
	x2_new_range = []
	y2_new_range = []

	while val >= 400.0 and yiter < len(y_range):
		x2_new_range.append(val)
		y2_new_range.append(y_range[yiter])
		val -= x2_step
		yiter += 1

#	print len(x2_new_range), len(y2_new_range)

	finalXRange = x1_new_range + x2_new_range
	finalYRange = y1_new_range + y2_new_range


	plt.figure(imgFileName + "_3")
	plt.plot(finalXRange, finalYRange, linewidth=1)
	#plt.axvline(x=2000.0)
	plt.xlim(finalXRange[0], finalXRange[-1])  # decreasing wave num
	plt.xlabel('Wave Number ($cm^{-1}$)')
	plt.ylabel('Transmittance (%T)')
	plt.title("Normalized FT-IR spectrum " + str(len(finalXRange)) + " points")
	plt.grid(True)

	plt.savefig(imgSavePath + imgFileName + '_normalized.png')

	############################## Normalization and Interpolation Ends #####################################



	#########################################################################################################
	####################################### Generate csv file ###############################################

	waveNumLst = finalXRange
	transmittance = finalYRange

	# Convert the wave nos and transmittance to a DataFrame i.e. a table
	df = DataFrame({'Wave Number (cm-1)': waveNumLst, 'Transmittance (%T)': transmittance})

	# adjust the columns of DataFrame
	dfCols = ['Wave Number (cm-1)', 'Transmittance (%T)']
	df = df[dfCols]

	# Directory Where csv data files are stored
	csvData = "csvData/"

	# check if a directory exists to store csv files
	if not os.path.exists(csvData):
		os.makedirs(csvData)		# if not, create one

	# Convert DatqFrame to make it compatible to csv with a "," separator
	df.to_csv(csvData + imgFileName + ".csv", sep=',', encoding='utf-8', index=False)

	################################## Generate csv file Ends ###############################################


	print "Image " + imgFileName + " extraction complete. " + str(len(waveNumLst)) + " data points available" 


######################## show the graphs

#plt.show()