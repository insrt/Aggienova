import matplotlib.pyplot as plt
from matplotlib.widgets import CheckButtons
import numpy as np
import math
import pandas as pd

#'ASASSN-15lh',
SNlist=['ASASSN-15lh','SN2011by','SN2007af', 'SN2005ke', 'SN2005hk', 'SN2012dn','SN2007Y',  'SN2006jc',   'SN2006aj',  'SN2011dh', 'SN2012aw', 'SN2008aw',   'SN2010jl', 'SN2008es', 'SN2018hna']
SNlegend=['ASASSN-15lh - SLSN I?','SN2011by - Ia-blue','SN2007af - Ia-red', 'SN2005ke - Ia-91bg', 'SN2005hk - Ia-02cx','SN2012dn - Ia SC','SN2007Y - Ib',  'SN2006jc - Ibn',  'SN2006aj - Ic/GRB',  'SN2011dh - IIb', 'SN2012aw - IIP', 'SN2008aw - IIL',   'SN2010jl - IIn', 'SN2008es - SLSNII', "SN2018hna"]
SNexplosiondates=[57178.5-40, 55672.0,54157.12, 53685.77, 53685.1-15.0,56114,54163.3-20.0,  54018,  53784.1,  55712.5, 56002.0, 54526.5,   55470, 54590, 58411]
SNexplosionrefs=['40 days before peak Dong_etal_2016','Nugent_etal_2011','Brown_etal_2012a','Brown_etal_2012a', 'Phillips_etal_2007','Brown et al. 2013 peak -20','Stritzinger_etal_2009 B peak - 20', 'Nakana_etal_2006 a few days before discovery',   'Campana_etal_2006 GRB trigger time',  'Arcavi_etal_2011', 'Fraser_etal_2012', 'discover - 1.5 days',   'ten days before first observation Stoll et al. 2011', 'ten days before first observation']

#grab vel and dm constants for each supernova
filename = 'SwiftSNweblist.csv'
data = open(filename, 'r')
df = pd.read_csv(filename)

#dictionary to store vel_corr and dm values for each supernova
snconst = {}
for j in SNlist:
    sn = df.loc[df['SNname'] == j]
    vals = sn[sn.columns[22:41]].values.tolist()[0]
    col_names = sn.columns[22:41].values.tolist()
  
    for i in range(len(vals)):
        if(isinstance(vals[i],float) or vals[i] == " "):
            vals[i] = np.nan
        elif(not("ref" in col_names[i])):
            vals[i] = float(vals[i])
    snconst[j] = vals


#set up plot
plt.ion() #turns on interactive plotting
fig = plt.figure()
ax = fig.add_subplot(111)


graphs = []
sndata = {} #dictionary to store uvm2mjd, uvm2mag, uvm2magerr for later use
for i in range(len(SNlist)):
	snname = SNlist[i]
	filename =  snname + '_uvotB15.1.dat'
	savename = '_pylightcurve.jpg'
	data = open('data/' + filename, 'r')

	explosion_mjd = SNexplosiondates[i]
	#find hubble flow & distance modulus
	host_vel_array = snconst[snname][2:4]
	h0 = 72.0
	h0err = 5.0
	
	distance_mod_cor = 5*math.log(float(host_vel_array[0])/h0,10)+25 #hubble flow
	distance_mod_cor_err = math.sqrt(((5*float(host_vel_array[1]))/(float(host_vel_array[0])*math.log(10,10)))**2+((5*200)/(float(host_vel_array[0])*math.log(10,10)))**2 + ((5*5.0)/(h0*math.log(10,10)))**2)
	
	if(np.isnan(snconst[snname][13]) != True): 
		distance_mod_cor = float(snconst[snname][13])
		distance_mod_cor_err = float(snconst[snname][14])
	if(np.isnan(snconst[snname][10]) != True):
		distance_mod_cor = float(snconst[snname][10])
		distance_mod_cor_err = snconst[snname][9]
	if(np.isnan(snconst[snname][7]) != True):
		distance_mod_cor = float(snconst[snname][7])
		distance_mod_cor_err = float(snconst[snname][8])
	if(np.isnan(snconst[snname][4]) != True):
		distance_mod_cor = float(snconst[snname][4])
		distance_mod_cor_err = float(snconst[snname][5])

	if(snname == "SN2018hna"):
		distance_mod_cor = 30.52
		distance_mod_cor_err = 0.29
	# print("SN: %s : %0.2f" % (snname, distance_mod_cor))
	# print(distance_mod_cor_err)
	# print("\n")
	#Initializing lists needed to plot the different filters separately
	uvm2mjd = []
	uvm2mag = []
	uvm2magerr = []
	# bmjd = []
	# bmag = []
	# bmagerr = []
	# vmjd = []
	# vmag = []
	# vmagerr = []

	#Reading the data in from the file
	for line in data:
		if not line[0] == "#":
			continue
		lines = np.genfromtxt(data, dtype=[('filter','S20'),('mjd',float),('mag',None),('magerr',None)], usecols = (0,1,2,3), unpack=True)
	filters1 = lines['filter']
	mjd1 = lines['mjd']
	mag1 = lines['mag']
	magerr1 = lines['magerr']

	#I needed to get rid of the NULL values in the ...15.1.dat files, so the next several lines are to make sure the 
	#program doesn't shut down because of them
	filterslist = []
	mjdlist  = []
	maglist = []
	magerrlist = []
	for i in range(len(filters1)):
		if not np.isnan(mag1[i]):
			filterslist.append(filters1[i])
			maglist.append(mag1[i])
			mjdlist.append(mjd1[i])
			magerrlist.append(magerr1[i])
	filters = np.array(filterslist)
	mjd = np.array(mjdlist)
	mag = np.array(maglist)
	magerr = np.array(magerrlist)	



	#breaking up the filters, mjd, mag, and magerr arrays into separate arrays to make plotting easier
	for i in range(len(filters)):
		if filters[i] == 'UVM2':
			uvm2mjd.append(mjd[i]-explosion_mjd)
			uvm2mag.append(mag[i]-distance_mod_cor)
			uvm2magerr.append(magerr[i])


	#store the sn data for later use
	sndata[snname] = [uvm2mjd, uvm2mag, uvm2magerr, distance_mod_cor]
	
	if uvm2mjd > 0:
		line = ax.plot(uvm2mjd, uvm2mag, linestyle='-', linewidth=1, marker="*", label=snname)
		graphs.append(line)
		
		# ax.errorbar(uvm2mjd, uvm2mag, yerr=uvm2magerr, elinewidth=2, capthick = 2, marker="*")

	##########This begins plotting portion##########

#initialize with values of first SN
max_mjd = max(sndata[SNlist[1]][0])
min_mjd = min(sndata[SNlist[1]][0])
max_mag = max(sndata[SNlist[1]][1])
min_mag = min(sndata[SNlist[1]][1])
for i in range(len(SNlist)):
#The next several lines are just to make sure the program doesn't shut down if one of the filters is missing 
	if len(sndata[SNlist[i]][0]) > 0:
		# plt.scatter(sndata[SNlist[i]][0], sndata[SNlist[i]][1], linewidth=1, marker="*",facecolors='none', s=50, label=SNlist[i])
		# plt.errorbar(sndata[SNlist[i]][0], sndata[SNlist[i]][1], yerr=sndata[SNlist[i]][2], elinewidth=2, capthick = 2)

		if(max(sndata[SNlist[i]][0]) > max_mjd):
			max_mjd = max(sndata[SNlist[i]][0])
		if(min(sndata[SNlist[i]][0]) < min_mjd):
			min_mjd = min(sndata[SNlist[i]][0])
		if(max(sndata[SNlist[i]][1]) > max_mag):
			max_mag = max(sndata[SNlist[i]][1])
		if(min(sndata[SNlist[i]][1]) < min_mag):
			min_mag = min(sndata[SNlist[i]][1])
	


#If you have labels specified in the plotting commands above, all you have to do is add this line to make a legend
#Also, apparently you can set the location to 'best' so it'll place it wherever python thinks it's the best
ax.set_xlabel('Days Since Explosion')

ax.set_ylabel('Absolute Magnitude (uvm2)')

# ax.set_title(snname, fontsize = 16)
#plt.title('UVOT Light Curves', fontsize = 12) 
#there aren't title and subtitle commands so instead I used suptitle and title to get a makeshift main title and subtitle

# ax.axis([int(mjd.min()-explosion_mjd)-1, math.ceil(np.amax(mjd-explosion_mjd))+1,math.ceil(np.amax(mag-distance_mod_cor)+np.amax(magerr)), int(mag.min()-distance_mod_cor)-0.25])
ax.axis([int(min_mjd)-1, math.ceil(max_mjd)+1,math.ceil(max_mag)+1, int(min_mag)-0.25])

'''
In order to make the axis format correctly, I used the minimum and maxiumum values for the MJD +/- 1 (in order to give the plots some extra
space) on the x-axis, and the min and max values for the magnitude. The functions int() and math.ceil() are used to round the numbers to the
nearest whole number to make the ends of the graph neater.
'''
plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
#The line above makes sure that the x axis is formatted correctly. The Offset option in matplotlib is automatically activated, so for big
#numbers like the MJD, it cuts it off in a strange way that is awful to look at. This makes it readable.

leg = ax.legend(loc='best')
lined = dict()
for legline, origline in zip(leg.get_lines(),graphs):
	legline.set_picker(5)
	lined[legline] = origline

def onpick(event):
	legline = event.artist
	origline = lined[legline][0]
	print(origline)
	vis = not origline.get_visible()
	origline.set_visible(vis)
	if vis:
		legline.set_alpha(1.0)
	else:
		legline.set_alpha(0.2)
	fig.canvas.draw()
fig.canvas.mpl_connect('pick_event', onpick)


plt.show(block=False)

input("Press <ENTER> to continue\n")
fig.savefig(savename) #you can uncomment this line if you want to save the figure to the SN file
fig.show()