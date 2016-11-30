# Here I am importing all of the necessary packages for this.
import numpy as np 
import netCDF4
import pandas as pd 
import xarray as xr 
from pylab import *
import cartopy.crs as ccrs
from matplotlib import pyplot as plt
import seaborn as sns
import os
from scipy.interpolate import griddata
from datetime import datetime
from datetime import timedelta
from dateutil import parser


class PostProcess:
	def __init__(self,date,number_of_ensembles,number_of_obs,obs_precision,time_periods,time_period_format):
		
		#	date = the date you want to run the Post-Processing on, YYYYMMDD
		#	number_of_ensembles = number of ensembles you're doing analysis on (e.g.,5)
		#	number_of_obs = number of observations you're doing analysis on
		#	time_periods_ens  = which indices are you going to be comparing against? (it's a list of indices)
		#	time_periods_obs = which indices are you going to be comparing vs ens against (a list of indices) ?
		# 	obs_precisions. 0.25 degree = 4, 0.5 degree = 2...& so on.

		#   time periods = list of string times, in hours, in which we are running the analysis.  e.g. ["00h","06h","12h"]
		#   time_period_format = a tuple that will convert the time periods into the hours.


		self.date = date
	#Using path and date, create Obs_List and Ens_List 
		self.dt = parser.parse(date)
		self.dt_folder = (self.dt).strftime('%b%d')
		self.path_to_ens = "Ensembles/" + self.dt_folder
		self.DS=[]
		self.Obs_set=[]
		self.tmp_obs = "TMP_surface"
		self.ens_timeindex=0
		self.obs_timeindex=0
		self.obs_precision = obs_precision
		self.time_periods = time_periods
		self.time_period_format = time_period_format


		for i in range(number_of_ensembles):
			fn_start = "wrfout_ens"
			fn_end = self.dt.strftime('%Y-%m-%d') + "_00"
			file = fn_start + str(i) +"_"+ fn_end
			filepath = os.path.join(self.path_to_ens,file)
			print filepath
			self.DS.append(xr.open_dataset(filepath))

		for i in range(number_of_obs):
			gfsfile = "fin_gfs_" + (str((i*6)).zfill(2))
			gfs_path = os.path.join(self.path_to_ens,gfsfile)
			print gfs_path
			self.Obs_set.append(xr.open_dataset(gfs_path))

	def WRF_bounds(self,dataset):
			wrf_size = (dataset.XLAT[0].shape[0]) - 1
			lat_max = dataset.XLAT[0][wrf_size].values[0]
			lat_min = dataset.XLAT[0][0].values[0]
			lon_max = dataset.XLONG[0][0].values[-1]
			lon_min = dataset.XLONG[0][0].values[0]
			dims = [[lat_min,lat_max],[lon_min,lon_max]]
			return dims

		#this rounds the numbers
	def round_bound(self,coord):
		return round((coord * self.obs_precision) / self.obs_precision)



	#Here we extract only the part of the observation dataset that we care about (ie what is associated with the WRF output)
	def Obs_rescale(self,obs,dims,feature,timepd,lon_adj,obs_timeindex):
		MyGrid = obs[obs_timeindex][feature][timepd]
		lat_min = dims[0][0]
		lat_max = dims[0][1]
		lon_min = dims[1][0]+lon_adj
		lon_max = dims[1][1]+lon_adj
		latrange = np.linspace(lat_min,lat_max,((lat_max-lat_min)*self.obs_precision+1))
		lonrange = np.linspace(lon_min,lon_max,((lon_max-lon_min)*self.obs_precision+1))
		MyGrid = MyGrid.sel(latitude=latrange)
		MyGrid = MyGrid.sel(longitude=lonrange)

		return MyGrid,latrange,lonrange
		


	#Here we make a grid of the points of the coordinates (lat,lon) that are affiliated with the WRF ensemble
	def CoordGrid(self,dataset):
		lat_list=[]
		lon_list=[dataset["XLONG"][0][0].values]
		lon_list = [x for x in lon_list[0]]
		for lat in dataset["XLAT"][0]:
			lat_list.append((lat[0].values+0))	
		coord_grid = []
		for lat in lat_list:
			for lon in lon_list:
				coord_grid.append([lat,lon])
		return coord_grid


	#This is a grid of the values associated with the above latlon coordinates in the dataset
	def ValuesGrid(self,dataset):
		values_grid=[]	
		for x in dataset:
			for y in x:
				values_grid.append((y.values+0))
		return values_grid

	#Create xi
	def create_xi(self,latrange,lonrange):
		xi=[]
		for i in latrange:
			for j in lonrange:
				xi.append([i,j])
		return xi



	# Take the griddata output and put it back into gridded form (instead of a list)
	def backToAGrid(self,newgrid,width):
		output_grid = []
		count = len(newgrid)/width
		for i in range(count):
			range_1 = i*width
			range_2 = (i+1)*width
			latlist = newgrid[range_1:range_2]
			output_grid.append(latlist)
		return output_grid


	def Ensemble_Differences(self,ens,obs,lon_adj,field,ens_timeindex,obs_timeindex):
		print "running ensemble diffs calc, please wait..."
		print "output summary: index 0 = regridded data, index 1 = rescaled gfs, index 2 = ens diffs"
		#Find latlon boundaries
		dims = self.WRF_bounds(ens)
		#Round the boundaries
		rd_dims = [map(self.round_bound,dims[0]),map(self.round_bound,dims[1])]
		#Narrow the scope of the GFS dataset
		GFS_Rescale_Output = self.Obs_rescale(obs,rd_dims,field,0,lon_adj,obs_timeindex)
		Rescaled_GFS = GFS_Rescale_Output[0]
		latrange = GFS_Rescale_Output[1]
		lonrange = GFS_Rescale_Output[2]
		width = len(lonrange)
		# Create the point grid
		point_grid = self.CoordGrid(ens.T2)
		# Ca[reate the value grid
		value_grid = self.ValuesGrid(ens.T2[ens_timeindex])
		# Create the Xi grid
		xi_grid = self.create_xi(latrange,lonrange-lon_adj)
		# Create the regridded_data
		Regrid_Intermediate = griddata(point_grid,value_grid,xi_grid,method = 'linear')
		Regridded_Data = np.asarray(self.backToAGrid(Regrid_Intermediate,width))
		Ens_Obs_Diff = Regridded_Data - Rescaled_GFS
		return Regridded_Data, Rescaled_GFS, Ens_Obs_Diff, point_grid, value_grid, xi_grid, Regrid_Intermediate, lonrange

	def Ens_Deux(self,ens_diff):
		lonrange = ens_diff[7]
		regridded_data = griddata(ens_diff[3],ens_diff[4],ens_diff[5],method = 'linear')
		regrid_2 = np.asarray(self.backToAGrid(regridded_data,len(lonrange)))
		Rescaled_GFS = ens_diff[1]
		Regridded_Data = regrid_2
		Ens_Obs_Diff = Regridded_Data - Rescaled_GFS
		return Regridded_Data, Rescaled_GFS, Ens_Obs_Diff


	def Create_Diffs(self,DS_Mat,ens_timeindex,obs_timeindex):
		Diffs_Mat = []
		Ens_Mat=[]
		GFS_Mat=[]
		for item in DS_Mat:
			Ens_process=self.Ens_Deux(self.Ensemble_Differences(item,self.Obs_set,0,self.tmp_obs,ens_timeindex,obs_timeindex))
			diffs = Ens_process[2]
			gfs = Ens_process[1]
			ensembles = Ens_process[0]
			Ens_Mat.append(ensembles)
			GFS_Mat.append(gfs)
			Diffs_Mat.append(diffs)
		return Diffs_Mat, Ens_Mat, GFS_Mat


	def timesToIndices(self):
		time_periods = self.time_periods
		time_period_format = self.time_period_format
		indexList=[]
		numlist = []
		enslist = []
		obslist = []
		for item in time_periods:
			numlist.append(int(item[:2]))
		for time in numlist:
			enslist.append(time/time_period_format[0])
			obslist.append(time/time_period_format[1])
		return enslist,obslist

	def PostProcess(self):
		PostProcessDict = {}
		IndexTimes = self.timesToIndices()
		time_periods = self.time_periods
		for i,period in enumerate(time_periods):
			PostProcessDict[period] = self.Create_Diffs(self.DS,IndexTimes[0][i],IndexTimes[1][i])
		return PostProcessDict


		

def Run_PostProcess(startdate,n_of_days):
	date_list = []
	date_list.append(startdate)
	for day in range(n_of_days):
		start = parser.parse(startdate)
		newday = (start + timedelta(days = day)).strftime('%Y-%m-%d')
		date_list.append(newday)
	PProcess_Array=[]
	for date in date_list:
		A = PostProcess(date,4,4,4,["00h","06h","12h","18h"],[3,6])
		PProcess_Array.append(A)
	Results_Array =[]
	for item in PProcess_Array:
		Results_Array.append(item.PostProcess())
	return Results_Array

	
def avg_err(hours,mat):
    count = 0
    summat= 0
    for item in mat:
        for i in range(4):
            summat=summat+ item[hours][0][i]
            count = count+1
    return summat/count

def std_error(hours, mat):
	emptymat = np.zeros(shape(mat[0][hours][0][1]),dtype=object)
	for i, one in enumerate(emptymat):
		for j,two in enumerate(one):
			emptymat[i,j]= []
	for i, lat in enumerate(mat[0][hours][0][1]):
		for j, lon in enumerate(lat):
			emptymat[i,j].append(float(lon))
	for day in mat:
		for i in day[hours][0]:
			for j,lat in enumerate(i):
				for k,lon in enumerate(lat):
					emptymat[j,k].append(float(lon))
	std_err_mat = np.zeros(shape(mat[0][hours][0][1]),dtype=object)
	get_std = lambda x,y: np.std(emptymat[x,y])
	for i, one in enumerate(std_err_mat):
		for j,two in enumerate(std_err_mat):
			std_err_mat[i,j] = get_std(i,j)

	return emptymat, std_err_mat









			



Nov20 = PostProcess("11-20-2016",4,4,4,["00h","06h","12h","18h"],[3,6])



