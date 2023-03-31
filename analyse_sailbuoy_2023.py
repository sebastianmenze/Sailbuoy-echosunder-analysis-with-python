# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 12:20:34 2023

@author: Administrator
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 10:53:49 2021

@author: Administrator
"""


from echolab2.instruments import EK80, EK60

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import glob 
import os

from scipy.ndimage.filters import uniform_filter1d

from scipy.signal import convolve2d
from scipy.interpolate import interp1d

from echopy import transform as tf
from echopy import resample as rs
from echopy import mask_impulse as mIN
from echopy import mask_seabed as mSB
from echopy import get_background as gBN
from echopy import mask_signal2noise as mSN
from echopy import mask_range as mRG
from echopy import mask_shoals as mSH


raw_files_folder=r'C:\Users\a5278\Documents\postdoc_krill\2023_cruise\Sailbuoy\EK80 data\raw_files'

os.chdir(raw_files_folder)
rawfiles= np.sort( glob.glob(  '*.raw'  ) )    
# rawfiles= np.sort( glob.glob(  os.path.join(raw_files_folder, '*.raw')  ) )    


krill_time=np.array([],dtype='datetime64[ms]')
krill_nasc=np.array([])
krill_lat=np.array([])
krill_lon=np.array([])

# rawfile=rawfiles[0]

df_wbt = pd.read_csv(r'C:\Users\a5278\Documents\postdoc_krill\2023_cruise\Sailbuoy\EK80 data\WBT.TXT', delim_whitespace=True,header=32)

a = df_wbt.iloc[:,5] + ' '+     df_wbt.iloc[:,6] 
df_wbt['time'] = pd.to_datetime(a,format='%d.%m.%Y %H:%M:%S')

dfp=pd.DataFrame([])

dfp_sw=pd.DataFrame([])

for rawfile in rawfiles:  
    # if not os.path.isfile( rawfile[:-4]+'.npy' ):
        try:
            
            
           raw_obj = EK80.EK80()
           raw_obj.read_raw(rawfile)
            
           print(raw_obj)
            
           raw_data = raw_obj.raw_data['EKA 268629-07 ES200-7CDK-Split'][0]
            
           cal_obj = raw_data.get_calibration()
            # Get sv values
           sv_obj = raw_data.get_sv(calibration = cal_obj)
            # Get sv as depth
            #sv_obj_as_depth = raw_data.get_sv(calibration = cal_obj,
            #    return_depth=True)
            
           # positions = raw_obj.nmea_data.interpolate(sv_obj, 'GLL')[1]
            # positions['latitude']
                    
            # meter_dif=geopy.distance.distance( (lsss_sv['Latitude'].iloc[0],lsss_sv['Longitude'].iloc[0]), (lsss_sv['Latitude'].iloc[-1],lsss_sv['Longitude'].iloc[-1]) ).km * 1000     
            
        
            # Get frequency label
           freq = sv_obj.frequency
            
            # Expand sv values into a 3d object
           data3d = np.expand_dims(sv_obj.data, axis=0)
           sv= np.flip( np.rot90( 10*np.log10( data3d[0,:,:] ) ) )
            
           # plt.figure(0)
           # plt.clf()
            
           # plt.imshow( sv  )
           # plt.clim([-90,-30])
           # plt.colorbar()
           
           Sv120=sv
           r120=sv_obj.range
           
           distancecovered_guess=(2*0.514444 * (60*10)) /1000
           km120=np.linspace(0,distancecovered_guess,sv.shape[1]) # guess 10min rec 2knot speed (2*0.514444 * (60*10)) /1000
           t120=sv_obj.ping_time
           
           # get loc
           t_ek80 = pd.Series(t120)
           ix_start = np.argmin( np.abs(df_wbt['time'] -  t_ek80.min() )  )        
           ix_end = np.argmin( np.abs(df_wbt['time'] -  t_ek80.max() )  )   
           
           lat=   df_wbt.iloc[ix_start,7]
           lon=   df_wbt.iloc[ix_start,8]
           
         # alpha120=40*np.ones(r120.shape)  # absoprtion coef
         
           #--------------------------------------------------------------------------       
           # Clean impulse noise      
           # Sv120in, m120in_ = mIN.wang(Sv120, thr=(-70,-30), erode=[(3,3)],
           #                          dilate=[(7,7)], median=[(7,7)])
           # Sv120in, m120in_ = mIN.wang(Sv120, thr=(-70,-40), erode=[(3,3)],
           #                    dilate=[(7,7)], median=[(7,7)])
           Sv120in=sv.copy()  
        # -------------------------------------------------------------------------
         # estimate and correct background noise       
           p120           = np.arange(len(t120))                
           s120           = np.arange(len(r120))  
         
                       
           bn120, m120bn_ = gBN.derobertis(Sv120, s120, p120, 5, 20, r120, 0.044926) # whats correct absoprtion?
           Sv120clean     = tf.log(tf.lin(Sv120in) - tf.lin(bn120))
          
           # plt.figure(num=1)
           # plt.clf()
           # plt.subplot(311)
           # plt.imshow((Sv120),aspect='auto')
           # plt.clim([-82,-40])
           # plt.colorbar()
           # plt.subplot(312)
           # plt.imshow((Sv120in),aspect='auto')
           # plt.clim([-82,-40])
           # plt.colorbar()
           # plt.subplot(313)
           # plt.imshow((Sv120clean),aspect='auto')
           # plt.clim([-82,-40])
           # plt.colorbar()          
         
         # -------------------------------------------------------------------------
         # mask low signal-to-noise 
           m120sn             = mSN.derobertis(Sv120clean, bn120, thr=12)
           Sv120clean[m120sn] = np.nan
           
           # Sv120clean=sv.copy()
           # Sv120clean[Sv120clean<-55]=-999
           
        # get mask for seabed
           m120sb = mSB.ariza(Sv120, r120, r0=20, r1=1000, roff=0,
                              thr=-38, ec=1, ek=(3,3), dc=10, dk=(5,15))
          # m120sb = mSB.ariza(Sv120, r120, r0=20, r1=1000, roff=0,
          #                   thr=-38, ec=1, ek=(3,3), dc=10, dk=(3,7))
           # m120sb = mSB.ariza(Sv120, r120, r0=20, r1=1000, roff=0,
           #                 thr=-38, ec=1, ek=(3,5), dc=10, dk=(3,10))           
           
           Sv120clean[m120sb]=-999
                
         
           # -------------------------------------------------------------------------
           # get swarms mask
           k = np.ones((3, 3))/3**2
           Sv120cvv = tf.log(convolve2d(tf.lin(Sv120clean), k,'same',boundary='symm'))   
         
         
           # plt.figure(num=1)
           # plt.clf()
           # plt.subplot(211)
           # plt.imshow((sv),aspect='auto')
           # plt.clim([-82,-50])
           # plt.colorbar()
           # plt.subplot(212)
           # plt.imshow((Sv120cvv),aspect='auto')
           # plt.clim([-82,-50])
           # plt.colorbar()
          
           # p120           = np.arange(np.shape(Sv120cvv)[1])                
           # s120           = np.arange(np.shape(Sv120cvv)[0])           
                 
           m120sh, m120sh_ = mSH.echoview(Sv120cvv, s120, p120, thr=-70,
                                     mincan=(3,10), maxlink=(3,15), minsho=(3,15))
          
         
    
           
         # -------------------------------------------------------------------------
         # get Sv with only swarms
           Sv120sw                    = Sv120clean.copy()
           Sv120sw[~m120sh] = np.nan
          
           ixdepthvalid=(r120<250) & (r120>20)
           Sv120sw[~ixdepthvalid,:]=np.nan
           
           td=t_ek80.max() - t_ek80.min()

           nasc_time=t120
           cellthickness=np.abs(np.mean(np.diff( r120) )) 
           # nasc=4*np.pi*1852**2 * np.nansum( np.power(10, Sv120sw /10)*cellthickness ,axis=0)    
           int_sv= np.nansum( np.power(10, Sv120sw /10) ,axis=0)    
          
           a =np.power(10, Sv120sw /10)
           a[np.isnan(a)]=0
           profile_sw= np.nanmean(a,axis=1  )
           profile_sw[profile_sw<0]=0
 
           a =np.power(10, Sv120clean /10)
           # a[~ixdepthvalid,:]=np.nan
           a[np.isnan(a)]=0
           profile= np.nanmean(a,axis=1  )
           profile[profile<0]=0
           

           # profile_sw= np.nanmean(Sv120sw,axis=1  )
           # profile_sw[np.isnan(profile_sw)]=0
 
           # profile= np.nanmean(Sv120clean,axis=1  )
           # profile[np.isnan(profile)]=0
           
           # dfp=pd.DataFrame([])
           # dfp['depth']=r120
           # dfp['NASC']=profile
           dfp=pd.concat([dfp,pd.DataFrame(profile)],axis=1)
           dfp_sw=pd.concat([dfp_sw,pd.DataFrame(profile_sw)],axis=1)
           
           # fig=plt.figure(num=1)
           # fig.set_size_inches(10,10)
           # plt.clf()
           # plt.subplot(211)
           # plt.imshow((Sv120),aspect='auto',extent=[0, td.total_seconds(), -r120.max(),0 ])
           # plt.clim([-80,-40])
           # plt.xlabel('Time in seconds')
           # plt.ylabel('Depth in m')
           # plt.title(rawfile)
           
           # plt.colorbar(label='s_V in dB')
           
           # plt.subplot(223)
           # plt.imshow((Sv120sw),aspect='auto',extent=[0, td.total_seconds(), -r120.max(),0 ])
           # plt.clim([-80,-40])
           # plt.xlabel('Time in seconds')
           # plt.ylabel('Depth in m')
           # plt.title('Swarm detection')
           # plt.colorbar(label='s_V in dB')
           
           # plt.subplot(224)
           # plt.plot( 10*np.log10(profile),-r120  ,'-k',label='all')
           # plt.plot(  10*np.log10(profile_sw),-r120  ,'-b',label='swarm detection')
           # plt.xlabel('Avg. s_V')
           # plt.ylabel('Depth in m')
           # plt.title('Avg. NASC profile')
           # plt.grid()
           # plt.legend()
           # plt.ylim([-200,0])
 
           
           # plt.tight_layout()

           # plt.savefig(rawfile[0:-3]+'png' )
               
        
                         
           
           krill_time=np.append(krill_time,nasc_time[0])
           krill_nasc=np.append(krill_nasc,np.mean(int_sv)    )
           
           krill_lat=np.append(krill_lat,lat)
           krill_lon=np.append(krill_lon,lon)     
           
           del raw_obj
        except Exception as e:
            print(e)
  
#%%


plt.figure(num=0)
plt.clf()
plt.plot( krill_time , krill_nasc,'.b')

nasc_filt = uniform_filter1d(krill_nasc, size=100)
plt.plot(krill_time,nasc_filt,'-r')

#%%

df =pd.DataFrame( np.transpose(dfp) )

df.index= pd.to_datetime(krill_time)
df.columns= r120

# df.to_csv('volumne_backscatter_profiles.csv')

track=pd.DataFrame([])
track['lat']=krill_lat
track['lon']=krill_lon
# track['integrated_sv']=krill_nasc
track.index=pd.to_datetime(krill_time)


ix_range= df.columns>20
track['integrated_sv']=  df.iloc[:,ix_range].sum(axis=1)  


# track.to_csv('sailbuoy_track_and_sv.csv')

#%%

plt.figure(num=0)
plt.clf()

plt.plot(track['integrated_sv'],'.k' )


plt.plot(track['integrated_sv'].resample('4h').mean() ,'-r',label='4h avg.')

t1=pd.to_datetime('2023-01-13 00:00:00')
t2=pd.to_datetime('2023-01-28 12:00:00')
plt.xlim([t1,t2  ])

plt.ylabel('Depth integrated volume backscatter')

plt.legend()

# plt.savefig('integrated_sv_timeseries.png')

#%%
# td= pd.to_timedelta( krill_time.max() - krill_time.min() )
# td.total_seconds()

# days=td.total_seconds() / (24*60*60)

# c=  10*np.log10(dfp)
# c[np.isinf(c)]=-999

# plt.figure(num=5)
# plt.clf()

# plt.imshow(c,aspect='auto',extent=[0,days,-r120.max(),0])
# plt.clim([-80,-40])

# plt.colorbar(label='Avg. volume backscatter')

#%%

import matplotlib.dates as mdates

t1=pd.to_datetime('2023-01-13 00:00:00')
t2=pd.to_datetime('2023-01-28 12:00:00')

mt = mdates.date2num((t1,t2))


ix = (df.index>=t1) & (df.index<=t2)

c=  10*np.log10( np.transpose( df.iloc[ix,:] )  )
c[np.isinf(c)]=-999

plt.figure(num=5)
plt.clf()

thr=16

ax=plt.subplot()


plt.imshow(c,aspect='auto',extent=[mt[0],mt[-1],-r120.max(),0])
plt.clim([-80,-40])

plt.colorbar(label='Avg. volume backscatter in dB')
plt.ylabel('Depth in m')
ax.xaxis_date()

# plt.savefig('sv_timeseries_sailbuoy_2023.jpg', dpi=300)

#%%

td= pd.to_timedelta( krill_time.max() - krill_time.min() )

days=td.total_seconds() / (24*60*60)

c=  10*np.log10(dfp_sw)
c[np.isinf(c)]=-999

plt.figure(num=5)
plt.clf()

plt.imshow(c,aspect='auto',extent=[0,days,-r120.max(),0])
plt.clim([-80,-40])

plt.colorbar(label='s_v')


#%% load maps

from netCDF4 import Dataset
import pandas as pd
import cartopy
import cartopy.crs as ccrs
import numpy as np
import matplotlib.pyplot as plt
import utm

# load data and slice out region of interest

# read mapdata

latlim=[-60.7,-60]
lonlim=[-47.5,-43.5]

spacer=1
gebcofile=r"C:\Users\a5278\Documents\gebco_2020_netcdf\GEBCO_2020.nc"
gebco = Dataset(gebcofile, mode='r')
g_lons = gebco.variables['lon'][:]
g_lon_inds = np.where((g_lons>=lonlim[0]) & (g_lons<=lonlim[1]))[0]
# jump over entries to reduce data
g_lon_inds=g_lon_inds[::spacer]

g_lons = g_lons[g_lon_inds].data
g_lats = gebco.variables['lat'][:]
g_lat_inds = np.where((g_lats>=latlim[0]) & (g_lats<=latlim[1]))[0]
# jump over entries to reduce data
g_lat_inds=g_lat_inds[::spacer]

g_lats = g_lats[g_lat_inds].data
d = gebco.variables['elevation'][g_lat_inds, g_lon_inds].data
gebco.close()

#%%

# plt.savefig('map_1.jpg',dpi=300)

#%%


fig=plt.figure(num=3)
plt.clf()

central_lon= lonlim[0]+(lonlim[1]-lonlim[0])/2
central_lat = latlim[0]+(latlim[1]-latlim[0])/2
extent = [lonlim[0],lonlim[1], latlim[0],latlim[1]]
#ax = plt.axes(projection=ccrs.PlateCarree(central_longitude= lonlim[0]+(lonlim[1]-lonlim[0])/2 ))

ax = plt.axes(projection=ccrs.Orthographic(central_lon, central_lat))

ax.set_extent(extent)
  
ax.gridlines(draw_labels=True)
#ax.coastlines(resolution='50m')
#ax.add_feature(cartopy.feature.LAND)

d_plot=d
d_plot[d<-4000]=-4000

plt.contourf(g_lons, g_lats, d_plot, np.arange(-4000,0,100),cmap='Blues_r',
                  linestyles=None, transform=ccrs.PlateCarree())
CS=plt.contour(g_lons, g_lats, d, [-2000,-1000,-500],colors='k',linewidth=.1,
                  linestyles='-', transform=ccrs.PlateCarree())
plt.clabel(CS, inline=True, fontsize=10, fmt='%i')

CS=plt.contourf(g_lons, g_lats, d, [0,8000],colors='silver',linewidth=1,
                  linestyles='-', transform=ccrs.PlateCarree())

CS=plt.contour(g_lons, g_lats, d, [0],colors='k',linewidth=1,
                  linestyles='-', transform=ccrs.PlateCarree())





plt.scatter(krill_lon,krill_lat,10,krill_nasc,transform=ccrs.PlateCarree() )
plt.clim([0,500])
plt.colorbar(label='NASC')

plt.tight_layout()


#%%


