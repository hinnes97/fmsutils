import numpy as np
import xarray as xr
from cartopy.util import add_cyclic_point
from fmsutils.fort_interp import fort_interp_mod as fim
from scipy.integrate import quad
from scipy.interpolate import interp1d

class grid:
    def __init__(self, data):
        self.lon = data['grid_xt'][:]
        self.lat = data['grid_yt'][:]
        self.time = data['time'][:]
        self.pk = data['pk'][:]
        self.bk = data['bk'][:]
        self.ps = data['ps'][:][:,:,:]

        self.data = data
        
        self.npx = len(self.lon)
        self.npy = len(self.lat)
        self.npz  = len(self.pk) - 1
        self.nt = len(self.time)
        
        self.eq = np.searchsorted(self.lat,0)
        self.sub = np.searchsorted(self.lon,0)
        self.anti = np.searchsorted(self.lon,180)
        self.wterm =  np.searchsorted(self.lon,90)
        self.eterm = np.searchsorted(self.lon, 270)

        self.plevels()

        self.data = self.data.rename({'grid_yt':'lat'})
        self.data = self.data.rename({'grid_xt':'lon'})

    def plevels(self):
        """ Work out the 3D pfull and phalf grids and set up the levels to be interpolated to"""
        npx = self.npx
        npy = self.npy
        npz = self.npz

        # Blocks below calculate pfull and phalf grid and define them on new dimensions ph_level
        # and pf_level, which is just a level index
        self.ph_levels = np.arange(0.5,npz+2-0.5)
        self.pf_levels = np.arange(1,npz+1)

        self.phalf = self.pk + self.ps*self.bk
        self.phalf= self.phalf.transpose("time", "phalf", "grid_yt", "grid_xt")
        
        self.phalf = self.phalf.swap_dims({"phalf":"ph_level"})
        del self.phalf["phalf"]
        self.phalf["ph_level"] = ("ph_level",self.ph_levels)

        self.pfull = self.phalf.diff("ph_level")/np.log(self.phalf).diff("ph_level")

        self.pfull["pf_level"] = ("ph_level", self.pf_levels)
        self.pfull = self.pfull.swap_dims({"ph_level":"pf_level"})
        del self.pfull['ph_level']

        self.data['p_half'] = self.phalf
        self.data['p_full'] = self.pfull

        # Levels to interpolate to (maybe change to be a(k) + b(k)*p_ref instead of log levels?)
        pfmin = np.amin(self.data['p_full'].data)
        pfmax = np.amax(self.data['p_full'].data)
        phmin = np.amin(self.data['p_half'].data)
        phmax = np.amax(self.data['p_half'].data)
        
        self.pfi = np.logspace(np.log10(pfmin), np.log10(pfmax), self.data['p_full'].shape[1])
        self.phi = np.logspace(np.log10(phmin), np.log10(phmax), self.data['p_half'].shape[1])
        
        self.data = self.data.assign_coords({'ph_int':self.phi, 'pf_int':self.pfi})
        
    def interp_var(self,var):
        """Interpolate variables to pressure levels pfi and phi"""
        try:
            variable = self.data[var]
        except:
            print('var not in dataset')
            
        if 'pfull' in variable.coords:
            variable.data = fim.interp_data(variable.data,self.data['p_full'].data, self.pfi)
            variable = variable.swap_dims({'pfull':'pf_int'})
            variable['pf_int'] = ('pf_int', self.pfi)
            del variable['pfull']
            self.data[var] = variable

        elif 'phalf' in variable.coords:
            variable.data = fim.interp_data(variable.data,self.data['p_half'].data, self.phi)
            variable = variable.swap_dims({'phalf':'ph_int'})
            variable['ph_int'] = ('ph_int', self.phi)
            del variable['phalf']
            self.data[var] = variable
            
    def interp_all(self):
        for var in self.data.data_vars:
            if ('pfull' in self.data[var].coords) or ('phalf' in self.data[var].coords):
                if (var != 'pk') and (var != 'bk'):
                    self.interp_var(var)
                    print(f"{var} interpolated")


    def wrap(self, var):
        """ Use cartopy function to wrap longitude coordinate so that we don't get an annoying white line at 0 degrees"""
        if 'lon' in self.data[var].coords:
            lon_idx = self.data[var].dims.index('lon')
            wrap_dat,wrap_lon = add_cyclic_point(self.data[var].data, coord=self.data['lon'].data, axis=lon_idx)
            
            attrs = self.data[var]['lon'].attrs
            coords = self.data[var].coords
            del coords['lon']
            coords = {**coords, 'wlon':('wlon',wrap_lon)}               

            dims = (*self.data[var].dims[:-1], 'wlon')
                
            new_da = xr.DataArray(data=wrap_dat,
                                  coords=coords,
                                  attrs = self.data[var].attrs,
                                  dims=dims)
            new_da['wlon'].attrs = attrs
            self.data[var] = new_da
            
    def unwrap(self, var):
        """Remove the wrapped longitude coordinate"""
        if 'wlon' in self.data[var].coords:
            lon_idx = self.data[var].dims.index('wlon')
            newlon = self.data['wlon'][:-1] 
            newdat=np.take(self.data[var].data,
                           axis=lon_idx,indices=range(self.npx))
            attrs = self.data[var]['wlon'].attrs
            coords = self.data[var].coords
            del coords['wlon']
            coords = {**coords, 'lon':('lon',newlon)}
            dims = (*self.data[var].dims[:-1], 'lon')
            new_da = xr.DataArray(data=newdat,
                                  coords=coords,
                                  attrs = self.data[var].attrs,
                                  dims=dims)
            new_da['lon'].attrs = attrs
            self.data[var] = new_da
    
    def height(self, R, g):
        """Calculate the height of pressure layers"""

        h = fim.height(self.data['temp'].data,self.data['p_half'].data,R,g)

        coords = self.data['temp'].coords
        del coords['pfull']
        coords = {**coords, 'phalf':self.data['phalf']}
        dims = ('time', 'phalf', 'lat', 'lon')
        h = xr.DataArray(data=h,
                         coords=coords,
                         dims=dims,
                          )
        self.data['h'] = h
