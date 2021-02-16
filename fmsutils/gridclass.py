import numpy as np
import xarray as xr
from cartopy.util import add_cyclic_point
from fmsutils.fort_interp import fort_interp_mod as fim
from scipy.integrate import quad
from scipy.interpolate import interp1d
from windspharm.standard import VectorWind

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
    
    def calc_divrot(self,R_p):
        """Calculate the rotational and divergent components 
        of the velocity field"""
        
        u = self.data['ucomp'].data
        v = self.data['vcomp'].data
        
        f_udiv = np.zeros_like(u)
        f_vdiv = np.zeros_like(u)
        f_urot = np.zeros_like(u)
        f_vrot = np.zeros_like(u)
        
        if 'time' in self.data['ucomp'].coords:            
            for l in range(self.npz):
                for t in range(self.nt):
                    f_vw = VectorWind(u[t,l], v[t,l], rsphere=R_p)
                    f_udiv_l, f_vdiv_l, f_urot_l, f_vrot_l = f_vw.helmholtz(truncation=21)
                    f_udiv[t,l] = f_udiv_l
                    f_vdiv[t,l] = f_vdiv_l
                    f_urot[t,l] = f_urot_l
                    f_vrot[t,l] = f_vrot_l
            
        else:
            for l in range(self.npz):
                f_vw = VectorWind(f_u[l], f_v[l], rsphere=R_p)
                f_udiv_l, f_vdiv_l, f_urot_l, f_vrot_l = f_vw.helmholtz(truncation=21)
                f_udiv[l] = f_udiv_l
                f_vdiv[l] = f_vdiv_l
                f_urot[l] = f_urot_l
                f_vrot[l] = f_vrot_l                
                
        udiv = xr.DataArray(data   = f_udiv,
                            coords = self.data['ucomp'].coords,
                            dims   = self.data['ucomp'].dims)
        
        urot = xr.DataArray(data   = f_urot,
                            coords = self.data['ucomp'].coords,
                            dims   = self.data['ucomp'].dims)
        
        vdiv = xr.DataArray(data   = f_vdiv,
                            coords = self.data['vcomp'].coords,
                            dims   = self.data['vcomp'].dims)
        
        vrot = xr.DataArray(data   = f_vrot,
                            coords = self.data['vcomp'].coords,
                            dims   = self.data['vcomp'].dims)
        
        self.data['udiv'] = udiv; self.data['urot'] = urot
        self.data['vdiv'] = vdiv; self.data['vrot'] = vrot
        
    def eddy(self, dataarray, dim):
        return dataarray - dataarray.mean(dim)
    
    def mpsi(self, R_p, g, tav=False):
        """Calculate the meridional mass streamfunction, returns streamfunction
           along with the contribution function to this integral"""
        if tav == True:
            vmean = self.data['vcomp'].mean(('time','lon'))
        else:
            vmean = self.data['vcomp'].mean('lon')
            
        dp  = self.data['ph_int'].diff('ph_int')
        dp  = dp.swap_dims({'ph_int':'pf_int'})
        cos = xr.ufuncs.cos(xr.ufuncs.radians(self.data['lat']))*2*np.pi*R_p/g
        
        contrib = cos*mean*dp
        result  = contrib.cumsum('pf_int')
        
        result = result.pad({'pf_int':(1,0)})
        
        result = result.swap_dims({'pf_int':'ph_int'})
        
        self.data['mpsi'] = result
        self.data['mpsi_contrib'] = contrib
                  
        
        

        
        
        