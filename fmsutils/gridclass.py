import numpy as np
import xarray as xr
from cartopy.util import add_cyclic_point
from fmsutils.fort_interp import fort_interp_mod as fim
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.interpolate import griddata
from windspharm.standard import VectorWind

class grid:
    def __init__(self, data, lon_ss=0):
        
        self.data = data

        if 'grid_xt' in data.coords:
            self.data = self.data.rename({'grid_xt':'lon'})
            self.npx = len(self.data['lon'])
        elif 'lon_TL' in data.coords:
            self.npx = len(self.data['lon_TL'])
                           
        if 'grid_yt' in data.coords:
            self.data = self.data.rename({'grid_yt':'lat'})
            self.npy = len(self.data['lat'])
        elif 'lat_TL' in data.coords:
            self.npy = len(self.data['lat_TL'])
            
        if 'time' in data.coords:
            self.nt = len(data['time'])
        
        self.pk = self.data['pk'][:]
        self.bk = self.data['bk'][:]
        self.ps = self.data['ps'][:][:,:,:]
        self.lon_ss = lon_ss

        self.npz  = len(self.pk) - 1

        self.plevels()
        
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

        if 'lon' in self.data.coords and 'lat' in self.data.coords:
            self.phalf= self.phalf.transpose("time", "phalf", "lat", "lon")
        elif 'lon_TL' in self.data.coords and 'lat_TL' in self.data.coords:
            self.phalf=self.phalf.transpose('time', 'phalf','lat_TL','lon_TL')
        
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
        
    def interp_vars(self,*variables):
        
        for var in variables:
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
            l = 'lon'
        elif 'lon_TL' in self.data[var].coords:
            l = 'lon_TL' 
        else:
            raise KeyError('No lon or lon_TL in coords')
            
        lon_idx = self.data[var].dims.index(l)
        wrap_dat,wrap_lon = add_cyclic_point(self.data[var].data, coord=self.data[l].data, axis=lon_idx)

        attrs = self.data[var][l].attrs
        coords = self.data[var].coords
        del coords[l]
        coords = {**coords, 'w'+l:('w'+l,wrap_lon)}               

        dims = (*self.data[var].dims[:-1], 'w'+l)

        new_da = xr.DataArray(data=wrap_dat,
                              coords=coords,
                              attrs = self.data[var].attrs,
                              dims=dims)
        new_da['w'+l].attrs = attrs
        self.data[var] = new_da
            
    def unwrap(self, var):
        """Remove the wrapped longitude coordinate"""
        if 'lon' in self.data[var].coords:
            l = 'lon'
        elif 'lon_TL' in self.data[var].coords:
            l = 'lon_TL' 
        else:
            raise KeyError('No lon or lon_TL in coords')
            
        if 'w'+l in self.data[var].coords:
            lon_idx = self.data[var].dims.index('w'+l)
            newlon = self.data['w'+l][:-1] 
            newdat=np.take(self.data[var].data,
                           axis=lon_idx,indices=range(self.npx))
            attrs = self.data[var]['w'+l].attrs
            coords = self.data[var].coords
            del coords['w'+l]
            coords = {**coords, l:(l,newlon)}
            dims = (*self.data[var].dims[:-1], l)
            new_da = xr.DataArray(data=newdat,
                                  coords=coords,
                                  attrs = self.data[var].attrs,
                                  dims=dims)
            new_da[l].attrs = attrs
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
        
        # N.B. latitude has to be ordered North to South so flip, then reflip
        lataxis = self.data['ucomp'].dims.index('lat')
        
        u = np.flip(self.data['ucomp'].data,axis=lataxis)
        v = np.flip(self.data['vcomp'].data,axis=lataxis)
        
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
                f_vw = VectorWind(u[l], v[l], rsphere=R_p)
                f_udiv_l, f_vdiv_l, f_urot_l, f_vrot_l = f_vw.helmholtz(truncation=21)
                f_udiv[l] = f_udiv_l
                f_vdiv[l] = f_vdiv_l
                f_urot[l] = f_urot_l
                f_vrot[l] = f_vrot_l                
                
        udiv = xr.DataArray(data   = np.flip(f_udiv,axis=lataxis),
                            coords = self.data['ucomp'].coords,
                            dims   = self.data['ucomp'].dims)
        
        urot = xr.DataArray(data   = np.flip(f_urot,axis=lataxis),
                            coords = self.data['ucomp'].coords,
                            dims   = self.data['ucomp'].dims)
        
        vdiv = xr.DataArray(data   = np.flip(f_vdiv,axis=lataxis),
                            coords = self.data['vcomp'].coords,
                            dims   = self.data['vcomp'].dims)
        
        vrot = xr.DataArray(data   = np.flip(f_vrot,axis=lataxis),
                            coords = self.data['vcomp'].coords,
                            dims   = self.data['vcomp'].dims)
        
        self.data['udiv'] = udiv; self.data['urot'] = urot
        self.data['vdiv'] = vdiv; self.data['vrot'] = vrot
        
        #return udiv, urot, vdiv, vrot
        
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
                        
        
    def new_TL_grid(self, lon_TL, lat_TL, *variables):
            """Returns a grid class with all the variables in tidally 
               locked coords"""
            
            # Set up an empty dataset
            ds = xr.Dataset()
            
            filled_ds = self.tidal_locked(self, lon_TL, lat_TL, ds, *variables, 'ps')
            
            filled_ds['bk'] = self.data['bk']
            filled_ds['pk'] = self.data['pk']
               
            return self.create_TL_grid(filled_ds)
    
    def add_to_TL_grid(self, latlon_grid, *variables):
        """Converts variables from latlon grid and adds to TL 
           grid class instance"""
        ds = self.data
        
        lon_TL, lat_TL = self.data['lon_TL'].data, self.data['lat_TL'].data
        
        filled_ds = self.tidal_locked(latlon_grid, lon_TL, lat_TL, ds, *variables)
        
        self.data = filled_ds
    
    @classmethod
    def tidal_locked(cls, grd, lon_TL, lat_TL, ds, *variables):
        lat, lon = grd.data['lat'].data, grd.data['lon'].data
            
        lat_TL_orig, lon_TL_orig = cls.transform_latlon_to_TL(lat, lon,lon_ss=grd.lon_ss)

        # Treat vector winds differently
        vector_quantities = ['ucomp', 'vcomp', 'udiv','vdiv','urot','vrot']

        for var in variables:
            if var not in vector_quantities:
                if ('lon' in grd.data[var].coords) and ('lat' in grd.data[var].coords):
                    raw_array = cls.interpolate_to_TL_ndim(lat,lon,lat_TL[:,None],lon_TL[None,:],grd.data[var].data,lon_ss=grd.lon_ss,method="nearest")

                    coords = grd.data[var].coords
                    del coords['lat']; del coords['lon']
                    coords = {**coords, 'lat_TL':lat_TL, 'lon_TL':lon_TL}
                    coords = list(coords.items())

                    da = xr.DataArray(data=raw_array, 
                                      coords=coords)

                    ds[var] = da

                else:
                    ds[var] = self.data[var]
           
        for wind in [('ucomp','vcomp'), ('udiv','vdiv'), ('urot','vrot')]:
            if wind[0] in variables and wind[1] in variables:
                u,v = grd.data[wind[0]].data, grd.data[wind[1]].data
                u_TL, v_TL = cls.transform_velocities_to_TL_interp(u,v,
                                                  lat,lon,lat_TL[:,None],
                                                  lon_TL[None,:],
                                           lon_ss=grd.lon_ss,method="nearest")

                coords = grd.data[wind[0]].coords
                del coords['lat']; del coords['lon']
                coords = {**coords, 'lat_TL':lat_TL, 'lon_TL':lon_TL}
                coords = list(coords.items())

                uda = xr.DataArray(data=u_TL,
                                   coords=coords)

                vda = xr.DataArray(data=v_TL,
                                   coords=coords)

                ds[wind[0]] = uda
                ds[wind[1]] = vda

            elif (wind[0] in variables) ^ (wind[1] in variables):
                print('Warning: require both u and v velocities to transform to TL coordinates')
            
            return ds

    @classmethod
    def create_TL_grid(cls,dataset):
        return cls(dataset)              
    
    @classmethod   
    def transform_latlon_to_TL(cls,lat_tmp,lon_tmp,lon_ss=0.):
        if lat_tmp.ndim==1 or lon_tmp.ndim==1:
            lon,lat = np.meshgrid(lon_tmp,lat_tmp)
        else:
            lat,lon = lat_tmp,lon_tmp
        x = np.cos(lat) * np.cos(lon)
        y = np.cos(lat) * np.sin(lon)
        z = np.sin(lat)

        lon_ss = -lon_ss
        xprime = np.cos(lon_ss)*x - np.sin(lon_ss)*y
        yprime = np.sin(lon_ss)*x + np.cos(lon_ss)*y

        lat_TL = np.arcsin(xprime)
        lon_TL = np.arctan2(yprime,z)
        # change from lon in {-pi,pi} to lon in {0,2pi}:
        lon_TL[lon_TL<0.] = lon_TL[lon_TL<0.] + 2.*np.pi

        return lat_TL, lon_TL
    
    @classmethod
    def transform_velocities_to_TL(cls,u,v,lat_tmp,lon_tmp,lon_ss=0.):
        # This returns TL velocities but on irregular grid.
        if lat_tmp.ndim==1 or lon_tmp.ndim==1:
            lon,lat = np.meshgrid(lon_tmp,lat_tmp)
        else:
            lat,lon = lat_tmp,lon_tmp
        lat_TL,lon_TL = cls.transform_latlon_to_TL(lat,lon,lon_ss)

        lon_ss = -lon_ss

        # (tedious algebra - got these via Mathematica)
        Dlon_TL_Dlon = 1./( 1./np.cos(lon+lon_ss)*np.tan(lat)+np.cos(lat)/np.sin(lat)*np.sin(lon+lon_ss)*np.tan(lon+lon_ss) )
        Dlon_TL_Dlat = 8.*np.sin(lon+lon_ss) / ( -6.+2.*np.cos(2.*lat)+np.cos(2.*(lat-lon-lon_ss))+\
                                                  2.*np.cos(2*(lon+lon_ss))+np.cos(2.*(lat+lon+lon_ss)) )
        Dlat_TL_Dlon = -np.cos(lat)*np.sin(lon+lon_ss)/( np.sqrt(1.-np.cos(lat)**2*np.cos(lon+lon_ss)**2) )
        Dlat_TL_Dlat = -np.cos(lon+lon_ss)*np.sin(lat)/( np.sqrt(1.-np.cos(lat)**2*np.cos(lon+lon_ss)**2) )

        u_TL = Dlon_TL_Dlon * np.cos(lat_TL)/np.cos(lat)*u + Dlon_TL_Dlat*np.cos(lat_TL)*v
        v_TL = Dlat_TL_Dlon * u/np.cos(lat) + Dlat_TL_Dlat*v
        return u_TL,v_TL
    
    @classmethod
    def transform_velocities_to_TL_interp(cls, u,v,lat,lon,lat_TL,
                                          lon_TL,lon_ss=0.,method="nearest"):
        # This returns TL velocities, interpolated onto regular grid.
        u_TL,v_TL = cls.transform_velocities_to_TL(u,v,lat,lon,lon_ss)
        u_TL_i = cls.interpolate_to_TL_ndim(lat,lon,lat_TL,lon_TL,u_TL,lon_ss,method=method)
        v_TL_i = cls.interpolate_to_TL_ndim(lat,lon,lat_TL,lon_TL,v_TL,lon_ss,method=method)
        return u_TL_i,v_TL_i
    
    @classmethod
    def interpolate_to_TL(cls,lat,lon,lat_TL,lon_TL,data,
                          lon_ss=0.,method="nearest"):
        # Assume data is 2D!
        # Interpolate data, given on an earth-like lat-lon grid, to a TL lat-lon grid.
        # First, transform the given lat-lon points into TL coords.
        # Second, use the given points as basis for interpolation.
        lat_TL_given,lon_TL_given = cls.transform_latlon_to_TL(lat,lon,lon_ss)
        data_interp = griddata( (lat_TL_given.ravel(),lon_TL_given.ravel()),\
                                data.ravel(),(lat_TL,lon_TL),method=method )
        return data_interp
    
    @classmethod
    def interpolate_to_TL_ndim(cls,lat,lon,lat_TL,lon_TL,data,lon_ss=0.,method="nearest"):
        # Uses above, but for N-dim data.
        # Also assume lat/lon are last two dims.
        if data.ndim==2:
            data_interp = cls.interpolate_to_TL(lat,lon,\
                                            lat_TL[None,:],\
                                            lon_TL[:,None],data,lon_ss,method).T
        elif data.ndim==3:
            data_interp = np.zeros((data.shape[0],lat_TL.size,lon_TL.size))
            for i in range(data.shape[0]):
                data_interp[i] = \
                np.squeeze(cls.interpolate_to_TL(lat,lon,lat_TL[:,None],lon_TL[None,:],\
                                      data[i],lon_ss,method))
        elif data.ndim==4:
            data_interp = np.zeros((data.shape[0],data.shape[1], \
                                lat_TL.size,lon_TL.size))
            for i in range(data.shape[0]):
                for j in range(data.shape[1]):
                    data_interp[i,j,...] = \
                                     cls.interpolate_to_TL(lat,lon,\
                                        lat_TL[None,:],lon_TL[:,None],\
                                        data[i,j,...],lon_ss,method)
        else:
            # fix this later...x
            print("(interpolate_to_TL_ndim) Error! expected 4 or fewer dimensions")
            pass

        return data_interp 
        
        

        