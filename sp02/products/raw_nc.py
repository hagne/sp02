#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 11:18:16 2022

@author: hagen

Hier should be everthing related to the raw netcdf file. In particulare the
instances that take the raw and make more out of it
"""
import pandas as pd
import numpy as np
import sp02.calibration as calib
import atmPy.general.measurement_site as atmms

class SP02RawData(object):
    def __init__(self, dataset, site = None, langley_fit_settings = None):
        assert(False), 'deprecated, this is now in atmPy.radiation.observations.spectral_irradiance'
        self.raw_data = dataset
        if isinstance(site, type(None)):
            assert('site' in dataset.attrs.keys()), 'If site is None, then the dataset has have lat,lon,site, site_name, attributes'
            self.site = atmms.Station(lat= dataset.attrs['site_latitude'], 
                                      lon = dataset.attrs['site_longitude'], 
                                      alt = dataset.attrs['site_elevation'], 
                                      name = dataset.attrs['site_name'], 
                                      abbreviation = dataset.attrs['site'],)
        else:
            self.site = site
        self.langley_fit_settings = langley_fit_settings
        self._sun_position = None
        self._am = None
        self._pm = None
        self._transmission = None
        # self._langleys_am = None
        # self._langleys_pm = None
        # self._langley_fitres_am = None
        # self._langley_fitres_pm = None
    
    @property
    def transmission(self):
        if isinstance(self._transmission, type(None)):
            #### load calibrations
            calibrations = calib.load_calibration_history()
            cal = calibrations[int(self.raw_data.serial_no.values)]
            # use the mean and only the actual channels, other channels are artefacts
            cal = cal.results['mean'].loc[:,self.raw_data.channle_wavelengths.values].sort_index()
            
            #### interpolate and resample calibration (V0s)
            dt = self.raw_data.datetime.to_pandas()
            calib_interp = pd.concat([cal,dt]).drop([0], axis = 1).sort_index().interpolate().reindex(dt.index)
            
            #### correct VOs for earth sun distance see functions above
            calib_interp_secorr = calib_interp.divide(self.sun_position.sun_earth_distance**2, axis = 0)
            
            #### match channels for operation
            channels = self.raw_data.channle_wavelengths.to_pandas()
            raw_data = self.raw_data.raw_data.to_pandas().rename(columns = channels)
            raw_data.columns.name = 'wl'
            
            #### get transmission
            self._transmission = raw_data/calib_interp_secorr
        return self._transmission
    
    @property
    def sun_position(self):
        if isinstance(self._sun_position, type(None)):
            self._sun_position = self.site.get_sun_position(self.raw_data.datetime)
        return self._sun_position
    
    @property
    def am(self):
        if isinstance(self._am, type(None)):
            self._get_langley_from_raw() 
        return self._am
    
    @property
    def pm(self):
        if isinstance(self._pm, type(None)):
            self._get_langley_from_raw() 
        return self._pm
    
    def tp_get_rdl(self):
        raw_df = self.raw_data.raw_data.to_pandas()
        
        # changing to local time
        raw_df_loc = raw_df.copy()
        index_local = raw_df.index + pd.to_timedelta(self.site.time_zone[1], 'h')
        raw_df_loc.index = index_local
        self.raw_df_loc = raw_df_loc
        
    
    def _get_langley_from_raw(self):
        raw_df = self.raw_data.raw_data.to_pandas()
        
        #### changing to local time
        raw_df_loc = raw_df.copy()
        index_local = raw_df.index + pd.to_timedelta(self.site.time_zone['diff2UTC_of_standard_time'], 'h')
        raw_df_loc.index = index_local
        # self.tp_rdl = raw_df_loc.copy()
        
        ##### getting the one day
        sunpos = self.sun_position.copy()
        start = raw_df_loc.index[0]
        if sunpos.iloc[0].airmass > 0:
            start = pd.to_datetime(f'{start.year}{start.month:02d}{start.day:02d}') + pd.to_timedelta(1,'d')
        end = start + pd.to_timedelta(1, 'd')
        raw_df_loc = raw_df_loc.truncate(start, end)

        #### localize and cut day for sunposition
        sunpos.index = index_local
        sunpos = sunpos.truncate(start, end)

        #### remove the night
        sunpos[sunpos.airmass < 0] = np.nan

        #### get the minimum airmass befor I start cutting it out
        noon = sunpos.airmass.idxmin()

        #### normalize to the sun_earth_distance
        raw_df_loc = raw_df_loc.multiply(sunpos.sun_earth_distance**2, axis=0)
    
        # langleys are the natural logarith of the voltage over the AMF ... -> log
        # to avoid warnings and strange values do some cleaning before log
        raw_df_loc[raw_df_loc <= 0] = np.nan
#         self.tp_raw_df = raw_df.copy()
        raw_df_loc = np.log(raw_df_loc)    
    
        # keep only what is considered relevant airmasses
        amf_min = 2.2 
        amf_max = 4.7
        sunpos[sunpos.airmass < amf_min] = np.nan
        sunpos[sunpos.airmass > amf_max] = np.nan

        sunpos_am = sunpos.copy()
        sunpos_pm = sunpos.copy()

        sunpos_am[sunpos.index > noon] = np.nan
        sunpos_pm[sunpos.index < noon] = np.nan
        

        langley_am = raw_df_loc.copy()
        langley_pm = raw_df_loc.copy()

        self.tp_sp_am = sunpos_am
        self.tp_sp_pm = sunpos_pm
        self.tp_df_am = langley_am[~sunpos_am.airmass.isna()].copy()
        self.tp_df_pm = langley_am[~sunpos_pm.airmass.isna()].copy()

        langley_am.index = sunpos_am.airmass
        langley_pm.index = sunpos_pm.airmass

        self._am = calib.Langley(self,langley_am[~langley_am.index.isna()], langley_fit_settings = self.langley_fit_settings)
        self._pm = calib.Langley(self,langley_pm[~langley_pm.index.isna()], langley_fit_settings = self.langley_fit_settings)
        return True
    
    