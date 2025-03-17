import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl
import sage_data_client
import datetime
import pandas as pd
import xarray as xr
import numpy as np

from zoneinfo import ZoneInfo

from metpy.calc import dewpoint_from_relative_humidity, wet_bulb_temperature
from metpy.units import units

from great_tables import GT, html, loc, exibble, style
import polars as pl
import polars.selectors as cs

pd.options.mode.copy_on_write = True

def wxt_table(ds, site_args):
    """
    Generate a research quality table to display hourly mean
    observations for the WXT

    Input
    -----
    ds : Xarray Dataset
        Xarray dataset containing ingested WXT dataset

    Output
    ------
    gt : great_tables display object
        Table display printed to screen
    """
    # Resample to hourly averages
    ds = ds.resample(time='1H').mean() 
    # Drop calculated variables from the WXT ingested data
    ds = ds.drop_vars(['wetbulb', 'humidity', 'wind_max_10s'])

    # generate table subtitle
    node_subtitle = site_args['site_ID'] + " (" + site_args['WSN'] + ") - Beehive Generated Table"

    # Great_Tables requires a dataframe, round decimals
    ##kat = ds.to_pandas().head().round(decimals=2)
    kat = ds.to_pandas().round(decimals=2).reset_index()
    # Apply utc timezone to convert to local time
    kat['local'] = kat['time'].dt.tz_localize('utc').dt.tz_convert('America/Chicago')
    # Add separate columns for date and hours for display
    kat['date'] = pd.to_datetime(kat['time']).dt.date
    kat['hour'] = pd.to_datetime(kat['time']).dt.time
    kat['chi_hour'] = pd.to_datetime(kat['local']).dt.time

    # Since we separated out full datetime, drop these
    kat = kat.drop(['time', 'local'], axis=1)

    # Convert to polar DataFrames
    pl_kat = pl.DataFrame(kat)

    # Select columns
    wind_mean = cs.starts_with("wind_mean")
    wind_dir = cs.starts_with("wind_dir")
    rain = cs.starts_with("rain")
    temp = cs.starts_with("temp")
    dew  = cs.starts_with("dew")
    pres = cs.starts_with("pres")

    # Create Great_Tables object
    gt = (
        GT(pl_kat)
        .tab_header(
            title="CROCUS Daily WXT - Hourly Mean",
            subtitle=node_subtitle
        )
        .tab_spanner(
            label="Time",
            columns=["date", "hour", "chi_hour"]
        )
        .tab_spanner(
            label="Measurement",
            columns=["temperature", "dewpoint", "pressure", "rainfall", "wind_dir_10s", "wind_mean_10s"]
        )
        .cols_move_to_start(columns=["date", "hour", "chi_hour"])
        .cols_label(
            date = html("<b>Date</b><br>YYYY-MM-DD"),
            hour = html("<b>Hour</b><br>UTC"),
            chi_hour = html("<b>Hour</b><br>Local"),
            temperature = html("<b>Temp</b><br>\N{DEGREE SIGN}C"),
            dewpoint = html("<b>Dewpt</b><br>\N{DEGREE SIGN}C"),
            pressure = html("<b>Press</b><br>hPa"),
            rainfall = html("<b>Precip</b><br>mm"),
            wind_dir_10s = html("<b>WindDir</b><br>\N{DEGREE SIGN}"),
            wind_mean_10s = html("<b>WindSpd</b><br>m/s")
        )
        # style ----
        .tab_style(
            style=style.fill(color="aliceblue"),
            locations=loc.body(columns=wind_mean),
        )
        .tab_style(
            style=style.fill(color="papayawhip"),
            locations=loc.body(columns=rain),
        )
        .tab_style(
            style=style.fill(color="mistyrose"),
            locations=loc.body(columns=temp),
        )
        .tab_style(
            style=style.fill(color="honeydew"),
            locations=loc.body(columns=dew),
        )
        .tab_style(
            style=style.fill(color="thistle"),
            locations=loc.body(columns=wind_dir)
        )
    )

    # Center align all values within the table
    gt = gt.cols_align(align="center")

    # free up memory
    del ds, kat, pl_kat

    return gt

def aqt_table(ds, site_args):
    """
    Generate a research quality table to display hourly mean
    observations for the AQT

    Input
    -----
    ds : Xarray Dataset
        Xarray dataset containing ingested AQT dataset

    Output
    ------
    gt : great_tables display object
        Table display printed to screen
    """
    # Resample to hourly averages
    ds = ds.resample(time='1H').mean() 
    # Rename the particle variable names for use in the table
    ds = ds.rename_vars({'pm2.5':"pm25", "pm1.0":"pm1", "pm10.0":"pm10"})
    # Drop calculated variables from the WXT ingested data
    ds = ds.drop_vars(['humidity'])

    # generate table subtitle
    node_subtitle = site_args['site_ID'] + " (" + site_args['WSN'] + ") - Beehive Generated Table"

    # Great_Tables requires a dataframe, round decimals
    ##kat = ds.to_pandas().head().round(decimals=2)
    kat = ds.to_pandas().round(decimals=2).reset_index()
    # Apply utc timezone to convert to local time
    kat['local'] = kat['time'].dt.tz_localize('utc').dt.tz_convert('America/Chicago')
    # Add separate columns for date and hours for display
    kat['date'] = pd.to_datetime(kat['time']).dt.date
    kat['hour'] = pd.to_datetime(kat['time']).dt.time
    kat['chi_hour'] = pd.to_datetime(kat['local']).dt.time

    # Since we separated out full datetime, drop these
    kat = kat.drop(['time', 'local'], axis=1)

    # Convert to polar DataFrames
    pl_kat = pl.DataFrame(kat)

    # Select columns
    temp = cs.starts_with("temp")
    dew  = cs.starts_with("dew")
    pm = cs.starts_with("pm")
    ozone = cs.starts_with("o3")
    pm25 = cs.starts_with("pm25")

    # Create Great_Tables object
    gt = (
        GT(pl_kat)
        .tab_header(
            title="CROCUS Daily AQT - Hourly Mean",
            subtitle=node_subtitle
        )
        .tab_spanner(
            label="Time",
            columns=["date", "hour", "chi_hour"]
        )
        .tab_spanner(
            label="Air Quality",
            columns=["pm25", "pm1", "pm10", "no", "o3", "no2", "co"]
        )
        .tab_spanner(
            label="Environment",
            columns=["temperature", "dewpoint", "pressure"]
        )
        .cols_move_to_start(columns=["date", 
                                     "hour", 
                                     "chi_hour", 
                                     "pm10", 
                                     "pm25", 
                                     "pm1"])
        .cols_label(
            date = html("<b>Date</b><br>YYYY-MM-DD"),
            hour = html("<b>Hour</b><br>UTC"),
            chi_hour = html("<b>Hour</b><br>Local"),
            temperature = html("<b>Temp</b><br>\N{DEGREE SIGN}C"),
            dewpoint = html("<b>Dewpt</b><br>\N{DEGREE SIGN}C"),
            pressure = html("<b>Press</b><br>hPa"),
            pm10 = html("<b>PM10</b><br>ug/m3"),
            pm25 = html("<b>PM2.5</b><br>ug/m3"),
            pm1 = html("<b>PM1.0</b><br>ug/m3"),
            no = html("<b>NO</b><br>ppm"),
            o3 = html("<b>O3</b><br>ppm"),
            no2 = html("<b>NO2</b><br>ppm"),
            co = html("<b>CO</b><br>ppm"),
        )
        # style ----
        .tab_style(
            style=style.fill(color="mistyrose"),
            locations=loc.body(columns=temp),
        )
        .tab_style(
            style=style.fill(color="honeydew"),
            locations=loc.body(columns=dew),
        )
        .tab_style(
            style=style.fill(color="thistle"),
            locations=loc.body(columns=pm)
        )
        .tab_style(
            style=style.fill(color="aliceblue"),
            locations=loc.body(columns=ozone),
        )
        .tab_style(
            style=style.fill(color="papayawhip"),
            locations=loc.body(columns=pm25),
        )
    )

    # Center align all values within the table
    gt = gt.cols_align(align="center")

    # free up memory
    del ds, kat, pl_kat

    return gt

def ingest_wxt(start, 
               end,
               global_attrs, 
               var_attrs):
    """
        Ingest from CROCUS WXTs using the Sage Data Client. 
    
        Parameters
        ----------
        start : str 
            String contiaining datetime in YYYY-MM-DDTHH:MM:SS in UTC

        end : str
            String containing datetime in YYYY-MM-DDTHH:MM:SS in UTC

        global_attrs : dict
            Attributes that are specific to the site.
        
        var_attrs : dict
            Attributes that map variables in Beehive to
            CF complaint netCDF valiables.
        
    
        Returns
        -------
        ds : xarray DataSet
            Xarray dataset containing xarray data
    
    """
    print("ingest_wxt: ", datetime.datetime.now(ZoneInfo(waggle_timezone)).strftime('%Y-%m-%dT%H:%M:00Z'))
    df_temp = sage_data_client.query(start=start,
                                     end=end, 
                                     filter={"name" : 'wxt.env.temp|wxt.env.humidity|wxt.env.pressure|wxt.rain.accumulation',
                                             "plugin" : global_attrs['wxt-plugin'],
                                             "vsn" : global_attrs['WSN'],
                                             "sensor" : "vaisala-wxt536"
                                     }
    )
    winds = sage_data_client.query(start=start,
                                   end=end, 
                                   filter={"name" : 'wxt.wind.speed|wxt.wind.direction',
                                           "plugin" : global_attrs['wxt-plugin'],
                                           "vsn" : global_attrs['WSN'],
                                           "sensor" : "vaisala-wxt536"
                                   }
    )
    
    hums = df_temp[df_temp['name']=='wxt.env.humidity']
    temps = df_temp[df_temp['name']=='wxt.env.temp']
    pres = df_temp[df_temp['name']=='wxt.env.pressure']
    rain = df_temp[df_temp['name']=='wxt.rain.accumulation']

    npres = len(pres)
    nhum = len(hums)
    ntemps = len(temps)
    nrains = len(rain)
    minsamps = min([nhum, ntemps, npres, nrains])

    temps['time'] = pd.DatetimeIndex(temps['timestamp'].values)

    vals = temps.set_index('time')[0:minsamps]
    vals['temperature'] = vals.value.to_numpy()[0:minsamps]
    vals['humidity'] = hums.value.to_numpy()[0:minsamps]
    vals['pressure'] = pres.value.to_numpy()[0:minsamps]
    vals['rainfall'] = rain.value.to_numpy()[0:minsamps]

    direction = winds[winds['name']=='wxt.wind.direction']
    speed = winds[winds['name']=='wxt.wind.speed']

    nspeed = len(speed)
    ndir = len(direction)
    minsamps = min([nspeed, ndir])

    speed['time'] = pd.DatetimeIndex(speed['timestamp'].values)
    windy = speed.set_index('time')[0:minsamps]
    windy['speed'] = windy.value.to_numpy()[0:minsamps]
    windy['direction'] = direction.value.to_numpy()[0:minsamps]

    winds10mean = windy.resample('10s').mean(numeric_only=True).ffill()
    winds10max = windy.resample('10s').max(numeric_only=True).ffill()
    dp = dewpoint_from_relative_humidity( vals.temperature.to_numpy() * units.degC, 
                                         vals.humidity.to_numpy() * units.percent)

    vals['dewpoint'] = dp
    vals10 = vals.resample('10s').mean(numeric_only=True).ffill() #ffil gets rid of nans due to empty resample periods
    wb = wet_bulb_temperature(vals10.pressure.to_numpy() * units.hPa,
                              vals10.temperature.to_numpy() * units.degC,
                              vals10.dewpoint.to_numpy() * units.degC)

    vals10['wetbulb'] = wb
    vals10['wind_dir_10s'] = winds10mean['direction']
    vals10['wind_mean_10s'] = winds10mean['speed']
    vals10['wind_max_10s'] = winds10max['speed']
    _ = vals10.pop('value')
     
    vals10xr = xr.Dataset.from_dataframe(vals10)
    vals10xr = vals10xr.sortby('time')
    
    vals10xr = vals10xr.assign_attrs(global_attrs)
    
    return vals10xr

def ingest_aqt(start,
               end,
               global_attrs, 
               var_attrs):
    """
        Ingest from CROCUS AQTs using the Sage Data Client. 

        Ingests a whole day of AQT data and saves it as a NetCDF to odir
    
        Parameters
        ----------
        start : str 
            String contiaining datetime in YYYY-MM-DDTHH:MM:SS in UTC

        end : str
            String containing datetime in YYYY-MM-DDTHH:MM:SS in UTC

        global_attrs : dict
            Attributes that are specific to the site.
        
        var_attrs : dict
            Attributes that map variables in Beehive to
            CF complaint netCDF valiables.
        
    
        Returns
        -------
        ds : xarray DataSet
            Xarray dataset containing xarray data
    
    """
    print("ingest_aqt: ", datetime.datetime.now(ZoneInfo(waggle_timezone)).strftime('%Y-%m-%dT%H:%M:00Z'))
    df_aq = sage_data_client.query(start=start,
                                   end=end, 
                                   filter={"plugin" : global_attrs['aqt-plugin'],
                                           "vsn" : global_attrs['WSN'],
                                           "sensor" : "vaisala-aqt530"
                                   }
    )
     # Rename specific column names
    pm25 = df_aq[df_aq['name']=='aqt.particle.pm2.5']
    pm10 = df_aq[df_aq['name']=='aqt.particle.pm1']
    pm100 = df_aq[df_aq['name']=='aqt.particle.pm10']
    no = df_aq[df_aq['name']=='aqt.gas.no']
    o3 = df_aq[df_aq['name']=='aqt.gas.ozone']
    no2 = df_aq[df_aq['name']=='aqt.gas.no2']
    co = df_aq[df_aq['name']=='aqt.gas.co']
    aqtemp = df_aq[df_aq['name']=='aqt.env.temp']
    aqhum = df_aq[df_aq['name']=='aqt.env.humidity']
    aqpres = df_aq[df_aq['name']=='aqt.env.pressure']

    # Convert instrument timestamp to Pandas Datatime object
    pm25['time'] = pd.DatetimeIndex(pm25['timestamp'].values)

    # Remove all meta data descriptions besides the index
    aqvals = pm25.loc[:, pm25.columns.intersection(["time"])]

    # Add all parameter to the output dataframe
    aqvals['pm2.5'] = pm25.value.to_numpy().astype(float)
    aqvals['pm1.0'] = pm10.value.to_numpy().astype(float)
    aqvals['pm10.0'] = pm100.value.to_numpy().astype(float)

    aqvals['no'] = no.value.to_numpy().astype(float)
    aqvals['o3'] = o3.value.to_numpy().astype(float)
    aqvals['no2'] = no2.value.to_numpy().astype(float)
    aqvals['co'] = co.value.to_numpy().astype(float)
    aqvals['temperature'] =  aqtemp.value.to_numpy().astype(float)
    aqvals['humidity'] =  aqhum.value.to_numpy().astype(float)
    aqvals['pressure'] =  aqpres.value.to_numpy().astype(float)

    # calculate dewpoint from relative humidity
    dp = dewpoint_from_relative_humidity(aqvals.temperature.to_numpy() * units.degC, 
                                         aqvals.humidity.to_numpy() * units.percent
    )
    aqvals['dewpoint'] = dp
    
    # Define the index
    aqvals = aqvals.set_index("time")
    
    valsxr = xr.Dataset.from_dataframe(aqvals)
    valsxr = valsxr.sortby('time')
    
    # Assign the global attributes
    valsxr = valsxr.assign_attrs(global_attrs)
    # Assign the individual parameter attributes
    for varname in var_attrs.keys():
        valsxr[varname] = valsxr[varname].assign_attrs(var_attrs[varname])

    return valsxr
    

def main(args, site_args, instr_attrs):

    # General information for cron
    print("\nWXT/AQT Beehive - Daily Table\n")

    # Display site and date information for cron
    print('Site:  ', site_args['site_ID'])
    print('Node:  ', site_args['WSN'])
    print('Start: ', args.start_date)
    print('End:   ', args.end_date, "\n")

    # retrieve WXT data
    wxt_ds = ingest_wxt(args.start_date, 
                        args.end_date,
                        site_args,
                        instr_attrs['wxt'])
    # retrieve the AQT data
    aqt_ds = ingest_aqt(args.start_date,
                        args.end_date,
                        site_args,
                        instr_attrs['aqt'])
    
    # define the output filename
    nout = (site_args['site_ID'] + '_' + site_args['WSN'] + '_' + 
            args.start_date.split('T')[0] + '_' + args.start_date.split('T')[1][:-1].replace(':', '') +
            '_' +
            args.end_date.split('T')[0] + '_' + args.end_date.split('T')[1][:-1].replace(':', '')
    )
    
    # Create daily tables
    #----WXT
    wxt_table_out = args.outdir + nout + "_wxt_hourly_table.png"
    gt_wxt = wxt_table(wxt_ds, site_args)
    gt_wxt.save(wxt_table_out)
    #----AQT
    aqt_table_out = args.outdir + nout + "_aqt_hourly_table.png"
    gt_aqt = aqt_table(aqt_ds, site_args)
    gt_aqt.save(aqt_table_out)

    # clean up memory
    del wxt_ds, aqt_ds, gt_wxt, gt_aqt

if __name__ == '__main__':

    waggle_timezone = "UTC"
    local_timezone = "America/Chicago"
     
    # Site attributes
    global_NEIU = {'conventions': "CF 1.10",
                   'site_ID' : "NEIU",
                   'CAMS_tag' : "CMS-WXT-002",
                   'datastream' : "CMS_wxt536_NEIU_a1",
                   'datalevel' : "a1",
                   "wxt-plugin" : "registry.sagecontinuum.org/jrobrien/waggle-wxt536:0.*",
                   "aqt-plugin" : "registry.sagecontinuum.org/jrobrien/waggle-aqt:0.23.5.*",
                   'WSN' : 'W08D',
                   'latitude' : 41.9804526,
                   'longitude' : -87.7196038}
    
    global_NU = {'conventions': "CF 1.10",
                 'WSN':'W099',
                 'site_ID' : "NU",
                 'CAMS_tag' : "CMS-WXT-005",
                 'datastream' : "CMS_wxt536_NU_a1",
                 'wxt-plugin' : "registry.sagecontinuum.org/jrobrien/waggle-wxt536:0.*",
                 "aqt-plugin" : "registry.sagecontinuum.org/jrobrien/waggle-aqt:0.23.5.*",
                 'datalevel' : "a1",
                 'latitude' : 42.051469749,
                 'longitude' : -87.677667183}
    
    global_CSU = {'conventions': "CF 1.10",
                  'WSN':'W08E',
                  'site_ID' : "CSU",
                  'CAMS_tag' : "CMS-WXT-003",
                  'datastream' : "CMS_wxt536_CSU_a1",
                  'wxt-plugin' : "registry.sagecontinuum.org/jrobrien/waggle-wxt536:0.*",
                  "aqt-plugin" : "registry.sagecontinuum.org/jrobrien/waggle-aqt:0.23.5.*",
                  'datalevel' : "a1",
                  'latitude' : 41.71991216,
                  'longitude' : -87.612834722}
    
    global_UIC = {'conventions': "CF 1.10",
                  'WSN':'W096',
                  'site_ID' : "UIC",
                  'CAMS_tag' : "CMS-WXT-006",
                  'datastream' : "CMS_wxt536_UIC_a1",
                  'wxt-plugin' : "registry.sagecontinuum.org/jrobrien/waggle-wxt536:0.*",
                  "aqt-plugin" : "registry.sagecontinuum.org/jrobrien/waggle-aqt:0.23.5.*",
                  'datalevel' : "a1",
                  'latitude' : 41.86943346,
                  'longitude' : -87.645871814}
    
    global_ATMOS = {'conventions': "CF 1.10",
                    'WSN':'W0A4',
                    'site_ID' : "ATMOS",
                    'CAMS_tag' : "CMS-WXT-001",
                    'datastream' : "CMS_wxt536_ATMOS_a1",
                    'wxt-plugin' : "registry.sagecontinuum.org/jrobrien/waggle-wxt536:0.*",
                    "aqt-plugin" : "registry.sagecontinuum.org/jrobrien/waggle-aqt:0.23.5.*",
                    'datalevel' : "a1",
                    'latitude' : 41.701556533,
                    'longitude' : -87.99507543}
    
    global_ADM = {'conventions': "CF 1.10",
                  'WSN':'W09E',
                  'site_ID' : "ADM",
                  'CAMS_tag' : "CMS-WXT-013",
                  'datastream' : "CMS_wxt536_ADM_a1",
                  'wxt-plugin' : "registry.sagecontinuum.org/jrobrien/waggle-wxt536:0.*",
                  "aqt-plugin" : "registry.sagecontinuum.org/jrobrien/waggle-aqt:0.23.5.*",
                  'datalevel' : "a1",
                  'latitude' : 41.867571129,
                  'longitude' : -87.649593711}

    global_CCICS = {'conventions': "CF 1.10",
                    'WSN':'W08B',
                    'site_ID' : "CCICS",
                    'CAMS_tag' : "CMS-WXT-099",
                    'datastream' : "CMS_wxt536_CCICS_a1",
                    'wxt-plugin' : "registry.sagecontinuum.org/jrobrien/waggle-wxt536:0.*",
                    "aqt-plugin" : "registry.sagecontinuum.org/jrobrien/waggle-aqt:0.23.5.*",
                    'datalevel' : "a1",
                    'latitude' : 41.822953674,
                    'longitude' : -87.609452418}
    
    global_BIG   = {'conventions': "CF 1.10",
                    'WSN':'W0A0',
                    'site_ID' : "BIG",
                    'CAMS_tag' : "CMS-WXT-014",
                    'datastream' : "CMS_wxt536_BIG_a1",
                    'wxt-plugin' : "registry.sagecontinuum.org/jrobrien/waggle-wxt536:0.*",
                    "aqt-plugin" : "registry.sagecontinuum.org/jrobrien/waggle-aqt:0.23.5.*",
                    'datalevel' : "a1",
                    'latitude' : 41.777009821,
                    'longitude' : -87.609746965}
        
    global_HUM  =  {'conventions': "CF 1.10",
                    'WSN':'W0A1',
                    'site_ID' : "HUM",
                    'CAMS_tag' : "CMS-WXT-010",
                    'datastream' : "CMS_wxt536_HUM_a1",
                    'wxt-plugin' : "registry.sagecontinuum.org/jrobrien/waggle-wxt536:0.*",
                    "aqt-plugin" : "registry.sagecontinuum.org/jrobrien/waggle-aqt:0.23.5.*",
                    'datalevel' : "a1",
                    'latitude' : 41.905513206,
                    'longitude' : -87.703525713}
    
    global_DOWN = {'conventions': "CF 1.10",
                    'WSN':'W09D',
                    'site_ID' : "DOWN",
                    'CAMS_tag' : "CMS-WXT-008",
                    'datastream' : "CMS_wxt536_DOWN_a1",
                    'wxt-plugin' : "registry.sagecontinuum.org/jrobrien/waggle-wxt536:0.*",
                    "aqt-plugin" : "registry.sagecontinuum.org/jrobrien/waggle-aqt:0.23.5.*",
                    'datalevel' : "a1",
                    'latitude' : 41.701476659,
                    'longitude' : -87.9953044}
    
    global_SHEDD = {'conventions': "CF 1.10",
                    'WSN':'W09E',
                    'site_ID' : "SHEDD",
                    'CAMS_tag' : "CMS-WXT-007",
                    'datastream' : "CMS_wxt536_SHEDD_a1",
                    'wxt-plugin' : "registry.sagecontinuum.org/jrobrien/waggle-wxt536:0.*",
                    "aqt-plugin" : "registry.sagecontinuum.org/jrobrien/waggle-aqt:0.23.5.*",
                    'datalevel' : "a1",
                    'latitude' : 41.867918965,
                    'longitude' : -87.613535027}

    global_VILLA = {'conventions': "CF 1.10",
                    'WSN':'W095',
                    'site_ID' : "VLPK",
                    'CAMS_tag' : "CMS-WXT-004",
                    'datastream' : "CMS_wxt536_VLPK_a1",
                    'wxt-plugin' : "registry.sagecontinuum.org/jrobrien/waggle-wxt536:0.*",
                    "aqt-plugin" : "registry.sagecontinuum.org/jrobrien/waggle-aqt:0.23.5.*",
                    'datalevel' : "a1",
                    'latitude' : 41.884884633495616,
                    'longitude' : -87.97871741056426}
    
    #put these in a dictionary for accessing
    global_sites = {'NU' : global_NU, 
                    'CSU': global_CSU,
                    'NEIU' : global_NEIU,
                    'UIC' : global_UIC,
                    'ATMOS' : global_ATMOS,
                    'ADM' : global_ADM,
                    'CCICS' : global_CCICS,
                    'BIG' : global_BIG,
                    'HUM': global_HUM,
                    "DOWN": global_DOWN,
                    "SHEDD": global_SHEDD,
                    "VLPK": global_VILLA}
    
    #Variable attributes
    var_attrs_wxt = {'temperature': {'standard_name' : 'air_temperature',
                           'units' : 'celsius'},
                    'humidity': {'standard_name' : 'relative_humidity',
                           'units' : 'percent'},
                    'dewpoint': {'standard_name' : 'dew_point_temperature',
                           'units' : 'celsius'},
                    'pressure': {'standard_name' : 'air_pressure',
                           'units' : 'hPa'},
                    'wind_mean_10s': {'standard_name' : 'wind_speed',
                           'units' : 'celsius'},
                    'wind_max_10s': {'standard_name' : 'wind_speed',
                           'units' : 'celsius'},
                    'wind_dir_10s': {'standard_name' : 'wind_from_direction',
                           'units' : 'degrees'},
                    'rainfall': {'standard_name' : 'precipitation_amount',
                           'units' : 'kg m-2'}}
    
    var_attrs_aqt = {'pm2.5' : {'standard_name' : 'mole_concentration_of_pm2p5_ambient_aerosol_particles_in_air',
                                'units' : 'ug/m^3'},
                    'pm10.0' : {'standard_name' : 'mole_concentration_of_pm10p0_ambient_aerosol_particles_in_air',
                                'units' : 'ug/m^3'},
                    'pm1.0' : {'standard_name' : 'mole_concentration_of_pm1p0_ambient_aerosol_particles_in_air',
                               'units' : 'ug/m^3'},
                    'no' : {'standard_name' : 'mole_fraction_of_nitrogen_monoxide_in_air',
                            'units' : 'Parts Per Million'},
                    'o3' : {'standard_name' : 'mole_fraction_of_ozone_in_air',
                            'units' : 'Parts Per Million'},
                    'co' : {'standard_name' : 'mole_fraction_of_carbon_monoxide_in_air',
                            'units' : 'Parts Per Million'},
                    'no2' : {'standard_name' : 'mole_fraction_of_nitrogen_dioxide_in_air',
                            'units' : 'Parts Per Million'},
                    'temperature': {'standard_name' : 'air_temperature',
                            'units' : 'celsius'},
                    'humidity': {'standard_name' : 'relative_humidity',
                            'units' : 'percent'},
                    'dewpoint': {'standard_name' : 'dew_point_temperature',
                            'units' : 'celsius'},
                    'pressure': {'standard_name' : 'air_pressure',
                            'units' : 'hPa'}}
    
    instr_attrs = {'wxt' : var_attrs_wxt,
                   'aqt' : var_attrs_aqt}


    descript = ("Generation of WXT/AQT QC Quicklooks for Selected CROCUS Node" +
                " utilizing data from Beehive (a0 level)"
    )
    parser = argparse.ArgumentParser(description=descript)
 
    parser.add_argument("--node",
                        type=str,
                        dest='node',
                        default="NEIU",
                        help="Specific CROCUS node to display data from"
    )
    parser.add_argument("--start_date",
                        type=str,
                        dest='start_date',
                        default=(datetime.datetime.now(ZoneInfo(waggle_timezone)) - 
                                 datetime.timedelta(days=1, hours=1)).strftime('%Y-%m-%dT%H:%M:00Z'),
                        help="Date to Start Quicklook Creation in YYYY-MM-DDThh:mm:ssZ format"
    )
    parser.add_argument("--end_date",
                        type=str,
                        dest="end_date",
                        default=(datetime.datetime.now(ZoneInfo(waggle_timezone)) - 
                                 datetime.timedelta(hours=1)).strftime('%Y-%m-%dT%H:%M:00Z'),
                        help="Date to End Quicklook Creation in YYYY-MM-DDThh:mm:ssZ format"
    )
    parser.add_argument("--outdir",
                        type=str,
                        dest="outdir",
                        default="./",
                        help="Directory to output figures. Default is GCE home directory")

    args = parser.parse_args()

    site_args = global_sites[args.node]

    main(args, site_args, instr_attrs)
