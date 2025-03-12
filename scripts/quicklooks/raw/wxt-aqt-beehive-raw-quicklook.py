import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl
import sage_data_client
import datetime
import pandas as pd
import xarray as xr
import numpy as np

from metpy.calc import dewpoint_from_relative_humidity, wet_bulb_temperature
from metpy.units import units
from zoneinfo import ZoneInfo
from scipy.stats.mstats import pearsonr

pd.options.mode.copy_on_write = True

def convert_mpl_times(ticks):
    """
    Convert matplotlib tickmarks (which are timezone naive) to account for local timezone (Chicago)

    Input
    -----
    ticks : list (floats)
        List containing floats of MPL numpy Datatime64 values

    Calls
    -----
    correct_timezone

    Output
    ------
    new_ticks : list (floats)
        List containing floats of MPL datetime64 values (subtracting out 
        the seconds to account for the timezone
    """
    def correct_timezone(dobj):
        """
        Set the datetime object to our Chicago time, check the utc offset. 

        Subtract the hours from the UTC datetime object and return
        """
        if dobj.astimezone(ZoneInfo("America/Chicago")).strftime("%z") == "-0600":
            out_date = dobj - datetime.timedelta(hours=6)
        elif dobj.astimezone(ZoneInfo("America/Chicago")).strftime("%z") == "-0500":
            out_date = dobj - datetime.timedelta(hours=5)
        else:
            out_date = dobj

        return out_date
        
    # convert np datetime64 to datetime objects with mpl
    utc_ticks = [mpl.dates.num2date(x) for x in ticks]
    # remove the timezone from datetime objects
    local_ticks = [correct_timezone(x) for x in utc_ticks]
    # convert the datetime objects back to np datetime64 ticks
    correct_ticks = [mpl.dates.date2num(x) for x in local_ticks]
    
    return correct_ticks

def timeseries(args, wxt_ds, aqt_ds, DATE=None):
    """
    With resampled and matched WXT/AQT DataSets, create a time series quicklook
    that displays temperature, wind speed/direction, particulate matter and 
    gaseous compositions for the selected date. 

    Input
    -----
    wxt_ds : Xarray DataSet
        Vaisala WXT data resampled to 1 min frequency for the selected dates
    aqt_ds : Xarray DataSet
        Vaisala AQT data resampled to 1 min frequency for the selected dates

    Keywords
    --------
    DATE : str
        Desired Date to display and use to select data

    Returns
    -------
    axarr : Matplotlib.pyplot Figure 
        Figure instance containing the desired quicklook; to be saved outside 
        of function
    """

    fig, axarr = plt.subplots(4, 1, figsize=[10, 14])
    fig.subplots_adjust(hspace=0.4)

    # --------------------------
    # Air & Dewpoint Temperature
    # --------------------------
    wxt_ds.temperature.plot(ax=axarr[0], label="WXT-Temp")
    aqt_ds.temperature.plot(ax=axarr[0], label="AQT-Temp")

    wxt_ds.dewpoint.plot(ax=axarr[0], label="WXT-Dewpoint", color="red")
    aqt_ds.dewpoint.plot(ax=axarr[0], label="AQT-Dewpoint", color='green')
    
    if DATE is not None:
        newDate = DATE[0:4] + '/' + DATE[4:6] + '/' + DATE[6:8]
    else:
        newDate = args.start_date.split('T')[0] + ' - ' + args.end_date.split('T')[0]

    axarr[0].set_title(args.node + ' RAW/BEEHIVE - ' + newDate)
    axarr[0].set_ylabel("Temperature \n [C]")
    axarr[0].set_xlabel("Time [UTC]")
    axarr[0].legend(loc="upper left")
    axarr[0].grid(True)
    axarr[0].xaxis.set_major_formatter(mpl.dates.DateFormatter('%Y-%m-%d \n %H:%M:%S'))

    # -------------------------
    # Wind Direction and Speed
    # -------------------------
    ax2 = axarr[1].twinx()
    wxt_ds.wind_dir_10s.plot(ax=axarr[1], 
                            label='Wind Direction',
                            color="tab:blue"
    )
    wxt_ds.wind_mean_10s.plot(ax=ax2, label="Wind Speed", color="tab:orange")

    axarr[1].set_ylim([0, 360])
    axarr[1].grid(True)
    axarr[1].set_xlabel("Time [UTC]")
    axarr[1].set_ylabel(r'Wind Direction [$^\circ$]', color="tab:blue")
    ax2.set_ylabel(r"Wind Speed [m/s]", color="tab:orange")
    axarr[1].xaxis.set_major_formatter(mpl.dates.DateFormatter('%Y-%m-%d \n %H:%M:%S'))

    # -----------------------
    # Particle Concentrations
    # -----------------------
    ax3 = axarr[2].twinx()
    aqt_ds["pm1.0"].plot(ax=axarr[2], label="PM 1.0")
    aqt_ds["pm2.5"].plot(ax=axarr[2], label="PM 2.5")
    aqt_ds["pm10.0"].plot(ax=axarr[2], label="PM 10")

    axarr[2].set_ylabel(r"Conc [ug/$m^3$]")
    axarr[2].set_xlabel("Time [UTC]")
    axarr[2].legend(loc="upper left")
    axarr[2].grid(True)
    wxt_ds.rainfall.plot(ax=ax3, color="tab:red")
    ax3.set_ylabel(r"Rainfall Accumulation [mm]", color="tab:red")
    axarr[2].xaxis.set_major_formatter(mpl.dates.DateFormatter('%Y-%m-%d \n %H:%M:%S'))
                                                      
    # -----------------------
    # Gas Concentrations
    # -----------------------
    aqt_ds["no2"].plot(ax=axarr[3], label="Nitrogen Dioxide")
    aqt_ds["o3"].plot(ax=axarr[3], label="Ozone")
    aqt_ds["no"].plot(ax=axarr[3], label="Nitrogen Oxide")
    new_co = aqt_ds["co"].data / 10
    axarr[3].plot(aqt_ds["time"], new_co, label="Carbon Monoxide (/10)")

    axarr[3].set_ylabel(r"Conc [ppm]")
    axarr[3].set_xlabel("Time [UTC]")
    axarr[3].legend(loc="upper left")
    axarr[3].grid(True)
    axarr[3].xaxis.set_major_formatter(mpl.dates.DateFormatter('%Y-%m-%d \n %H:%M:%S'))

    # Define the twin axis. 
    axB = axarr[3].twiny()
    # Define the tick locations. 
    new_tick_locations = []
    for tick in axarr[3].get_xticks():
        new_tick_locations.append(tick)
    # Define the UTC date in Local Time
    newtimes = convert_mpl_times(new_tick_locations)
    # Move twinned axis ticks and label from top to bottom
    axB.xaxis.set_ticks_position("bottom")
    axB.xaxis.set_label_position("bottom")
    # Offset the twin axis below the host. 
    axB.spines["bottom"].set_position(("axes",-0.4))
    # Turn on only the bottom spine. 
    axB.set_frame_on(True)
    axB.patch.set_visible(False)
    # Set the locations. 
    axB.xaxis.set_major_locator(mpl.dates.AutoDateLocator())
    axB.xaxis.set_major_formatter(mpl.dates.DateFormatter('%Y-%m-%d \n %H:%M:%S'))
    axB.set_xticks(newtimes[:])
    # Define the dual axis labels. 
    axB.set_xlabel("Time [US/Chicago]")
    # Define the dual axis limits. 
    axB.set_xlim(min(newtimes),max(newtimes))
    # Set the limits to each subplot
    axarr[0].set_xlim(min(new_tick_locations), max(new_tick_locations))
    axarr[1].set_xlim(min(new_tick_locations), max(new_tick_locations))
    axarr[2].set_xlim(min(new_tick_locations), max(new_tick_locations))
    axarr[3].set_xlim(min(new_tick_locations), max(new_tick_locations))
    
    return axarr

def environment(args, wxt_ds, aqt_ds, DATE=None):
    fig, axarr2 = plt.subplots(3, 2, figsize=[18, 14])
    fig.subplots_adjust(hspace=0.25, wspace=0.1)

    if DATE is not None:
        newDate = DATE[0:4] + '/' + DATE[4:6] + '/' + DATE[6:8]
    else:
        newDate = args.start_date.split('T')[0] + ' - ' + args.end_date.split('T')[0]

    plt.suptitle(args.node + 'RAW/BEEHIVE - ' + newDate)

    # ----------------
    # Average the Data
    # ----------------
    wxt_1min_ds = wxt_ds.resample(time="1Min").mean()
    aqt_1min_ds = aqt_ds.resample(time="1Min").mean()

    # --------------------------
    # Air Temp / Wind Direction
    # --------------------------
    var_min = np.min([wxt_1min_ds.sel(time=slice(DATE, DATE)).temperature.min(skipna=True).data,
                     aqt_1min_ds.sel(time=slice(DATE, DATE)).temperature.min(skipna=True).data]
    )
    var_max = np.max([np.round(wxt_1min_ds.sel(time=slice(DATE, DATE)).temperature.max(skipna=True).data),
                     np.round(aqt_1min_ds.sel(time=slice(DATE, DATE)).temperature.max(skipna=True).data)]
    )

    scc = axarr2[0, 0].scatter(wxt_1min_ds.sel(time=slice(DATE, DATE)).temperature.data,
                               aqt_1min_ds.sel(time=slice(DATE, DATE)).temperature.data,
                               c=wxt_1min_ds.sel(time=slice(DATE, DATE)).wind_dir_10s.data,
                               cmap="coolwarm"
    )

    axarr2[0, 0].set_xlim([var_min - 1, var_max + 1])
    axarr2[0, 0].set_ylim([var_min - 1, var_max + 1])
    axarr2[0, 0].set_title('Ambient Air Temperature [C]')
    axarr2[0, 0].set_xlabel('WXT')
    axarr2[0, 0].set_ylabel('AQT')
    axarr2[0, 0].grid(True)
    # determine 1:1 ratio line
    ratio = np.linspace(var_min-1, var_max+1)
    axarr2[0, 0].plot(ratio, ratio, color="tab:red", linestyle="--", linewidth=2.0)
    # plot the colorbar
    cbar = plt.colorbar(scc)
    cbar.ax.set_ylabel(r"WXT Wind Direction [$^o$]")

    # Determine the best fit line
    aqt = np.ma.masked_invalid(aqt_1min_ds.sel(time=slice(DATE, DATE)).temperature.data)
    wxt = np.ma.masked_invalid(wxt_1min_ds.sel(time=slice(DATE, DATE)).temperature.data)
    z = np.ma.polyfit(wxt, aqt, 1)
    p = np.poly1d(z)
    
    # Plot the best fit line
    axarr2[0, 0].plot(wxt, p(wxt), 'k', linewidth=2)
    # Display the line equation
    axarr2[0, 0].text(var_min-0.75, var_max, f"y = {z[0]:.3f}x + ({z[1]:.3f})", color='k', fontsize=12)
    # Calculate Pearson Correlation Coefficient
    cc_conc = pearsonr(wxt, aqt)
    # Display the Pearson CC
    axarr2[0, 0].text(var_min-0.75, var_max-1.0, "Pearson CC: %.2f" % (cc_conc[0]), fontsize=12)
    # Display the total number of samples
    axarr2[0, 0].text(var_min-0.75, var_max-2.0, "N = %.0f" % (wxt_1min_ds.sel(time=slice(DATE, DATE)).temperature.data.shape[0]), fontsize=12)

    # ---------------------
    # Air Temp / Wind Speed
    # ---------------------
    scc = axarr2[0, 1].scatter(wxt_1min_ds.sel(time=slice(DATE, DATE)).temperature.data,
                               aqt_1min_ds.sel(time=slice(DATE, DATE)).temperature.data,
                               c=wxt_1min_ds.sel(time=slice(DATE, DATE)).wind_mean_10s.data,
                               cmap="coolwarm"
    )

    axarr2[0, 1].set_xlim([var_min - 1, var_max + 1])
    axarr2[0, 1].set_ylim([var_min - 1, var_max + 1])
    axarr2[0, 1].set_title('Ambient Air Temperature [C]')
    axarr2[0, 1].set_xlabel('WXT')
    axarr2[0, 1].set_ylabel('AQT')
    axarr2[0, 1].grid(True)
    # determine 1:1 ratio line
    ratio = np.linspace(var_min-1, var_max+1)
    axarr2[0, 1].plot(ratio, ratio, color="tab:red", linestyle="--", linewidth=2.0)
    # plot the colorbar
    cbar = plt.colorbar(scc)
    cbar.ax.set_ylabel(r"WXT Mean Wind Speed [m/s]")

    # Plot the best fit line
    axarr2[0, 1].plot(wxt, p(wxt), 'k', linewidth=2)
    # Display the line equation
    axarr2[0, 1].text(var_min-0.75, var_max+0.0, f"y = {z[0]:.3f}x + ({z[1]:.3f})", color='k', fontsize=12)
    # Calculate Pearson Correlation Coefficient
    cc_conc = pearsonr(wxt, aqt)
    # Display the Pearson CC
    axarr2[0, 1].text(var_min-0.75, var_max-1.0, "Pearson CC: %.2f" % (cc_conc[0]), fontsize=12)
    # Display the total number of samples
    axarr2[0, 1].text(var_min-0.75, var_max-2.0, "N = %.0f" % (wxt_1min_ds.sel(time=slice(DATE, DATE)).temperature.data.shape[0]), fontsize=12)
  

    # --------------------------
    # Humidity / Wind Direction
    # --------------------------
    var_min = np.min([wxt_1min_ds.sel(time=slice(DATE, DATE)).humidity.min(skipna=True).data,
                     aqt_1min_ds.sel(time=slice(DATE, DATE)).humidity.min(skipna=True).data]
    )
    var_max = np.max([np.round(wxt_1min_ds.sel(time=slice(DATE, DATE)).humidity.max(skipna=True).data),
                     np.round(aqt_1min_ds.sel(time=slice(DATE, DATE)).humidity.max(skipna=True).data)]
    )

    scc = axarr2[1, 0].scatter(wxt_1min_ds.sel(time=slice(DATE, DATE)).humidity.data,
                               aqt_1min_ds.sel(time=slice(DATE, DATE)).humidity.data,
                               c=wxt_1min_ds.sel(time=slice(DATE, DATE)).wind_dir_10s.data,
                               cmap="coolwarm"
    )

    axarr2[1, 0].set_xlim([var_min - 1, 100])
    axarr2[1, 0].set_ylim([var_min - 1, 100])
    axarr2[1, 0].set_title('Humidity [%]')
    axarr2[1, 0].set_xlabel('WXT')
    axarr2[1, 0].set_ylabel('AQT')
    axarr2[1, 0].grid(True)
    # determine 1:1 ratio line
    ratio = np.linspace(var_min-1, 100)
    axarr2[1, 0].plot(ratio, ratio, color="tab:red", linestyle="--", linewidth=2.0)
    # plot the colorbar
    cbar = plt.colorbar(scc)
    cbar.ax.set_ylabel(r"WXT Wind Direction [$^o$]")

    # Determine the best fit line
    aqt = np.ma.masked_invalid(aqt_1min_ds.sel(time=slice(DATE, DATE)).dewpoint.data)
    wxt = np.ma.masked_invalid(wxt_1min_ds.sel(time=slice(DATE, DATE)).dewpoint.data)
    z = np.ma.polyfit(wxt, aqt, 1)
    p = np.poly1d(z)
    # Plot the best fit line
    axarr2[1, 0].plot(wxt, p(wxt), 'k', linewidth=2)
    # Display the line equation
    axarr2[1, 0].text(var_min-0.25, 96, f"y = {z[0]:.3f}x + ({z[1]:.3f})", color='k', fontsize=12)
    # Calculate Pearson Correlation Coefficient
    cc_conc = pearsonr(wxt, aqt)
    # Display the Pearson CC
    axarr2[1, 0].text(var_min-0.25, 93, "Pearson CC: %.2f" % (cc_conc[0]), fontsize=12)
    # Display the total number of samples
    axarr2[1, 0].text(var_min-0.25, 90, "N = %.0f" % (wxt_1min_ds.sel(time=slice(DATE, DATE)).dewpoint.data.shape[0]), fontsize=12)

    # ----------------------
    # Humidity / Wind Speed
    # ----------------------
    scc = axarr2[1, 1].scatter(wxt_1min_ds.sel(time=slice(DATE, DATE)).humidity.data,
                               aqt_1min_ds.sel(time=slice(DATE, DATE)).humidity.data,
                               c=wxt_1min_ds.sel(time=slice(DATE, DATE)).wind_mean_10s.data,
                               cmap="coolwarm"
    )

    axarr2[1, 1].set_xlim([var_min - 1, 100])
    axarr2[1, 1].set_ylim([var_min - 1, 100])
    axarr2[1, 1].set_title('Humidity [%]')
    axarr2[1, 1].set_xlabel('WXT')
    axarr2[1, 1].set_ylabel('AQT')
    axarr2[1, 1].grid(True)
    # determine 1:1 ratio line
    ratio = np.linspace(var_min-1 , 100)
    axarr2[1, 1].plot(ratio, ratio, color="tab:red", linestyle="--", linewidth=2.0)
    # plot the colorbar
    cbar = plt.colorbar(scc)
    cbar.ax.set_ylabel(r"WXT Mean Wind Speed [m/s]")

    # Plot the best fit line
    axarr2[1, 1].plot(wxt, p(wxt), 'k', linewidth=2)
    # Display the line equation
    axarr2[1, 1].text(var_min-0.25, 96, f"y = {z[0]:.3f}x + ({z[1]:.3f})", color='k', fontsize=12)
    # Calculate Pearson Correlation Coefficient
    cc_conc = pearsonr(wxt, aqt)
    # Display the Pearson CC
    axarr2[1, 1].text(var_min-0.25, 93, "Pearson CC: %.2f" % (cc_conc[0]), fontsize=12)
    # Display the total number of samples
    axarr2[1, 1].text(var_min-0.25, 90, "N = %.0f" % (wxt_1min_ds.sel(time=slice(DATE, DATE)).dewpoint.data.shape[0]), fontsize=12)

    # --------------------------
    # Pressure / Wind Direction
    # --------------------------
    var_min = np.min([wxt_1min_ds.sel(time=slice(DATE, DATE)).pressure.min(skipna=True).data,
                     aqt_1min_ds.sel(time=slice(DATE, DATE)).pressure.min(skipna=True).data]
    )
    var_max = np.max([np.round(wxt_1min_ds.sel(time=slice(DATE, DATE)).pressure.max(skipna=True).data),
                     np.round(aqt_1min_ds.sel(time=slice(DATE, DATE)).pressure.max(skipna=True).data)]
    )

    scc = axarr2[2, 0].scatter(wxt_1min_ds.sel(time=slice(DATE, DATE)).pressure.data,
                               aqt_1min_ds.sel(time=slice(DATE, DATE)).pressure.data,
                               c=wxt_1min_ds.sel(time=slice(DATE, DATE)).wind_dir_10s.data,
                               cmap="coolwarm"
    )

    axarr2[2, 0].set_xlim([var_min - 1, var_max + 1])
    axarr2[2, 0].set_ylim([var_min - 1, var_max + 1])
    axarr2[2, 0].set_title('Ambient Pressure [hPa]')
    axarr2[2, 0].set_xlabel('WXT')
    axarr2[2, 0].set_ylabel('AQT')
    axarr2[2, 0].grid(True)
    # determine 1:1 ratio line
    ratio = np.linspace(var_min-1, var_max+1)
    axarr2[2, 0].plot(ratio, ratio, color="tab:red", linestyle="--", linewidth=2.0)
    # plot the colorbar
    cbar = plt.colorbar(scc)
    cbar.ax.set_ylabel(r"WXT Wind Direction [$^o$]")

    # Determine the best fit line
    aqt = np.ma.masked_invalid(aqt_1min_ds.sel(time=slice(DATE, DATE)).pressure.data)
    wxt = np.ma.masked_invalid(wxt_1min_ds.sel(time=slice(DATE, DATE)).pressure.data)
    z = np.ma.polyfit(wxt, aqt, 1)
    p = np.poly1d(z)
    # Plot the best fit line
    axarr2[2, 0].plot(wxt, p(wxt), 'k', linewidth=2)
    # Display the line equation
    axarr2[2, 0].text(var_min-0.9, var_max+0.6, f"y = {z[0]:.3f}x + ({z[1]:.3f})", color='k', fontsize=12)
    # Calculate Pearson Correlation Coefficient
    cc_conc = pearsonr(wxt, aqt)
    # Display the Pearson CC
    axarr2[2, 0].text(var_min-0.9, var_max+0.2, "Pearson CC: %.2f" % (cc_conc[0]), fontsize=12)
    # Display the total number of samples
    axarr2[2, 0].text(var_min-0.9, var_max-0.2, "N = %.0f" % (wxt_1min_ds.sel(time=slice(DATE, DATE)).pressure.data.shape[0]), fontsize=12)

    # ---------------------
    # Pressure / Wind Speed
    # ---------------------

    scc = axarr2[2, 1].scatter(wxt_1min_ds.sel(time=slice(DATE, DATE)).pressure.data,
                               aqt_1min_ds.sel(time=slice(DATE, DATE)).pressure.data,
                               c=wxt_1min_ds.sel(time=slice(DATE, DATE)).wind_mean_10s.data,
                               cmap="coolwarm"
    )

    axarr2[2, 1].set_xlim([var_min - 1, var_max + 1])
    axarr2[2, 1].set_ylim([var_min - 1, var_max + 1])
    axarr2[2, 1].set_title('Ambient Pressure [hPa]')
    axarr2[2, 1].set_xlabel('WXT')
    axarr2[2, 1].set_ylabel('AQT')
    axarr2[2, 1].grid(True)
    # determine 1:1 ratio line
    ratio = np.linspace(var_min-1, var_max+1)
    axarr2[2, 1].plot(ratio, ratio, color="tab:red", linestyle="--", linewidth=2.0)
    # plot the colorbar
    cbar = plt.colorbar(scc)
    cbar.ax.set_ylabel(r"WXT Mean Wind Speed [m/s]")

    # Plot the best fit line
    axarr2[2, 1].plot(wxt, p(wxt), 'k', linewidth=2)
    # Display the line equation
    axarr2[2, 1].text(var_min-0.9, var_max+0.6, f"y = {z[0]:.3f}x + ({z[1]:.3f})", color='k', fontsize=12)
    # Calculate Pearson Correlation Coefficient
    cc_conc = pearsonr(wxt, aqt)
    # Display the Pearson CC
    axarr2[2, 1].text(var_min-0.9, var_max+0.2, "Pearson CC: %.2f" % (cc_conc[0]), fontsize=12)
    # Display the total number of samples
    axarr2[2, 1].text(var_min-0.9, var_max-0.2, "N = %.0f" % (wxt_1min_ds.sel(time=slice(DATE, DATE)).pressure.data.shape[0]), fontsize=12)

    # free up memory
    del wxt_1min_ds, aqt_1min_ds, z, p, aqt, wxt, cc_conc, ratio

    return axarr2

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
    print("\nWXT/AQT Beehive Daily Quicklooks\n")

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
    
    # Display the instrument statistics for cron 
    print("\nNumber of Obs, Expected Number of Obs, Percent of Expected")
    if wxt_ds:
        print("Number of WXT observations: ", 
              len(wxt_ds.temperature.data),
              ", 8640 ,",
              (len(wxt_ds.temperature.data) / 8640.)*100,
              " %"
        )
    else:
        print("Number of WXT observations: No Data Found")

    if aqt_ds:
        print("Number of AQT observations: ", 
              len(aqt_ds.temperature.data),
              ", 4320 ,",
              (len(aqt_ds.temperature.data) / 4320.)*100,
              " %\n"
        )
    else:
        print("Number of AQT observations: No Data Found")

    # Display the daily mean values as santity check
    print("Daily Mean - WXT Observations")
    for var in wxt_ds.data_vars:
        print(var + ' - ' + str(wxt_ds[var].resample(time="1D").mean().data[0].round(decimals=3)))

    print("\nDaily Mean - AQT Observations")
    for var in aqt_ds.data_vars:
        print(var + ' - ' + str(aqt_ds[var].resample(time="1D").mean().data[0].round(decimals=3)))


    # define the output filename
    nout = (site_args['site_ID'] + '_' + site_args['WSN'] + '_' + 
            args.start_date.split('T')[0] + '_' + args.start_date.split('T')[1][:-1].replace(':', '') +
            '_' +
            args.end_date.split('T')[0] + '_' + args.end_date.split('T')[1][:-1].replace(':', '')
    )
    # Plot the data
    time_fig = timeseries(args, wxt_ds, aqt_ds)
    plt.savefig(args.outdir + nout + '_timeseries.png')
    
    env_fig = environment(args, wxt_ds, aqt_ds)
    plt.savefig(args.outdir + nout + '_environment.png')

    # clean up memory
    del wxt_ds, aqt_ds, time_fig, env_fig

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
    
    #put these in a dictionary for accessing
    global_sites = {'NU' : global_NU, 
                    'CSU': global_CSU,
                    'NEIU': global_NEIU,
                    'UIC': global_UIC,
                    'ATMOS': global_ATMOS,
                    'ADM': global_ADM,
                    'CCICS': global_CCICS,
                    'BIG': global_BIG,
                    'HUM': global_HUM,
                    "DOWN": global_DOWN,
                    "SHEDD": global_SHEDD}
    
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
                                 datetime.timedelta(days=1)).strftime('%Y-%m-%dT%H:%M:00Z'),
                        help="Date to Start Quicklook Creation in YYYY-MM-DDThh:mm:ssZ format"
    )
    parser.add_argument("--end_date",
                        type=str,
                        dest="end_date",
                        default=datetime.datetime.now(ZoneInfo(waggle_timezone)).strftime('%Y-%m-%dT%H:%M:00Z'),
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
