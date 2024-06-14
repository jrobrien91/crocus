import argparse
import matplotlib.pyplot as plt

def timeseries(wxt_ds, aqt_ds, DATE=None):
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
    fig.subplots_adjust(hspace=0.5)

    # --------------------------
    # Air & Dewpoint Temperature
    # --------------------------
    wxt_ds.temperature.plot(ax=axarr[0], label="WXT-Temp")
    aqt_ds.temperature.plot(ax=axarr[0], label="AQT-Temp")

    wxt_ds.dewpoint.plot(ax=axarr[0], label="WXT-Dewpoint", color="red")
    aqt_ds.dewpoint.plot(ax=axarr[0], label="AQT-Dewpoint", color='green')

    newDate = DATE[0:4] + '/' + DATE[4:6] + '/' + DATE[6:8]

    axarr[0].set_title("NEIU - W08D - " + newDate)
    axarr[0].set_ylabel("Temperature \n [C]")
    axarr[0].set_xlabel("Time [UTC]")
    axarr[0].legend(loc="upper left")
    axarr[0].grid(True)

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
    axarr[1].set_ylabel(r'Wind Direction [$^\circ$]', color="tab:blue")
    ax2.set_ylabel(r"Wind Speed [m/s]", color="tab:orange")

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
    axarr[3].legend(loc="lower left")
    axarr[3].grid(True)

    #plt.savefig('NEIU_W08D_' + DATE + '_timeseries.png')
    return axarr

def main(args):
     print('stuff')
     print(args.node)

if __name__ == '__main__':
     descript = ("Generation of WXT/AQT QC Quicklooks for Selected CROCUS Node" +
                 " utilizing ingested (b1-level) files"
     )
     parser = argparse.ArgumentParser(description=descript)
 
     parser.add_argument("--node",
                         type=str,
                         dest='node',
                         default="W08D",
                         help="Specific CROCUS node to display data from"
                         )
     parser.add_argument("--start_date",
                         type=str,
                         dest='start_date',
                         default="20230601",
                         help="Date to Start Quicklook Creation"
                         )
     parser.add_argument("--end_date",
                         type=str,
                         dest="end_date",
                         default="20240401",
                         help="Date to End Quicklook Creation"
                         )

     args = parser.parse_args()

     main(args)