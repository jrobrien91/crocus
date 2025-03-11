"""
crocus_mrms_3d_ref.py

Using Amazon S3 storage of MRMS, download all 33 level files, merge, and display

HISTORY:
    10 March 2025 - Joe O'Brien <obrienj@anl.gov> - Written and tested with 
        4 March 2025 case. 
"""

import glob
import tempfile
import gzip
import os
import argparse
import datetime
from zoneinfo import ZoneInfo

import xarray as xr
import numpy as np

from cartopy import crs as ccrs, feature as cfeature
from cartopy.io.img_tiles import OSM
from matplotlib.transforms import offset_copy
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

import cfgrib
import cmweather
import fsspec
import imageio

def mrms_ref_mosaic(ds_merged,
                    slice_sites,
                    chi_box=(271.9, 272.5, 41.6, 42.15),
                    crocus_nodes=True,
                    time_index=0,
                    elevation=1.0,
                    site="NEIU"):
    """
    Create a display to visualize MRMS 3D Reflectivity across
     the Chicago-CROCUS domain, marking locations of the nodes,
     and taking latitudial/longitudal slices. 

    Input
    -----
    ds_merged - xarray dataset
        MRMS 3-D Reflectivity, previously uncompressed and merged

    global_sites : dict
        Information that are specific to the CROCUS sites.

    chi_box : list (of floats)
        Bounding Box to display MRMS CONUS subset over in 
            [Longitude Min, Longitude Max, Latitude Min, Latitude Max]
    
    crocus_nodes : bool (default True)
        if False, do not display the CROCUS node locations        
    
    time : int
        Time indice to display within the figure

    elevation : float
        Select elevation to display within the lat/lon regional display 
    
    site : str
        Identifer of the node for lat/lon slices

    Returns
    -------
    axarr : Matplotlib.pyplot Figure 
        Figure instance containing the desired quicklook; to be saved outside 
        of function
    """

    # Initial Subset for Time
    ds_merged = ds_merged.isel(time=time_index)

    #---------------------------------------------------
    # Define the GridSpec for Detailed Subplot Placement
    #---------------------------------------------------
    fig = plt.figure(figsize=(14, 8))
    tiler = OSM()

    gs0 = gridspec.GridSpec(1, 2, figure=fig)
    gs00 = gs0[1].subgridspec(2, 1, hspace=.32)

    # update the extent of the subplot
    gs0.update(top=.90, bottom=0.1, left=0.1, right=.95)

    # -----------------------------
    # Display the Gridded Elevation
    # -----------------------------
    ax = fig.add_subplot(gs0[0], projection=ccrs.PlateCarree())

    ## subset the data
    ds_merged.sel(
        elevation=elevation, method="nearest"
    ).DBZ.plot(
        transform=ccrs.PlateCarree(),
        vmin=-15,
        vmax=65,
        ax=ax,
        cmap="ChaseSpectral",
        cbar_kwargs={"location" : "bottom"}
    )

    # Add some various map elements to the plot to make it recognizable.
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.BORDERS)
    ax.add_image(tiler, 12, zorder=1, alpha=0.55)
    ax.gridlines(draw_labels=True)

    # Set plot bounds
    ax.set_extent(chi_box)

    # add in crosshairs to indicate the lat/lon slices
    ax.axhline(y=slice_sites[site]["latitude"], color="black", linestyle="--")
    ax.axvline(x=slice_sites[site]["longitude"], color="red", linestyle="--")

    # Display the location of the CROCUS nodes
    if crocus_nodes:
        for key in slice_sites:
            # Add a marker for the CROCUS sites.
            ax.plot(slice_sites[key]['longitude'],
                    slice_sites[key]['latitude'],
                    marker='o',
                    color='black',
                    markersize=10,
                    alpha=0.7,
                    transform=ccrs.PlateCarree())

            # Use the cartopy interface to create a matplotlib transform object
            # for the Geodetic coordinate system. We will use this along with
            # matplotlib's offset_copy function to define a coordinate system which
            # translates the text by 25 pixels to the left.
            geodetic_transform = ccrs.PlateCarree()._as_mpl_transform(ax)
            text_transform = offset_copy(geodetic_transform, units='dots', x=+50, y=+15)

            # SHEDD and UIC are to close to each other. If UIC, plot to the left.
            if key == "UIC":
                # Add text 25 pixels to the left of the volcano.
                ax.text(slice_sites[key]['longitude']-0.07,
                        slice_sites[key]['latitude'],
                        slice_sites[key]['site_ID'],
                        verticalalignment='center',
                        horizontalalignment='left',
                        transform=text_transform,
                        bbox=dict(facecolor='sandybrown',
                        alpha=0.5,
                        boxstyle='round'))
            elif key == "BIG" or key == "NU" or key == "CSU":
                # Add text 25 pixels to the left of the volcano.
                ax.text(slice_sites[key]['longitude']-0.05,
                        slice_sites[key]['latitude'],
                        slice_sites[key]['site_ID'],
                        verticalalignment='center',
                        horizontalalignment='left',
                        transform=text_transform,
                        bbox=dict(facecolor='sandybrown',
                        alpha=0.5,
                        boxstyle='round'))
            else:
                # Add text 25 pixels to the left of the volcano.
                ax.text(slice_sites[key]['longitude'],
                        slice_sites[key]['latitude'],
                        slice_sites[key]['site_ID'],
                        verticalalignment='center',
                        horizontalalignment='right',
                        transform=text_transform,
                        bbox=dict(facecolor='sandybrown',
                        alpha=0.5,
                        boxstyle='round'))

    # update the title of the display
    ax.set_title(
        np.datetime_as_string(ds_merged["valid_time"].data, unit="s")
        .replace("T", " - ") + "Z\n"
        f"{ds_merged['elevation'].data[0]} km MSL - MRMS Mosaic"
    )

    # ---------------------------
    # Display the Latitudal Slice
    # ---------------------------
    ax2 = fig.add_subplot(gs00[0])
    ds_merged.sel(
        latitude=slice_sites[site]['latitude'], method="nearest"
    ).sel(
        longitude=slice(chi_box[0] - 360, chi_box[1] - 360)
    ).DBZ.plot(
        y="elevation",
        vmin=-10,
        vmax=60,
        ax=ax2,
        cmap="ChaseSpectral"
    )
    if crocus_nodes:
        ax2.set_title(
            f"CROCUS Node - {slice_sites[site]['site_ID']}; "
            f"Latitude = {slice_sites[site]['latitude']}\N{DEGREE SIGN}"
        )
    else:
        ax2.set_title("Latitude = " + str(slice_sites[site]["latitude"]) + r"$\degree$")
    ax2.set_ylabel("Height \n [km above MSL]")
    ax2.set_xlabel(r"Longitude $\degree$")

    # ----------------------------
    # Display the Longitudal Slice
    # ----------------------------
    ax3 = fig.add_subplot(gs00[1])

    ds_merged.sel(
        longitude=slice_sites[site]['longitude'], method="nearest"
    ).sel(
        latitude=slice(chi_box[3], chi_box[2])
    ).DBZ.plot(
        y="elevation",
        vmin=-10,
        vmax=60,
        ax=ax3,
        cmap="ChaseSpectral"
    )

    if crocus_nodes:
        ax3.set_title(
            f"CROCUS Node - {slice_sites[site]['site_ID']}; "
            f"Longitude = {slice_sites[site]['longitude']}\N{DEGREE SIGN}",
            color="red"
        )
    else:
        ax3.set_title(
            f"Longitude = {slice_sites[site]['longitude']}\N{DEGREE SIGN}",
            color="red"
        )
    ax3.set_ylabel("Height \n [km above MSL]")
    ax3.set_xlabel(r"Latitude $\degree$")

    return fig

def main(nargs, in_sites):
    """
    crocus_mrms_3d_ref - main function

    INPUT
    -----
    nargs : dict
        command line defined arguments
    in_sites : dict
        CROCUS site information

    OUTPUT
    ------
    mrms - gif
        MRMS three panel display gif for user defined hour
    """

    # Initial connection to the Amazon S3 bucket
    fs = fsspec.filesystem("s3", anon=True)

    # to lazy to switch all of these
    DATE = nargs.date
    HOUR = nargs.hour
    B_BOX = nargs.bounding_box

    s3_bucket = [f"s3://noaa-mrms-pds/CONUS/MergedReflectivityQC_00.50/{DATE}/*{DATE}-{HOUR}*",
                 f"s3://noaa-mrms-pds/CONUS/MergedReflectivityQC_00.75/{DATE}/*{DATE}-{HOUR}*",
                 f"s3://noaa-mrms-pds/CONUS/MergedReflectivityQC_01.00/{DATE}/*{DATE}-{HOUR}*",
                 f"s3://noaa-mrms-pds/CONUS/MergedReflectivityQC_01.25/{DATE}/*{DATE}-{HOUR}*",
                 f"s3://noaa-mrms-pds/CONUS/MergedReflectivityQC_01.50/{DATE}/*{DATE}-{HOUR}*",
                 f"s3://noaa-mrms-pds/CONUS/MergedReflectivityQC_01.75/{DATE}/*{DATE}-{HOUR}*",
                 f"s3://noaa-mrms-pds/CONUS/MergedReflectivityQC_02.00/{DATE}/*{DATE}-{HOUR}*",
                 f"s3://noaa-mrms-pds/CONUS/MergedReflectivityQC_02.25/{DATE}/*{DATE}-{HOUR}*",
                 f"s3://noaa-mrms-pds/CONUS/MergedReflectivityQC_02.50/{DATE}/*{DATE}-{HOUR}*",
                 f"s3://noaa-mrms-pds/CONUS/MergedReflectivityQC_02.75/{DATE}/*{DATE}-{HOUR}*",
                 f"s3://noaa-mrms-pds/CONUS/MergedReflectivityQC_03.00/{DATE}/*{DATE}-{HOUR}*",
                 f"s3://noaa-mrms-pds/CONUS/MergedReflectivityQC_03.50/{DATE}/*{DATE}-{HOUR}*",
                 f"s3://noaa-mrms-pds/CONUS/MergedReflectivityQC_04.00/{DATE}/*{DATE}-{HOUR}*",
                 f"s3://noaa-mrms-pds/CONUS/MergedReflectivityQC_04.50/{DATE}/*{DATE}-{HOUR}*",
                 f"s3://noaa-mrms-pds/CONUS/MergedReflectivityQC_05.00/{DATE}/*{DATE}-{HOUR}*",
                 f"s3://noaa-mrms-pds/CONUS/MergedReflectivityQC_05.50/{DATE}/*{DATE}-{HOUR}*",
                 f"s3://noaa-mrms-pds/CONUS/MergedReflectivityQC_06.00/{DATE}/*{DATE}-{HOUR}*",
                 f"s3://noaa-mrms-pds/CONUS/MergedReflectivityQC_06.50/{DATE}/*{DATE}-{HOUR}*",
                 f"s3://noaa-mrms-pds/CONUS/MergedReflectivityQC_07.00/{DATE}/*{DATE}-{HOUR}*",
                 f"s3://noaa-mrms-pds/CONUS/MergedReflectivityQC_07.50/{DATE}/*{DATE}-{HOUR}*",
                 f"s3://noaa-mrms-pds/CONUS/MergedReflectivityQC_08.00/{DATE}/*{DATE}-{HOUR}*",
                 f"s3://noaa-mrms-pds/CONUS/MergedReflectivityQC_08.50/{DATE}/*{DATE}-{HOUR}*",
                 f"s3://noaa-mrms-pds/CONUS/MergedReflectivityQC_09.00/{DATE}/*{DATE}-{HOUR}*",
                 f"s3://noaa-mrms-pds/CONUS/MergedReflectivityQC_10.00/{DATE}/*{DATE}-{HOUR}*",
                 f"s3://noaa-mrms-pds/CONUS/MergedReflectivityQC_11.00/{DATE}/*{DATE}-{HOUR}*",
                 f"s3://noaa-mrms-pds/CONUS/MergedReflectivityQC_12.00/{DATE}/*{DATE}-{HOUR}*",
                 f"s3://noaa-mrms-pds/CONUS/MergedReflectivityQC_13.00/{DATE}/*{DATE}-{HOUR}*",
                 f"s3://noaa-mrms-pds/CONUS/MergedReflectivityQC_14.00/{DATE}/*{DATE}-{HOUR}*",
                 f"s3://noaa-mrms-pds/CONUS/MergedReflectivityQC_15.00/{DATE}/*{DATE}-{HOUR}*",
                 f"s3://noaa-mrms-pds/CONUS/MergedReflectivityQC_16.00/{DATE}/*{DATE}-{HOUR}*",
                 f"s3://noaa-mrms-pds/CONUS/MergedReflectivityQC_17.00/{DATE}/*{DATE}-{HOUR}*",
                 f"s3://noaa-mrms-pds/CONUS/MergedReflectivityQC_18.00/{DATE}/*{DATE}-{HOUR}*",
                 f"s3://noaa-mrms-pds/CONUS/MergedReflectivityQC_19.00/{DATE}/*{DATE}-{HOUR}*"
    ]

    # Download, Uncompress, Read, Convert the MRMS Data, store in list for concatenation
    ds_list = []
    for scan in s3_bucket:
        file_path = sorted(fs.glob(scan))
        elev_list = []
        for mrms in file_path:
            with fs.open(mrms, 'rb') as gzip_file:
                with tempfile.NamedTemporaryFile(suffix=".grib2") as f:
                    # Uncompress and read the file
                    f.write(gzip.decompress(gzip_file.read()))
                    ds = xr.load_dataset(f.name)
                    # To concatenate, need to add in elevation level
                    elevation_value = float(file_path[0].split('/')[-1].split('_')[-2])
                    ds = ds.assign_coords({"elevation" : elevation_value})
                    ds["elevation"].attrs["units"] = "Kilometers Above Mean Sea Level"
                    # Variables are `unknown` in these files, rename and assign units
                    ds = ds.rename({"unknown" : "DBZ"})
                    ds["DBZ"].attrs["units"] = "dBZ"
                    ds["DBZ"].attrs["long_name"] = "Reflectivity"
                    # Subset for the desired bounding box and take out all missing values
                    ds = (
                        ds.sel(latitude=slice(B_BOX[3], B_BOX[2]),
                               longitude=slice(B_BOX[0], B_BOX[1])
                        )
                        .where(ds.DBZ > -20)
                    )
                    # Transition to degrees west
                    ds['longitude'] = ds['longitude'].data - 360
                    ds['longitude'].attrs['units'] = "degrees_east"
                    # save for concatenation
                    elev_list.append(ds)
        # Concatenation per elevation
        ds_list.append(xr.concat(elev_list, dim="time"))
        # free up space
        del elev_list

    # concatenate the scans
    ds_merged = xr.concat(ds_list, dim="elevation")

    # free up more memory
    del ds_list

    # define a temporary directory to hold static images
    templocation = tempfile.mkdtemp()

    # Iterate through the files, uncompress and call the MRMS CREF function
    for i, value in enumerate(ds_merged.time.data):
        timestamp = np.datetime_as_string(ds_merged.time.data[i])
        date_part = timestamp.split('T')[0].replace('-', '')
        time_part = timestamp.split('T')[1].split('.')[0].replace(':', '')

        out_fig = mrms_ref_mosaic(ds_merged,
                                  in_sites,
                                  chi_box=tuple(B_BOX),
                                  time_index=i,
                                  crocus_nodes=nargs.override,
                                  elevation=nargs.elevation,
                                  site=nargs.node)
        out_fig.savefig(
            f"{templocation}/mrms-radar-3d-ref-{date_part}-{time_part}.png"
        )
        # free up space
        del out_fig

    # Define files created and define movie path
    map_images = sorted(glob.glob(templocation + "/mrms-radar-3d-ref*"))
    gif_title = nargs.outdir + "mrms-radar-cref-" + nargs.date + "-" + nargs.hour + ".gif"

    # Check to see if the file exists - if it does, delete it
    if os.path.exists(gif_title):
        os.remove(gif_title)

    # Loop through and create the gif
    with imageio.get_writer(gif_title, mode='I', duration=0.2) as writer:
        for filename in map_images:
            image = imageio.imread(filename)
            writer.append_data(image)

    # Free Up Memory
    del fs, ds_merged, gif_title, map_images, templocation

if __name__ == '__main__':

    WAGGLE_TIMEZONE = "UTC"
    LOCAL_TIMEZONE = "America/Chicago"

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
                  'longitude' : -87.645337665}

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
                   'latitude' : 41.795306577,
                   'longitude' : -88.006144107}

    global_SHEDD = {'conventions': "CF 1.10",
                    'WSN':'W09E',
                    'site_ID' : "SHEDD",
                    'CAMS_tag' : "CMS-WXT-007",
                    'datastream' : "CMS_wxt536_SHEDD_a1",
                    'plugin' : "registry.sagecontinuum.org/jrobrien/waggle-wxt536:0.*",
                    "aqt-plugin" : "registry.sagecontinuum.org/jrobrien/waggle-aqt:0.23.5.*",
                    'datalevel' : "a1",
                    'latitude' : 41.867953617,
                    'longitude' : -87.613603892}

    #put these in a dictionary for accessing
    global_sites = {'NU' : global_NU,
                    'CSU': global_CSU,
                    'NEIU': global_NEIU,
                    'UIC': global_UIC,
                    'ATMOS': global_ATMOS,
                    'CCICS': global_CCICS,
                    'BIG': global_BIG,
                    'HUM': global_HUM,
                    "DOWN": global_DOWN,
                    "SHEDD": global_SHEDD}

    DESCRIPT = ("Hourly Generation of Multi-Radar, Multi-Sensor (MRMS) 3-D " +
                "Reflectivity Quicklooks and GIFs for the CROCUS Chicago Domain."
    )
    BOX_DESCRIPT = (
        "Display Bounding Box in [Longitude Min, Longitude Max, "
        "Latitude Min, Latitude Max] format. "
        "Default is CROCUS-Chicago [271.9, 272.5, 41.6, 42.15]"
    )

    parser = argparse.ArgumentParser(
        description=DESCRIPT,
        usage=(
            "python crocus-mrms-3d-ref.py --date 20250304 "
            "--hour 13 --outdir /Users/jrobrien/dev/crocus/"
        )
    )

    parser.add_argument("--date",
                        type=str,
                        dest='date',
                        default=(datetime.datetime.now(ZoneInfo(WAGGLE_TIMEZONE)) -
                                 datetime.timedelta(days=1)).strftime('%Y%m%d'),
                        help="Date to Start Quicklook Creation in YYYYMMDD [UTC] format"
    )
    parser.add_argument("--hour",
                        type=str,
                        dest='hour',
                        default="12",
                        help="Hour to Start Quicklook Creation in HH [UTC] format"
    )
    parser.add_argument("--box",
                        type=list,
                        dest='bounding_box',
                        default=[271.9, 272.5, 41.6, 42.15],
                        help=BOX_DESCRIPT
    )
    parser.add_argument("--override",
                        type=bool,
                        dest='override',
                        default=True,
                        help="If False, overide display of the CROCUS node locations"
    )
    parser.add_argument("--node",
                        type=str,
                        dest='node',
                        default="ATMOS",
                        help="Specific CROCUS node to display lat/lon slices"
    )
    parser.add_argument("--elevation",
                        type=float,
                        dest='elevation',
                        default=1.0,
                        help="Default MRMS Elevation to Display"
    )
    parser.add_argument("--outdir",
                        type=str,
                        dest="outdir",
                        default="./",
                        help="Directory to output GIFs. Default is GCE home directory")

    args = parser.parse_args()

    main(args, global_sites)
