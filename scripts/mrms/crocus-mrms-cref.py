import cfgrib
import cmweather
import fsspec
import glob
import tempfile
import gzip
import imageio
import os
import argparse
import datetime

import xarray as xr
import numpy as np

from zoneinfo import ZoneInfo
from cartopy import crs as ccrs, feature as cfeature
from cartopy.io.img_tiles import GoogleTiles, OSM
from matplotlib.transforms import offset_copy
from matplotlib import pyplot as plt


def mrms_cref(nfile, global_sites):
    """
    Create a display to visualize MRMS Composite Reflectivity across
     the Chicago-CROCUS domain, marking locations of the node. 

    Input
    -----
    nfile - GRIB2 file (uncompressed)
        MRMS Composite Reflectivity grib2 file, which has been uncompressed

    global_sites : dict
        Information that are specific to the CROCUS sites.

    Returns
    -------
    axarr : Matplotlib.pyplot Figure 
        Figure instance containing the desired quicklook; to be saved outside 
        of function
    """
    
    # Read in the MRMS CREF file
    mrms = xr.load_dataset(nfile, engine="cfgrib")

    # Rename the dependent variable from unknown
    mrms = mrms.rename({'unknown' : 'CREF'})

    # Define CREF attributes
    mrms["CREF"].attrs["units"] = "dBZ"
    mrms["CREF"].attrs["long_name"] = "Equivalent Reflectivity"

    # Subset the data
    sdata = mrms.sel(latitude=slice(42.5, 41.5), longitude=slice(271.5, 273)).where(mrms.CREF > -10)

    # Set up the figure
    fig = plt.figure(figsize=(16, 8))
    # Initialize OpenStreetMap tile
    tiler = OSM()
    # Create a subplot and define projection
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

        # plot the CREF, making sure to transform for the projection
    # make sure to set the alpha to not obscure the street map beneath
    sdata.CREF.plot(transform=ccrs.PlateCarree(),
                    vmin=-10, 
                    vmax=65, 
                    ax=ax,
                    alpha=0.78,
                    cmap="ChaseSpectral")
    
    # Add some various map elements to the plot to make it recognizable.
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.STATES)
    ax.add_feature(cfeature.BORDERS)
    ax.add_image(tiler, 11)
    ax.gridlines(draw_labels=True)

    # Set CROCUS - Chicago Domain bounds
    ax.set_extent([271.9, 272.5, 41.6, 42.15])

    # Display the CROCUS node location and tags
    for key in global_sites:
        ax.plot(global_sites[key]['longitude'], 
                global_sites[key]['latitude'], 
                marker='o', 
                color='black', 
                markersize=10, 
                alpha=0.7, 
                transform=ccrs.PlateCarree())

        # Use the cartopy interface to create a matplotlib transform object
        # for the Geodetic coordinate system. We will use this along with
        # matplotlib's offset_copy function to define a coordinate system which
        # translates the text by 25 pixels to the left.
        # note - taken from cartopy examples
        geodetic_transform = ccrs.PlateCarree()._as_mpl_transform(ax)
        text_transform = offset_copy(geodetic_transform, units='dots', x=+50, y=+15)

        # Add text 25 pixels to the left of the volcano.
        ax.text(global_sites[key]['longitude'], 
                global_sites[key]['latitude'], 
                global_sites[key]['site_ID'], 
                verticalalignment='center', 
                horizontalalignment='right', 
                transform=text_transform,
                bbox=dict(facecolor='sandybrown', 
                alpha=0.5, 
                boxstyle='round'))
        
    # update the title of the display
    ax.set_title(np.datetime_as_string(sdata['valid_time'].data, unit='s').replace("T", " - ") + 
                 "Z\n" + "Composite Reflectivity Mosasic")

    # free up memory
    del mrms, sdata, tiler

    return ax


def main(args, global_sites):

    # Initial connection to the Amazon S3 bucket
    fs = fsspec.filesystem("s3", anon=True)

    # Define the potential files for CREF
    file_path = sorted(fs.glob(f"s3://noaa-mrms-pds/CONUS/MergedReflectivityQCComposite_00.50/{args.date}/*"))

    # define a temporary directory to hold static images
    templocation = tempfile.mkdtemp()

    # Iterate through the files, uncompress and call the MRMS CREF function
    for nfile in file_path:
        with fs.open(nfile, 'rb') as gzip_file:
            with tempfile.NamedTemporaryFile(suffix=".grib2") as f:
                # Uncompress the GRIB2 file
                f.write(gzip.decompress(gzip_file.read()))
                # Call the plotting function
                image = mrms_cref(f.name, global_sites)
                # save the figure to a temporary location for later GIF creation
                plt.savefig(templocation + 
                            '/' + 
                            'mrms-radar-cref-' + 
                            nfile.split('_')[-1].split('.')[0] +
                            '.png')
                # free up memory
                del image
                
    # Define files created and define movie path
    map_images = sorted(glob.glob(templocation + f"/mrms-radar-cref*"))
    gif_title = args.outdir + f"mrms-radar-cref-" + args.date + ".gif"

    # Check to see if the file exists - if it does, delete it
    if os.path.exists(gif_title):
        os.remove(gif_title)

    # Loop through and create the gif
    with imageio.get_writer(gif_title, mode='I', duration=0.2) as writer:
        for filename in map_images:
            image = imageio.imread(filename)
            writer.append_data(image)

    # Free Up Memory
    del fs, file_path, templocation, gif_title, map_images
      
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
    
    #put these in a dictionary for accessing
    global_sites = {'NU' : global_NU, 
                    'CSU': global_CSU,
                    'NEIU' : global_NEIU,
                    'UIC' : global_UIC,
                    'ATMOS' : global_ATMOS,
                    'CCICS' : global_CCICS}
    
    descript = ("Generation of Multi-Radar, Multi-Sensor (MRMS) Composite " +
                "Reflectivity Quicklooks and GIFs for the CROCUS Chicago Domain."
    )
    parser = argparse.ArgumentParser(description=descript,
                                     usage="python crocus-mrms-cref.py --start_date 20241031 --outdir /Users/jrobrien/dev/crocus/")
 
    parser.add_argument("--date",
                        type=str,
                        dest='date',
                        default=(datetime.datetime.now(ZoneInfo(waggle_timezone)) - 
                                 datetime.timedelta(days=1)).strftime('%Y%m%d'),
                        help="Date to Start Quicklook Creation in YYYYMMDD [UTC] format"
    )
    parser.add_argument("--outdir",
                        type=str,
                        dest="outdir",
                        default="./",
                        help="Directory to output GIFs. Default is GCE home directory")

    args = parser.parse_args()

    main(args, global_sites)