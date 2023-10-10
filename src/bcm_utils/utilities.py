import fiona
import rasterio
import rasterio.mask

import numpy as np
import pandas as pd
import geopandas as gpd

from pathlib import Path
from shapely.geometry import Polygon


def create_bounding_polygon(area_of_interest, crs=None):
    """
    Create a bounding polygon for an area of interest

    Parameters
    ----------
    area_of_interest : str
        file path and name of a geopandas compatible spatial format

    crs : str or None
        coordinate reference system compatible with geopandas

    Returns
    -------
    gpd.GeoDataFrame
        GeoDataFrame containing the bounding polygon of the area of interest
    """
    # read area of interest into a geopandas GeoDataFrame
    aoi = gpd.read_file(area_of_interest)

    # check if a crs exists for the area of interest
    if not aoi.crs and crs:
        aoi.crs = crs

    # get the minimum and maximum x- and y- coordinates
    minx, miny, maxx, maxy = aoi.total_bounds

    # create a shapely polygon of the bounding box created by the combinations of minimum and maximum x- and y- coordinates
    bounding_polygon = Polygon([(minx, miny), (minx, maxy), (maxx, maxy), (maxx, miny)])

    # convert the shapely geometry to a geopandas GeoDataFrame
    bp = gpd.GeoDataFrame(geometry=[bounding_polygon], crs=aoi.crs)

    return bp


def clip_raster(
    raster: str,
    raster_crs: str,
    clip_boundary: str,
    out_path: str,
    out_ext: str | None = None,
) -> None:
    """
    Clip a raster based on a provided clip boundary

    Parameters
    ----------
    raster : str
        file path and name of raster to clip

    raster_crs: str
        coordinate reference system for input raster

    clip_boundary: str
        geopandas-compatible format boundary to clip raster

    out_path : str
        file path to write clipped raster

    out_ext : str | None, default None
        raster file extension.

        .. note valid formats can be determined via rasterio.drivers.raster_driver_extensions()

    Returns
    -------
    None
        writes the clipped raster to the output location
    """

    # get a list of geometries to use for clipping
    shapes = gpd.read_file(clip_boundary).geometry.tolist()

    with rasterio.open(raster, crs=raster_crs) as src:
        out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True)
        out_meta = src.meta

    if not out_ext:
        # if an output file extension is not provided use the input raster format to specify the driver
        driver = rasterio.drivers.driver_from_extension(raster)
        raster_name = Path(raster).name

    else:
        # get the raster format driver from rasterio based on the provided raster file extension
        raster_drivers = rasterio.drivers.raster_driver_extensions()
        driver = raster_drivers[out_ext]
        raster_name = f"{Path(raster).stem}.{out_ext}"

    out_meta.update(
        {
            "driver": driver,
            "height": out_image.shape[1],
            "width": out_image.shape[2],
            "transform": out_transform,
            "crs": raster_crs,
        }
    )

    # output raster path and name
    clipped_raster = Path(out_path).joinpath(raster_name)

    # write clipped raster to destination if it does not exist
    if not clipped_raster.exists():
        with rasterio.open(clipped_raster, "w", **out_meta) as dest:
            dest.write(out_image)


if __name__ == "__main__":
    raster = r"D:\INTERA\BCM\rch\rch1899dec.asc"
    bcm_crs = "EPSG:3310"
    bp = r"D:\INTERA\Projects\Borrego\GIS\bounding_polygon.shp"
    out_ext = "tif"
    out_path = r"D:\INTERA\Projects\Borrego\BCM\test"

    clip_raster(raster, bcm_crs, bp, out_path)
