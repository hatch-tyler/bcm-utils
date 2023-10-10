import pathlib
import fiona
import rasterio
import rasterio.mask

import numpy as np
import pandas as pd
import geopandas as gpd

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
