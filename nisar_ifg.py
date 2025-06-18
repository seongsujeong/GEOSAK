# Routines for interferogram computation

import numpy as np
import glob
import os
from osgeo import gdal
import matplotlib.pyplot as plt

GSLC_LAYER_PATH = '/science/LSAR/GSLC/grids/frequencyA/HH'
GUNW_WRAPPED_IFG_PATH = '/science/LSAR/GUNW/grids/frequencyA/wrappedInterferogram/HH/wrappedInterferogram'
GEOTIFF_OPTIONS = ["COMPRESS=DEFLATE", "PREDICTOR=1", "ZLEVEL=1"]

def compute_coherence(ref_gslc_path, sec_gslc_path, coh_out_path, nlooks_x=8, nlooks_y=8):
    ds_ref_gslc = gdal.Open(f'NETCDF:{ref_gslc_path}:{GSLC_LAYER_PATH}', gdal.GA_ReadOnly)
    ds_sec_gslc = gdal.Open(f'NETCDF:{sec_gslc_path}:{GSLC_LAYER_PATH}', gdal.GA_ReadOnly)

    ref_gslc_arr = ds_ref_gslc.ReadAsArray()
    sec_gslc_arr = ds_sec_gslc.ReadAsArray()

    gtf_gslc = ds_ref_gslc.GetGeoTransform()
    srs_gslc = ds_ref_gslc.GetSpatialRef()

    # multi-look
    multilook_shape = (ref_gslc_arr.shape[0] // nlooks_y,
                       ref_gslc_arr.shape[1] // nlooks_x)

    gtf_multilooked = (gtf_gslc[0],
                       gtf_gslc[1] * nlooks_x,
                       gtf_gslc[2],
                       gtf_gslc[3],
                       gtf_gslc[4],
                       gtf_gslc[5] * nlooks_y)

    coherence = compute_coherence_from_arr(ref_gslc_arr, sec_gslc_arr, nlooks_x, nlooks_y)

    # save the GSLC ifg
    drv_out = gdal.GetDriverByName('GTiff')
    ds_out = drv_out.Create(coh_out_path,
                            multilook_shape[1],
                            multilook_shape[0],
                            1,
                            gdal.GDT_Float32,
                            options=GEOTIFF_OPTIONS)

    ds_out.SetGeoTransform(gtf_multilooked)
    ds_out.SetProjection(srs_gslc.ExportToWkt())

    band = ds_out.GetRasterBand(1)
    band.WriteArray(coherence)
    band.FlushCache()

    ds_out = None


def compute_coherence_from_arr(ref_gslc_arr, sec_gslc_arr,nlooks_x=8, nlooks_y=8):
    # multi-look
    multilook_shape = (ref_gslc_arr.shape[0] // nlooks_y,
                       ref_gslc_arr.shape[1] // nlooks_x)

    # pre-allocate summation arrays
    s1s2conj = np.zeros(multilook_shape, dtype=np.complex64)
    s1abs = np.zeros(multilook_shape)
    s2abs = np.zeros(multilook_shape)


    for y0 in range(nlooks_y):
        for x0 in range(nlooks_x):
            s1sub = ref_gslc_arr[y0::nlooks_y, x0::nlooks_x]
            s2sub = sec_gslc_arr[y0::nlooks_y, x0::nlooks_x]

            if s1sub.shape != multilook_shape:
                # trim the sub-array of the full interferogram
                s1sub = s1sub[:multilook_shape[0], :multilook_shape[1]]
                s2sub = s2sub[:multilook_shape[0], :multilook_shape[1]]

            s1abs += np.abs(s1sub) ** 2
            s2abs += np.abs(s2sub) ** 2
            s1s2conj += s1sub * np.conj(s2sub)

    # Compute the coherence
    coherence = np.abs(s1s2conj) / np.sqrt(s1abs * s2abs)

    return coherence


def form_gslc_ifg(ref_gslc_path, sec_gslc_path, ifg_out_path, nlooks_x=8, nlooks_y=8,
                  save_coherence=False):
    '''
    Compute the interferogram from two GSLC products, and multi-look it.

    Parameters
    ----------
    ref_gslc_path : str
        Path to the reference GSLC product.
    sec_gslc_path : str
        Path to the secondary GSLC product.
    ifg_out_path : str
        Path to the output interferogram GeoTIFF file.
    nlooks_x : int, optional
        Number of looks in the range direction. Default is 8.
    nlooks_y : int, optional
        Number of looks in the azimuth direction. Default is 8.

    Returns
    -------
    None
    '''
    ds_ref_gslc = gdal.Open(f'NETCDF:{ref_gslc_path}:{GSLC_LAYER_PATH}', gdal.GA_ReadOnly)
    ds_sec_gslc = gdal.Open(f'NETCDF:{sec_gslc_path}:{GSLC_LAYER_PATH}', gdal.GA_ReadOnly)

    ref_gslc_arr = ds_ref_gslc.ReadAsArray()
    sec_gslc_arr = ds_sec_gslc.ReadAsArray()

    ifg_arr = ref_gslc_arr * np.conj(sec_gslc_arr)

    del ref_gslc_arr, sec_gslc_arr

    gtf_gslc = ds_ref_gslc.GetGeoTransform()
    srs_gslc = ds_ref_gslc.GetSpatialRef()

    # multi-look
    multilook_shape = (ifg_arr.shape[0] // nlooks_y,
                       ifg_arr.shape[1] // nlooks_x)
    ifg_arr_multilooked = np.zeros(multilook_shape, dtype=np.complex64)

    gtf_multilooked = (gtf_gslc[0],
                       gtf_gslc[1] * nlooks_x,
                       gtf_gslc[2],
                       gtf_gslc[3],
                       gtf_gslc[4],
                       gtf_gslc[5] * nlooks_y)

    for y0 in range(nlooks_y):
        for x0 in range(nlooks_x):
            ifg_sub = ifg_arr[y0::nlooks_y, x0::nlooks_x]

            if ifg_sub.shape != multilook_shape:
                # trim the sub-array of the full interferogram
                ifg_sub = ifg_sub[:multilook_shape[0], :multilook_shape[1]]

            ifg_arr_multilooked += ifg_sub

    # save the GSLC ifg
    drv_out = gdal.GetDriverByName('GTiff')
    ds_out = drv_out.Create(ifg_out_path,
                            ifg_arr_multilooked.shape[1],
                            ifg_arr_multilooked.shape[0],
                            1,
                            gdal.GDT_CFloat32,
                            options=GEOTIFF_OPTIONS)

    ds_out.SetGeoTransform(gtf_multilooked)
    ds_out.SetProjection(srs_gslc.ExportToWkt())

    band = ds_out.GetRasterBand(1)
    band.WriteArray(ifg_arr_multilooked)
    band.FlushCache()
    ds_out = None

    if save_coherence:
        coherence_arr = compute_coherence_from_arr(ref_gslc_arr, sec_gslc_arr, nlooks_x, nlooks_y)
        coh_out_path = ifg_out_path.replace('.tif', '_coherence.tif')
        ds_coherence = drv_out.Create(coh_out_path,
                                      ifg_arr_multilooked.shape[1],
                                      ifg_arr_multilooked.shape[0],
                                      1,
                                      gdal.GDT_Float32,
                                      options=GEOTIFF_OPTIONS)
        ds_coherence.SetGeoTransform(gtf_multilooked)
        ds_coherence.SetProjection(srs_gslc.ExportToWkt())

        band_coherence = ds_coherence.GetRasterBand(1)
        band_coherence.WriteArray(coherence_arr)
        band_coherence.FlushCache()
        ds_coherence = None


def form_gslc_ifg_skip(ref_gslc_path, sec_gslc_path, ifg_out_path, stride_x=2, stride_y=4):
    '''
    Form interferogram from GSLC with x / y stride

    Parameters
    ----------
    ref_gslc_path : str
        Path to the reference GSLC product.
    sec_gslc_path : str
        Path to the secondary GSLC product.
    ifg_out_path : str
        Path to the output interferogram GeoTIFF file.
    stride_x : int, optional
        Number of looks in the range direction. Default is 8.
    stride_y : int, optional
        Number of looks in the azimuth direction. Default is 8.

    Returns
    -------
    None

    '''
    ds_ref_gslc = gdal.Open(f'NETCDF:{ref_gslc_path}:{GSLC_LAYER_PATH}', gdal.GA_ReadOnly)
    ds_sec_gslc = gdal.Open(f'NETCDF:{sec_gslc_path}:{GSLC_LAYER_PATH}', gdal.GA_ReadOnly)

    ref_gslc_arr = ds_ref_gslc.ReadAsArray()
    sec_gslc_arr = ds_sec_gslc.ReadAsArray()

    ifg_arr = ref_gslc_arr * np.conj(sec_gslc_arr)

    del ref_gslc_arr, sec_gslc_arr

    gtf_gslc = ds_ref_gslc.GetGeoTransform()
    srs_gslc = ds_ref_gslc.GetSpatialRef()

    gtf_multilooked = (gtf_gslc[0],
                       gtf_gslc[1] * stride_x,
                       gtf_gslc[2],
                       gtf_gslc[3],
                       gtf_gslc[4],
                       gtf_gslc[5] * stride_y)

    ifg_arr_skip = ifg_arr[0::stride_y, 0::stride_x]

    # save the GSLC ifg
    drv_out = gdal.GetDriverByName('GTiff')
    ds_out = drv_out.Create(ifg_out_path,
                            ifg_arr_skip.shape[1],
                            ifg_arr_skip.shape[0],
                            1,
                            gdal.GDT_CFloat32,
                            options=GEOTIFF_OPTIONS)

    ds_out.SetGeoTransform(gtf_multilooked)
    ds_out.SetProjection(srs_gslc.ExportToWkt())

    band = ds_out.GetRasterBand(1)
    band.WriteArray(ifg_arr_skip)
    band.FlushCache()

    ds_out = None


def save_diff_phase_array(ddiff_arr, raster_out_path, reference_raster_path):
    ds_ref_gslc = gdal.Open(reference_raster_path, gdal.GA_ReadOnly)
    gtf_ref = ds_ref_gslc.GetGeoTransform()
    srs_ref = ds_ref_gslc.GetSpatialRef()

    drv_out = gdal.GetDriverByName('GTiff')
    ds_out = drv_out.Create(raster_out_path,
                            ddiff_arr.shape[1],
                            ddiff_arr.shape[0],
                            1,
                            gdal.GDT_Float32,
                            options=GEOTIFF_OPTIONS)

    ds_out.SetGeoTransform(gtf_ref)
    ds_out.SetProjection(srs_ref.ExportToWkt())

    band = ds_out.GetRasterBand(1)
    band.WriteArray(ddiff_arr)
    band.FlushCache()

    ds_out = None


def save_phase(raster_in, phase_raster_out, ifg_arr=None):
    ds_in = gdal.Open(raster_in, gdal.GA_ReadOnly)
    gtf_in = ds_in.GetGeoTransform()
    srs_in = ds_in.GetSpatialRef()

    if ifg_arr is None:
        ifg_arr = ds_in.ReadAsArray()

    drv_out = gdal.GetDriverByName('GTiff')
    ds_out = drv_out.Create(phase_raster_out,
                            ifg_arr.shape[1],
                            ifg_arr.shape[0],
                            1,
                            gdal.GDT_Float32,
                            options=GEOTIFF_OPTIONS)

    ds_out.SetGeoTransform(gtf_in)
    ds_out.SetProjection(srs_in.ExportToWkt())

    band = ds_out.GetRasterBand(1)
    band.WriteArray(np.angle(ifg_arr))
    band.FlushCache()

    ds_out = None


def save_amplitude(raster_in, amp_raster_out, db_scale=False,
                   nlooks_x=2, nlooks_y=4):
    ds_in = gdal.Open(raster_in, gdal.GA_ReadOnly)
    gtf_in = ds_in.GetGeoTransform()
    srs_in = ds_in.GetSpatialRef()

    slc_arr = ds_in.ReadAsArray()

    # multi-look
    multilook_shape = (slc_arr.shape[0] // nlooks_y,
                       slc_arr.shape[1] // nlooks_x)
    amp_arr_multilooked = np.zeros(multilook_shape, dtype=np.complex64)

    gtf_multilooked = (gtf_in[0],
                       gtf_in[1] * nlooks_x,
                       gtf_in[2],
                       gtf_in[3],
                       gtf_in[4],
                       gtf_in[5] * nlooks_y)

    for y0 in range(nlooks_y):
        for x0 in range(nlooks_x):
            slc_subsample = slc_arr[y0::nlooks_y, x0::nlooks_x]
            amp_arr_multilooked += np.abs(slc_subsample)

    drv_out = gdal.GetDriverByName('GTiff')
    ds_out = drv_out.Create(amp_raster_out,
                            amp_arr_multilooked.shape[1],
                            amp_arr_multilooked.shape[0],
                            1,
                            gdal.GDT_Float32,
                            options=GEOTIFF_OPTIONS)

    ds_out.SetGeoTransform(gtf_multilooked)
    ds_out.SetProjection(srs_in.ExportToWkt())

    band = ds_out.GetRasterBand(1)
    if db_scale:
        band.WriteArray(np.log10(amp_arr_multilooked**2) * 10)
    else:
        band.WriteArray(amp_arr_multilooked)
    band.FlushCache()

    ds_out = None

