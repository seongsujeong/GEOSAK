'''
GRWT - GDAL Raster Wrapper and Toolkit

Functionality
 - Simplify I/O procedure of geospatial raster data using GDAL
 - provide simple toolkit to visualize explore, and manipulate the raster

 Derived from grup.py, an internal code for UCI data processing pipeline
  - Major difference: On-demand data loading
    - More pythonic way of defining the class members (i.e. using @property decorator)
    - More modularized code structure
    - Added GAMMA raster support
'''

from numpy.core.numeric import outer
from osgeo import gdal, ogr, osr
import numpy as np
import os
import sys
import argparse
import subprocess


default_raster_options = {'GTiff':["COMPRESS=DEFLATE", "PREDICTOR=1", "ZLEVEL=1"]}
default_format = 'GTiff'

class raster:
    def __init__(self,filename_raster=None):
        #member definition
        self._filename = None
        self._rasterobj = None
        self._str_driver = None
        self._z = None  #image array
        self._nx = None  # Image width
        self._ny = None  # Image height
        self._nz = None  # Number of bands in the raster
        self._GeoTransform = None
        self._Projection = None
        self._nodata = None
        self._option = []
        self.isArea = True
        self.verbose = False  # make the output message a bit more informative. Could be useful for debugging.

        if type(filename_raster) == str:
            self.filename = filename_raster
            #self.load_rasterobj(self.filename)

    @property
    def filename(self):
        return self._filename

    @filename.setter  # NOTE: Setting the filename will cause re-loading the raster object
    def filename(self, filename_raster):
        self._filename = filename_raster
        #self.load_rasterobj()
        # NOTE: Loading the data in this stage might cause disruption when
        # trying to create a new raster file. Suggest not to call this member.


    def load_rasterobj(self):
        self._rasterobj = None  # de-reference the raster object if already loaded
        if self._filename != None:
            if os.path.isfile(self._filename) or \
                ('vsi' in self._filename[:10]) or \
                 ('vsigzip' in self._filename[:10]) or ('vsitar' in self._filename[:10]) or \
                (self._filename.startswith('NETCDF')):  # i.e. gdal-supported raster format
                try:
                    self._rasterobj = gdal.Open(self._filename)

                    #reset the on-demand members
                    self._str_driver = None
                    self._z = None
                    self._nx = None
                    self._ny = None
                    self._nz = None
                    self._GeoTransform = None
                    self._Projection = None
                    self._nodata = None

                except:
                    print('ERROR: Not supported by GDAL - '
                          f'grwt.raster.load_rasterobj(), filename={self._filename}')
                    #try to see if the file is npy array
            elif self._filename.startswith('gamma:'): #
                # NOTE: string format for GAMMA raster:
                # gamma:[filename]:[par file]:dataype[e.g. float32, complex64, uint16, ...]
                self._str_driver = 'gamma'

                filename_data = self._filename.split(':')[1]
                filename_par = self._filename.split(':')[2]

                #TODO: keep implementing the gamma loading routine

            elif self._filename == None:
                print(f'ERROR - Filename was not defined: grwt.raster.load_rasterobj()')

    @property
    def str_driver(self):
        if self._str_driver == None:
            if self._rasterobj == None:
                self.load_rasterobj()
            try:
                self._str_driver = self._rasterobj.GetDriver().GetDescription()
            except:
                self._str_driver = None

        return self._str_driver

    @str_driver.setter
    def str_driver(self, str_in):
        self._str_driver=str_in

    @property
    def nodata(self):
        if self._nodata == None:
            #try to update self._nodata
            if self._filename != None:
                if self._rasterobj == None:
                    self.load_rasterobj()
                arr_nodata = []
                for i in range(self.nz):
                    band_in_raster = self._rasterobj.GetRasterBand(i+1)
                    arr_nodata.append(band_in_raster.GetNoDataValue())
                self._nodata = tuple(arr_nodata)
        return self._nodata

    @nodata.setter
    def nodata(self, number_or_list_or_tuple):
        if type(number_or_list_or_tuple) == list:
            if len(number_or_list_or_tuple) == self.nz:
                self._nodata = tuple(number_or_list_or_tuple)
            else:
                print('ERROR - grwt.raster.nodata.setter - incorrect input list dimension')

        elif type(number_or_list_or_tuple) == tuple:
            if len(number_or_list_or_tuple) == self.nz:
                self._nodata = number_or_list_or_tuple
            else:
                print('ERROR - grwt.raster.nodata.setter - incorrect input tuple dimension')

        else:
            self._nodata = tuple([number_or_list_or_tuple] * self.nz)

    #on-demand processing of the metadata and the data
    @property
    def nx(self): #raster width
        if self._nx == None:
            if self._rasterobj == None:
                self.load_rasterobj()
            self._nx = self._rasterobj.RasterXSize
        return self._nx

    @nx.setter
    def nx(self, nx_in): #To be used when user wants to manually input the raster size;
        self._nx = nx_in

    @property
    def ny(self): #raster width
        if self._ny == None:
            if self._rasterobj == None:
                self.load_rasterobj()
            self._ny = self._rasterobj.RasterYSize
        return self._ny

    @ny.setter
    def ny(self, ny_in): #To be used when user wants to manually input the raster size;
        self._ny = ny_in

    @property
    def z(self):
        if self._z is None:
            if self._rasterobj == None:
                self.load_rasterobj()
            self._z = self._rasterobj.ReadAsArray()
        return self._z

    @z.setter
    def z(self, arr_in, reshape='auto'): #In case the user manually brings the data array
        #NOTE options in reshaping the multiband data: True, False, 'auto'
        if arr_in.ndim == 2: #single-band raster
            self._z = arr_in
            self.nx = self._z.shape[1]
            self.ny = self._z.shape[0]
            self._nz = 1

        elif arr_in.ndim == 3: #multi-band raster:
            #determine whether or not to reshape 'arr_in'
            vec_arr_shape = arr_in.shape
            # smaller 1st element in the shape vector - possible GDAL raster style
            if arr_in.shape[0] < (vec_arr_shape[1] + vec_arr_shape[2]) / 2:
                flag_need_to_reshape = False
            else:
                flag_need_to_reshape = True

            if reshape == True or (flag_need_to_reshape and reshape == 'auto'):
                #perform reshaping
                self._z = np.transpose(arr_in, (1, 2, 0))
            else:
                if flag_need_to_reshape:
                    print('NOTE: grwt.raster.z - Reshaping suggested but rejected by user flag.')
                self._z = arr_in

            self._nx = self._z.shape[2]
            self._ny = self._z.shape[1]
            self._nz = self._z.shape[0]

    @property
    def GeoTransform(self):
        if self._GeoTransform == None:
            if self._rasterobj == None:
                self.load_rasterobj()
            self._GeoTransform = self._rasterobj.GetGeoTransform()
        return self._GeoTransform

    @GeoTransform.setter
    # To be used when the user wants to set the geotransform parameter
    def GeoTransform(self, list_or_tuple):
        if type(list_or_tuple) == list:
            self._GeoTransform = tuple(list_or_tuple)
        else:
            self._GeoTransform = list_or_tuple

    @property
    def Projection(self):
        if self._Projection == None:
            if self._rasterobj == None:
                self.load_rasterobj()
            self._Projection = self._rasterobj.GetProjection()
        return self._Projection

    @Projection.setter
    def Projection(self, epsg_or_filename_or_wkt_or_raster):
        #Input: EPSG number [int]
        #       Another raster instance
        #       GDAL raster filename
        #       proj string
        if type(epsg_or_filename_or_wkt_or_raster) == int: #EPSG
            SR = osr.SpatialReference()
            SR.ImportFromEPSG(epsg_or_filename_or_wkt_or_raster)
            self._Projection = SR.ExportToWkt()
        elif type(epsg_or_filename_or_wkt_or_raster) == raster:  # Another raster instance
            self._Projection = epsg_or_filename_or_wkt_or_raster.Projection
        elif '[' in epsg_or_filename_or_wkt_or_raster: # WKT string
            # TODO: Think about better and more general way of detecting the WKT for projection
            self._Projection = epsg_or_filename_or_wkt_or_raster
        else:
            rasterobj_ref = open(epsg_or_filename_or_wkt_or_raster)
            self._Projection = rasterobj_ref.GetProjection()

    @property
    def nz(self): #number of bands in raster
        if self._nz == None:
            if self._filename == None:
                self._nz = self._z.shape[0]
            else:
                if self._rasterobj == None:
                    self.load_rasterobj()
                self._nz = self._rasterobj.RasterCount #try to give out the value by GDAL header info

        return self._nz

    @property
    def x(self): #grid in x coord.
        return (np.arange(self.nx)).astype(np.float32) * self.GeoTransform[1] + self.GeoTransform[0]

    @property
    def y(self): #grid in y coord.
        return (np.arange(self.ny)).astype(np.float32) * self.GeoTransform[5] + self.GeoTransform[3]

    @property
    def option(self):
        return self._option

    @option.setter
    def option(self, list_option):
        self._option = list_option

    def go_bigtiff(self, force=False):
        # Check if the raster needs to be "bigtiff mode" in case of geotiff
        try:
            size_raster = self.z.dtype.itemsize * self.nx * self.ny
        except:
            size_raster = -1

        if size_raster >= 4 * (2 ** 30) or force: # 4GB
            return True
        elif size_raster < 0: #cannot calculate the array size in bytes
            print('grup.raster.go_bigtiff() - cannot calculate the raster size in byte.')
            return None
        else:
            return False


    def write(self, filename_raster=None, overwrite=False):
        dict_corres_dtype_npy_gdal = {
            'uint8': gdal.GDT_Byte,
            'uint16': gdal.GDT_UInt16,
            'uint32': gdal.GDT_UInt32,
            'int16': gdal.GDT_Int16,
            'int32': gdal.GDT_Int32,
            'float32': gdal.GDT_Float32,
            'float64': gdal.GDT_Float64,
            'complex64': gdal.GDT_CFloat32,
            'complex128': gdal.GDT_CFloat64
        }

        if filename_raster is None:
            filename_out = self._filename
        else:
            filename_out = filename_raster

        if os.path.exists(filename_out) and (not overwrite):
            print('File exists ({}) - Overwriting denied by user'.format(filename_out))

        else:
            if os.path.exists(filename_out):
                os.remove(filename_out)

            if self.str_driver is None:
                self.str_driver = default_format

            drvout = gdal.GetDriverByName(self.str_driver)
            try:
                output = drvout.Create(filename_raster,
                                    self.nx, self.ny, self.nz,
                                    dict_corres_dtype_npy_gdal[str(self.z.dtype)],
                                    self.option)
            except:
                print('grwt.raster.write - Cannot create GDAL output object.')
                output = None

            #define the GeoTransform and Projection info
            if output != None:
                output.SetGeoTransform(self.GeoTransform)
                if self.Projection == None:
                    print('Warning: Projection information was not set.')
                output.SetProjection(self.Projection)

                #write the array and the file
                if self.nz == 1:
                    band_out = output.GetRasterBand(1)
                    if self.nodata is not None:
                        band_out.SetNoDataValue(self.nodata)
                    band_out.WriteArray(self.z)

                else: #multiple band raster
                    for i in range(self.nz):
                        #print('writing: band {} of {}'.format(i+1,self.nz))
                        if self.nodata is not None:
                            output.GetRasterBand(i + 1).SetNoDataValue(self.nodata[i])
                            output.GetRasterBand(i + 1).WriteArray(self.z[i])
                        else:
                            output.GetRasterBand(i + 1).WriteArray(self.z[i])

                #finalize
                output.FlushCache()
                output = None


    def copy(self):
        out = raster()
        out.z = self.z
        out.GeoTransform = self.GeoTransform
        out.Projection = self.Projection
        out.nodata = self.nodata
        out.str_driver = self.str_driver
        out.z = self.z

        return out


    # Miscellaneous members
    @property
    def boundary(self):
        #Return is tuple of (xmin, ymin, xmax, ymax)
        x0 = self.GeoTransform[0]
        y0 = self.GeoTransform[3]
        x1 = x0 + self.GeoTransform[1] * self.nx
        y1 = y0 + self.GeoTransform[5] * self.ny

        if x0 > x1:
            xmin = x1
            xmax = x0
        else:
            xmin = x0
            xmax = x1

        if y0 > y1:
            ymin = y1
            ymax = y0
        else:
            ymin = y0
            ymax = y1

        return (xmin, ymin, xmax, ymax)




######################################## Some utilities ########################################

#Load the magnitude and phase raster files separately, and give out the raster in complex number
#To be used in some DInSAR data processing results
def magphase2complex(magnitude_filename_or_raster, phase_filename_or_raster, isRadian=True):
    if type(magnitude_filename_or_raster) == str:
        raster_mag = raster(magnitude_filename_or_raster)
    else:
        raster_mag = magnitude_filename_or_raster

    if type(phase_filename_or_raster) == str:
        raster_phase = raster(phase_filename_or_raster)
    else:
        raster_phase = phase_filename_or_raster

    if raster_mag.nx == raster_phase.nx and raster_mag.ny == raster_phase.ny:
        raster_out = raster_mag.copy()
        raster_out.z = (np.cos(raster_phase.z) * raster_mag.z) + (np.sin(raster_phase.z) * raster_mag.z) * 1.0j
    else:
        print('ERROR: grwt.magphase2complex() - Input rasters\' dimensions are not the same: '
              f'{raster_mag.z.shape}(mag) vs. {raster_phase.z.shape}(phase)')
    return raster_out


#clip the raster object with the pixel coordinates provided
def clip_raster_by_imgcoord(raster_in, extent_px, raster_out=None):
    # extent_px: 4-element array (e.g. list, tuple, or numpy array) Order is the same as gdalwarp utility (xmin,ymin,xmax,ymax)
    # NOTE - leave raster_out as None when you want the function to return the clipped raster
    #      - Provide the filename into raster_out when you want the raster saved in file system
    #

    if type(raster_in) == str:
        raster_src = raster(raster_in)
    else:
        raster_src = raster_in

    raster_tgt = raster()

    raster_tgt.nx = extent_px[2] - extent_px[0] + 1
    raster_tgt.ny = extent_px[3] - extent_px[1] + 1
    #raster_tgt.nz=raster_src.nz #NOTE: It should be set automatically by the following section in this function

    if raster_src.nz > 1: #multiband
        raster_tgt.z = raster_src.z[:, extent_px[1]:extent_px[3]+1, extent_px[0]:extent_px[2]+1].copy()
    else: #singleband
        raster_tgt.z = raster_src.z[extent_px[1]:extent_px[3]+1, extent_px[0]:extent_px[2]+1].copy()

    #bring the projection information
    raster_tgt.Projection = raster_src.Projection
    dx = raster_src.GeoTransform[1]
    dy = raster_src.GeoTransform[5]
    x0 = raster_src.GeoTransform[0] + dx * extent_px[0]
    y0 = raster_src.GeoTransform[3] + dy * extent_px[1]

    raster_tgt.GeoTransform = (x0, dx, 0, y0, 0, dy)

    if raster_out is None:
        return raster_tgt

    else:
        raster_tgt.write(raster_out)


def clip_to_reference(filename_src, filename_ref, filename_out,
                      resampling='cubic', epsg_out=3031, dryrun=False):
    #TODO: Implement it in more elegant way (i.e. do not rely on shell commands)
    form_command_gdalwarp = ('gdalwarp -r {RS} -tr {RES_X} {RES_Y} -te '
                             '{XMIN} {YMIN} {XMAX} {YMAX} -t_srs epsg:{EPSG} {IN} {OUT}')
    raster_src = raster(filename_src)
    raster_ref = raster(filename_ref)

    xmin = raster_ref.GeoTransform[0]
    xmax = raster_ref.GeoTransform[0] + raster_ref.GeoTransform[1] * raster_ref.nx
    ymax = raster_ref.GeoTransform[3]
    ymin = raster_ref.GeoTransform[3] + raster_ref.GeoTransform[5] * raster_ref.ny

    str_command_gdalwarp = form_command_gdalwarp.format(RS=resampling,
                                                      RES_X=raster_ref.GeoTransform[1],
                                                      RES_Y=abs(raster_ref.GeoTransform[1]),
                                                      XMIN=xmin, YMIN=ymin, XMAX=xmax, YMAX=ymax,
                                                      EPSG=epsg_out) # TODO: soft-code "EPSG=epsg_out"

    if dryrun:
        print(str_command_gdalwarp)
        rtnval = None
    else:
        rtnval = subprocess.call(str_command_gdalwarp, shell=True)

    return None


######################################## GAMMA support ########################################


#parse GAMMA par file into dict
#NOTE:
# - The output dict has a pair whose key is 'header.' This contains the lines at the beginning of the par file
def par2dict(filename_par):

    fields_do_not_split=[
        'title',
        'datum_name',
        'ellipsoid_name',
        'datum_country_list'
    ]

    def parse_string(str_in): #nested function
        #functionality: Convert str_in into variable of appropriate data type i.e. int, float, and string
        if str_in.replace('-','').isnumeric():
            #str_in is integer
            var_out=int(str_in)
        elif str_in.replace('-','').replace('+','').replace('.','').replace('e','').replace('E','').isnumeric(): #regular floating point expression e.g. -24.5254
            #str_in is floating point
            var_out=float(str_in)
        else:
            var_out = str(str_in) #NOTE this casting is not that meaningful...
        return var_out

    dict_out = {}
    dict_out['header'] = []
    with open(filename_par, 'r') as fin:
        lines_in = fin.read().split('\n')

    for line in lines_in:
        if line.replace(' ', '') == '': #blank line
            continue
        elif 'END OF' in line: # "The last line" indicator
            continue
        elif ':' in line:
            #parse the string
            key_and_val = line.split(':')


            if not key_and_val[0] in dict_out.keys(): #create key if does not exist in the dict (i.e. first encounter)
                dict_out[key_and_val[0]] = None
            #Couple of special fields and values
            if key_and_val[0] in fields_do_not_split:
                dict_out[key_and_val[0]] = ':'.join(key_and_val[1:]).lstrip(' ')
            else:
                seg_val = key_and_val[1].split()
                if len(seg_val) == 1:
                    dict_out[key_and_val[0]] = parse_string(seg_val[0])
                else: #series of values
                    if dict_out[key_and_val[0]] == None:
                        dict_out[key_and_val[0]] = []
                        for val_str in seg_val:
                            dict_out[key_and_val[0]].append(parse_string(val_str))

        else:
            #put the string in header
            dict_out['header'].append(line)

    return dict_out


#Returns EPSG number from gcpar
#NOTE - only works for polar stereographic for now...
#TODO - Extend it for UTM and GCS (e.g. EGS84)
def epsg_from_gcpar(dict_par):
    epsg_out = None #Value to return when EPSG detection did not work
    if dict_par['DEM_projection'] == 'PS':
        if dict_par['PS_secant_lat'][0] == -71.0 and dict_par['PS_central_meridian'][0] == -0.0:
            epsg_out = 3031 #Polar stereographic South
        elif dict_par['PS_secant_lat'][0] == 70.0 and dict_par['PS_central_meridian'][0] == -45.0:
            epsg_out = 3413 #Polar stereographic South

    elif dict_par['DEM_projection'] == 'UTM':
        zone_number = dict_par['projection_zone']
        #Choose if the projection should be UTM-north or UTM-south
        #NOTE: In practical perspective, it is not so meaningful - UTM-north and UTM-south use the same reprojection parameters.
        #calculate the image extent
        #x0=dict_par['corner_east']
        ymax = dict_par['corner_north']
        #x1=x0+dict_par['post_east']*dict_par['width']
        ymin = ymax + dict_par['post_north'] * dict_par['nlines']
        if ymin > 0 and ymax > 0:   #northern
            epsg_out = 32600 + zone_number
        elif ymax < 0 and ymin < 0: #southern
            epsg_out = 32700 + zone_number
        else: #equator crosses the image
            if abs(ymin) > ymax:  #consider as southern
                epsg_out = 32700 + zone_number
            else:               #consider as northern
                epsg_out = 32600 + zone_number

    return epsg_out


#Loads GAMMA raster arrays (SLC, MLI, etc.) as grwt.raster instance.
def load_gamma_raster(filename_data, filename_par, dtype=None):
    # Fields in GAMMA par files that contains image shape information
    # format: (one of the lines in the header) - width, nlines, data format
    par_imageshapeinfo = {
        'Interferogram and Image Offset Parameter File': ['offset_estimation_range_samples',
                                                          'offset_estimation_azimuth_samples',
                                                          None],
        'Gamma Interferometric SAR Processor (ISP) - Image Parameter File': ['range_samples',
                                                                             'azimuth_lines',
                                                                             'image_format'],
        'Gamma DIFF&GEO DEM/MAP parameter file': ['width',
                                                  'nlines',
                                                  'data_format'],
        'GAMMA Differential Interferometry (DIFF) DEM parameter file': ['width',
                                                                        'nlines',
                                                                        'data_format'],
        'Gamma DIFF&GEO Processing Parameters': ['map_width',
                                                 'map_azimuth_lines',
                                                 None]
    }
    dict_datatype = {
        'REAL*4': np.float32,
        'REAL*8': np.float64,
        'FCOMPLEX': np.complex64
    }

    #retrieve the image size
    dict_par = par2dict(filename_par)
    fields_format = par_imageshapeinfo[dict_par['header'][-1]]
    nx = dict_par[fields_format[0]]
    ny = dict_par[fields_format[1]]

    if dtype == None:
        datatype_in = dict_datatype[dict_par[fields_format[2]]]
    else:
        datatype_in=dtype

    try:
        raster_out = raster()
        raster_out.str_driver = 'GAMMA'
        raster_out.z = np.fromfile(filename_data, dtype=datatype_in).byteswap().reshape((ny, nx))

        if 'Gamma DIFF&GEO DEM/MAP parameter file' in dict_par['header']: #special treatment for gcpar
            raster_out.Projection = epsg_from_gcpar(dict_par)
            raster_out.GeoTransform = (dict_par['corner_east'][0],
                                       dict_par['post_east'][0],
                                       0,
                                       dict_par['corner_north'][0],
                                       0,
                                       dict_par['post_north'][0])

        return raster_out
    except:
        print('ERROR: grwt.load_gamma_raster() - cannot load the data file')
        return None
