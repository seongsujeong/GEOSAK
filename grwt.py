#GRWT - GDAL Raster Wrapper and Toolkit


#Functionality
# - Simplify I/O procedure of geospatial raster data using GDAL
# - provide simple toolkit to visualize explore, and manipulate the raster
#
# Derived from grup.py, an internal code for UCI data processing pipeline
#  - Major difference: On-demand data loading
#


from numpy.core.numeric import outer
from osgeo import gdal,ogr,osr
import numpy as np
import os
import sys
import argparse



#Default raster options dfor each fiel format
default_raster_options={'GTiff':['COMPRESS=LZW']}
default_format='GTiff'

class raster:
    def __init__(self,filename_raster=None):
        #member definition
        self._filename=None
        self._rasterobj=None
        self._str_driver=None
        self._z=None #image array
        self._nx=None # Image width
        self._ny=None # Image height
        self._nz=None # Number of bands in the raster
        self._GeoTransform=None
        self._Projection=None
        self._nodata=None
        self._option=None
        self.isArea=True
        self.verbose=False #make the output message a bit more informative. Could be useful for debugging.

        if type(filename_raster)==str:
            self.filename=filename_raster
            #self.load_rasterobj(self.filename)

    @property
    def filename(self):
        return self._filename

    @filename.setter #NOTE: Setting the filename will cause re-loading the raster object
    def filename(self,filename_raster):
        self._filename=filename_raster
        self.load_rasterobj()


    def load_rasterobj(self):
        self._rasterobj=None #de-reference the raster object if alreade loaded
        if self._filename!=None:
            if os.path.isfile(self._filename) or \
                ('vsi' in self._filename[:10]) or ('vsigzip' in self._filename[:10]) or ('vsitar' in self._filename[:10]) or \
                (self._filename.startswith('NETCDF')): #i.e. gdal-supported raster format
                try:
                    self._rasterobj=gdal.Open(self._filename)
                    
                    #reset the on-demand menbers
                    self._str_driver=None
                    self._z=None
                    self._nx=None
                    self._ny=None
                    self._nz=None
                    self._GeoTransform=None
                    self._projection=None
                    self._nodata=None
                    
                    #self._str_driver=self._rasterobj.GetDriver().GetDescription()
                except:
                    print('ERROR: Not supported by GDAL - grwt.raster.load_rasterobj(), filename={}'.format(self._filename))
                    #try to see if the file is npy array
            elif self._filename.startswith('gamma:'): #
                #NOTE: string format for GAMMA raster: gamma:[filename]:[par file]:dataype[e.g. float32, complex64, uint16, ...]
                self._str_driver='gamma'
                
                filename_data=self._filename.split(':')[1]
                filename_par=self._filename.split(':')[2]

                #TODO: keep implementing the gamma loading routine
                
            elif self._filename==None:
                print(print('ERROR - Filename was not defined: grwt.raster.load_rasterobj()'.format(self._filename)))

    @property
    def str_driver(self):
        if self._str_driver==None:
            if self._rasterobj==None:
                self.load_rasterobj()
            try:
                self._str_driver=self._rasterobj.GetDriver().GetDescription()
            except:
                self._str_driver=None
            
        return self._str_driver
        
    @str_driver.setter
    def str_driver(self,str_in):
        self._str_driver=str_in

    @property
    def nodata(self):
        if self._nodata==None:
            #try to update self._nodata
            if self._filename!=None:
                if self._rasterobj==None:
                    self.load_rasterobj()
                arr_nodata=[]
                for i in range(self.nz):
                    band_in_raster=self._rasterobj.GetRasterBand(i+1)
                    arr_nodata.append(band_in_raster.GetNoDataValue())
                self._nodata=tuple(arr_nodata)
        return self._nodata
                
    @nodata.setter
    def nodata(self,number_or_list_or_tuple):
        if type(number_or_list_or_tuple)==list:
            if len(number_or_list_or_tuple)==self.nz:
                self._nodata=tuple(number_or_list_or_tuple)
            else:
                print('ERROR - grwt.raster.nodta.setter - incorrect input list dimension')
        
        elif type(number_or_list_or_tuple)==tuple:
            if len(number_or_list_or_tuple)==self.nz:
                self._nodata=number_or_list_or_tuple
            else:
                print('ERROR - grwt.raster.nodta.setter - incorrect input tuple dimension')

        else:
            self._nodata=tuple([number_or_list_or_tuple]*self.nz)
        




    #on-demand processing of the metadata and the data
    @property
    def nx(self): #raster width
        if self._nx==None:
            if self._rasterobj==None:
                self.load_rasterobj()
            self._nx=self._rasterobj.RasterXSize
        return self._nx
    
    @nx.setter
    def nx(self,nx_in): #To be used when user wants to manually input the raster size; 
        self._nx=nx_in


    @property
    def ny(self): #raster width
        if self._ny==None:
            if self._rasterobj==None:
                self.load_rasterobj()
            self._ny=self._rasterobj.RasterYSize
        return self._ny
    
    @ny.setter
    def ny(self,ny_in): #To be used when user wants to manually input the raster size; 
        self._ny=ny_in


    @property
    def z(self):
        if self._z is None:
            if self._rasterobj==None:
                self.load_rasterobj()
            self._z=self._rasterobj.ReadAsArray()
        return self._z

    @z.setter
    def z(self,arr_in,reshape='auto'): #In case the user manually brings the data array
        #NOTE options in reshaping the multiband data: True, False, 'auto'
        if len(arr_in.shape)==2: #single-band raster
            self._z=arr_in
            self.nx=self._z.shape[1]
            self.ny=self._z.shape[0]
            self.nz=1

        elif len(arr_in.shape)==3: #multi-band raster:
            #determine whether or not to reshape 'arr_in'
            
            vec_arr_shape=arr_in.shape
            if arr_in.shape[0]<(vec_arr_shape[1]+vec_arr_shape[2])/2: #smaller 1st element in the shape vector - possible GDAL raster style
                flag_need_to_reshape=False
            else:
                flag_need_to_reshape=True
            
            if reshape==True or (flag_need_to_reshape and reshape=='auto'):
                #perform reshaping
                self._z=np.transpose(arr_in,(1,2,0))
            else:
                if(flag_need_to_reshape):
                    print('NOTE: grwt.raster.z - Reshaping suggested but rejected by user flag.')
                self._z=arr_in

            self._nx=self._z.shape[2]
            self._ny=self._z.shape[1]
            self._nz=self._z.shape[0]


            
    @property
    def GeoTransform(self):
        if self._GeoTransform==None:
            if self._rasterobj==None:
                self.load_rasterobj()
            self._GeoTransform=self._rasterobj.GetGeoTransform()
        return self._GeoTransform

    @GeoTransform.setter
    def GeoTransform(self,list_or_tuple): #To be used when the user wants to set the geotransform parameter
        if type(list_or_tuple)==list:
            self._GeoTransform=tuple(list_or_tuple)
        else:
            self._GeoTransform=list_or_tuple

    @property
    def Projection(self):
        if self._Projection==None:
            if self._rasterobj==None:
                self.load_rasterobj()
            self._Projection=self._rasterobj.GetProjection()
        return self._Projection

    @Projection.setter
    def Projection(self,epsg_or_filename_or_wkt_or_raster):
        #Input: EPSG number [int]
        #       Another raster instance
        #       GDAL raster filename
        #       proj string
        if type(epsg_or_filename_or_wkt_or_raster)==int: #EPSG
            SR=osr.SpatialReference()
            SR.ImportFromEPSG(epsg_or_filename_or_wkt_or_raster)
            self._Projection=SR.ExportToWkt()
        elif type(epsg_or_filename_or_wkt_or_raster)==raster: #Another raster instance
            self._Projection=epsg_or_filename_or_wkt_or_raster.Projection
        elif 'PROJCS' in epsg_or_filename_or_wkt_or_raster: #WKT string
            self._Projection=epsg_or_filename_or_wkt_or_raster
        else:
            rasterobj_ref=open(epsg_or_filename_or_wkt_or_raster)
            self._Projection=rasterobj_ref.GetProjection()


    @property
    def nz(self): #number of bands in raster
        if self._nz==None:
            if self._filename==None:
                self._nz=self._z.shape[0]
            else:
                if self._rasterobj==None:
                    self.load_rasterobj()
                self._nz=self._rasterobj.RasterCount #try to give out the value by GDAL header info
                
        return self._nz

    @property
    def x(self): #grid in x coord. 
        return (np.arange(self.nx)).astype(np.float32)*self.GeoTransform[1]+self.GeoTransform[0]

    @property
    def y(self): #grid in y coord. 
        return (np.arange(self.ny)).astype(np.float32)*self.GeoTransform[5]+self.GeoTransform[3]

    @property
    def option(self):
        if self._option==None:
            return []
        else:
            return self._option

    def go_bigtiff(self,force=False): #Check if the raster needs to be "bigtiff mode" in case of geotiff
        try:
            size_raster=self.z.dtype.itemsize*self.nx*self.ny
        except:
            size_raster=-1
            
        if size_raster>= 4*(2**30) or force: #4GB
            return True
        elif size_raster<0: #cannot calculate the array size in bytes
            print('grup.raster.go_bigtiff() - cannot calculate the raster size in byte.')
            return None
        else:
            return False


    def write(self,filename_raster=None,overwrite=False):
        dict_corres_dtype_npy_gdal={
            'uint8':gdal.GDT_Byte,
            'uint16':gdal.GDT_UInt16,
            'uint32':gdal.GDT_UInt32,
            'int16':gdal.GDT_Int16,
            'int32':gdal.GDT_Int32,
            'float32':gdal.GDT_Float32,
            'float64':gdal.GDT_Float64,
            'complex64':gdal.GDT_CFloat32,
            'complex128':gdal.GDT_CFloat64
        }

        if filename_raster==None:
            filename_out=self._filename
        else:
            filename_out=filename_raster

        if os.path.exists(filename_out) and (not overwrite):
            print('File exists ({}) - Overwriting denied by user'.format(filename_out))
            
        else:
            if os.path.exists(filename_out):
                os.remove(filename_out)
            
            if self.str_driver==None:
                self.str_driver=default_format

            drvout=gdal.GetDriverByName(self.str_driver)
            try:
                output=drvout.Create(filename_raster,\
                                    self.nx,self.ny,self.nz,\
                                    dict_corres_dtype_npy_gdal[str(self.z.dtype)],
                                    self.option)
            except:
                print('grwt.raster.write - Cannot create GDAL output object.')
                output=None
            
            #define the GeoTransform and Projection info
            if output!=None:
                output.SetGeoTransform(self.GeoTransform)
                if self.Projection==None:
                    print('Warning: Projection information was not set.')
                output.SetProjection(self.Projection)
                
                #write the array and the file
                #if self.nz==1:
                #    band_out=output.GetRasterBand(1)
                #    if self.nodata!=None:
                #        band_out.SetNoDataValue(self.nodata)
                #    band_out.WriteArray(self.z)
                # 
                #else: #multiple band raster
                for i in range(self.nz):
                    print('writing: band {} of {}'.format(i+1,self.nz))
                    if self.nodata is not None:
                        output.GetRasterBand(i+1).SetNoDataValue(self.nodata[i])
                        output.GetRasterBand(i+1).WriteArray(self.z[i])
                    else:
                        output.GetRasterBand(i+1).WriteArray(self.z[i])

                #finalize
                output.FlushCache()
                output=None







    def copy(self):
        out=raster()
        out.z=self.z
        out.GeoTransform=self.GeoTransform
        out.Projection=self.Projection
        out.nodata=self.nodata
        out.str_driver=self.str_driver
        out.z=self.z

        return out
        

    # Miscellaneous members
    @property
    def boundary(self):
        #Return is tuple of (xmin, ymin, xmax, ymax)
        x0=self.GeoTransform[0]
        y0=self.GeoTransform[3]
        x1=x0+self.GeoTransform[1]*self.nx
        y1=y0+self.GeoTransform[5]*self.ny

        if x0>x1:
            xmin=x1
            xmax=x0
        else:
            xmin=x0
            xmax=x1

        if y0>y1:
            ymin=y1
            ymax=y0
        else:
            ymin=y0
            ymax=y1
        
        return (xmin,ymin,xmax,ymax)
        






if __name__=='__main__':
    ################ TEST CODE ###############
    i0=raster('/Users/seongsu/Desktop/GRWT_TEST/coco20160313_20160314-20160601_20160602.tif')
    print(i0.nodata)
    i0.write('/Users/seongsu/Desktop/GRWT_TEST/coco20160313_20160314-20160601_20160602_copy.tif')
    i1=i0.copy() 
    i1.write('/Users/seongsu/Desktop/GRWT_TEST/coco20160313_20160314-20160601_20160602_copy2.tif')










