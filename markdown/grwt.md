
# GRWT - GDAL Raster Wrapper and Toolkit

 (c) Seongsu Jeong @ UC Irvine

 1. Functionality:
    - Provides simplified basic I/O routine for GDAL raster files
    - manipulate the geocoding and projection information easily
    - On-demand loading of data array and other metadata
    - GAMMA raster I/O support *- currently under development*

 2. Denepdency:
    - You will only need to install gdal, which can be easily installed by using conda.

 3. How to install:
    - clone or download GSAK to the system.
    - add the GSAK location to $PYTHONPATH
    - import grwt in python code, console, or jupyter notebook
    - That's it!

 4. Note
    - This module is part of GSAK - GeoSpatial 
    - This code is currently under development as the code maintainer finds a new necessity of adding a new feature.
    - Default file format is Geotiff.

 5. Example

    5.1. Read the raster file
    - Provide the file name to the instance.
    - The data will be loaded on demand i.e. as soon as the user refers to the member in the instance.

    ```python
    import grwt
    i0=grwt.raster()
    i0.filename='input.tif'
    i1=grwt.raster('input.tif') #equivalent to the two lines codes right above
    
    print(i0.nx) #Raster width
    print(i0.ny) #Raster height
    print(i0.nz) #Number of bands in the raster

    print(i0.x) #x grid in map coordinate
    print(i1.y) #x grid in map coordinate

    print(i1.z.shape) # Should be equal to (i1.nz,i1.ny,i1.nx)
    print(i0.GeoTransform) #GeoTransform tuple
    print(i0.Projection) #Projection in WKT
    ```

    5.2. Save the raster file
    - use `raster.write()` to save file. Calling this member without any argument will attempt to write those filename is defined in `raster.filename`. Note that it will not work because overwriting is disabled by default. To enable overwriting, try `raster.write(overwrite=True)`.

    ```python
    i0.write() #It will usually not work unless overwrite=True is provided.
    i0.write(overwrite=True) #This should work
    i0.write('output.tif')
    ```

    5.3. Copy the raster
    - Doing this will copy every data from one to another, including array data, GeoTransformation, Projection etc.

    ```python
    i0=grwt.raster('input.tif')
    i1=i0.copy()
    i1.write('input_copy.tif')
    ```

    5.4. Manipulate the projection:
    - Three ways of doing it: See the examples below

    ```python
    i0.Projection(3031) #Set up by EPSG number
    i0.Projection('reference.tif') # Brings from another GDAL raster file
    i0.Projection(i1) # Brings from another grwt.raster instance
    ```
