# -*- coding: utf-8 -*-

"""
***************************************************************************
*                                                                         *
*   QGIS Snowslide Tool v.2.0.0-beta                                      *
#   QQGIS-Version: 3.28.15-Firenze                                        * 
*   pyAlpineRisk                                                          *
*   Nicole Kamp & Franz Langegger                                         *
*   niki.kamp@gmail.com                                                   *
*   www.nicolekamp.com                                                    *
*   December 2022                                                         *
*   revised 03/2024                                                       *
*                                                                         *
*   Bugs/ ToDos:                                                          *
*   1. rebuild processes with gdal and geopandas                          *
*   2. EPSG codes and transformations                                     *
*   3. remodel step 5                                                     *
*   4. Error, when using GDAL instead of PIL (pix_dtm)                     *
***************************************************************************
"""

from qgis.PyQt.QtCore import QCoreApplication
from qgis.PyQt.QtCore import QFileInfo
from qgis.core import (QgsProcessing,
                       QgsProject,
                       QgsVectorLayer,
                       QgsTextFormat,
                       QgsExpression,
                       QgsFeatureRequest,
                       QgsFeature,
                       QgsGeometry,
                       QgsPoint,
                       QgsPointXY,
                       QgsVectorFileWriter,
                       QgsRasterBandStats,
                       QgsColorRampShader,
                       QgsRasterTransparency,
                       QgsFillSymbol,
                       QgsLineSymbol,
                       QgsSymbol,
                       QgsRasterShader,
                       QgsSingleBandPseudoColorRenderer,
                       QgsRendererCategory,
                       QgsMarkerSymbol,
                       QgsCategorizedSymbolRenderer,
                       QgsSimpleFillSymbolLayer,
                       QgsWkbTypes,
                       QgsVectorLayerSimpleLabeling,
                       QgsRuleBasedLabeling,
                       QgsPalLayerSettings,
                       QgsFeatureSink,
                       QgsProcessingException,
                       QgsProcessingAlgorithm,
                       QgsProcessingParameterFeatureSource,
                       QgsProcessingParameterRasterLayer,
                       QgsProcessingParameterNumber,
                       QgsProcessingParameterString,
                       QgsProcessingParameterFile,
                       QgsRasterLayer,
                       QgsCoordinateReferenceSystem,
                       QgsProcessingParameterFolderDestination,
                       QgsProcessingParameterFeatureSink,
                       QgsProcessingParameterExtent,
                       QgsProcessingMultiStepFeedback,
                       QgsProcessingParameterEnum,
                       QgsMessageLog)
from qgis import processing
from qgis.analysis import QgsRasterCalculatorEntry, QgsRasterCalculator
from qgis.PyQt.QtWidgets import QApplication
from qgis.PyQt.QtGui import QIcon
from qgis.PyQt.QtGui import QColor
from qgis.utils import iface

from PIL import Image
Image.MAX_IMAGE_PIXELS = None
import osgeo.gdal as gdal
import numpy as np

from osgeo import ogr, osr
from shapely.geometry import Polygon
from shapely.geometry import Point
from shapely.ops import unary_union
import geopandas as gpd
import matplotlib
import subprocess
from subprocess import call

import string, os, sys, copy, shutil, math, numpy, time, datetime, glob
from time import *
from sys import *


##---------------------------------------------------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------------------------------------------------##
###FUNCTIONS
# Get Elevation Value
def get_elevation(x, y, band, gt):
    px = int((x - gt[0]) / gt[1])
    py = int((y - gt[3]) / gt[5])

    # Are coordinates within raster
    if px < 0 or py < 0 or px >= band.XSize or py >= band.YSize:
        return float('nan')  # Outside the raster
    
    elevation = band.ReadAsArray(px, py, 1, 1)
    if elevation is None:
        return float('nan')

    return elevation[0][0]
    
    
    
##---------------------------------------------------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------------------------------------------------##
class GleitschneeProcessingAlgorithm(QgsProcessingAlgorithm):
    """
GIS-gestützte Reichweiten Ermittlung von Schneerutschprozessen
Eingabe der Input- und Output-Parameter
    """
    INPUT_Shape = 'INPUT_SHP'
    INPUT_ALS = 'DGM'
    #Richtungsgrad = 'RG'
    S1 = 'S1'
    S2 = 'S2'
    S3 = 'S3'
    SELECTION = 'Selection'
    RESAMPLE = 'Aufloesung'
    TEMP = 'TEMP'
    
    SELECTIONS = ['Linie', 'Polygon']
    SELECTIONS2 = ['1m','2m','3m','5m']
    
    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return GleitschneeProcessingAlgorithm()

    def name(self):
        return 'snowslidetool'

    def displayName(self):
        return self.tr('Snowslide_Tool_v_2_0_0-beta')

    def group(self):
        return self.tr('pyAlpineRisk')

    def groupId(self):
        return 'pyAlpineRiskScripts'

    def shortHelpString(self):
        return self.tr("GIS-gestützte Reichweiten Ermittlung von Schneerutschprozessen (Langegger & Kamp, Version 1.2).\n"
        "\n\n"
        "ARBEITSSCHRITTE\n"
        "1. Abfragelinie oder Fläche erstellen (QGIS Werkzeuge „Neue Shapedatei“ oder QGIS Erweiterung „QDRAW“)\n"
        "2. SchneerutschTool starten (Bedienfeld – Verarbeitungswerkzeuge – Skripte)\n"
        "3. Input-Geometrie und Datensatz auswählen\n"
        "4. erforderliche Gefällsstufen festlegen\n"
        "5. Datengrundlage (ALS DGM) und Rasterauflösung festlegen\n"
        "6. Optional: Output-Ordner festlegen (Vorgabe: „C:/temp/GT“)\n"
        "7. Berechnung starten\n"
        "8. Das Ergebnis wird in der entsprechenden Darstellung automatisch zur aktuellen QGIS Oberfläche hinzugefügt\n"
        )
        

    def initAlgorithm(self, config=None):
        self.selections = [self.tr('Linie'),
                        self.tr('Polygon')]
        self.selections2 = [self.tr('1m'),
                        self.tr('2m'),
                        self.tr('3m'),
                        self.tr('5m')]

        self.addParameter(
            QgsProcessingParameterEnum(
                self.SELECTION,
                self.tr('Input-Geometrie'),
                self.selections,
                defaultValue=0
            )
        )
        
        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.INPUT_Shape,
                self.tr('Anriss'),
                [QgsProcessing.TypeVectorAnyGeometry]
            )
        )
        
        #self.addParameter(
            #QgsProcessingParameterNumber(
                #self.Richtungsgrad,
                #self.tr('Höhendifferenz zu Startpunkt (Beginn Richtungsverlängerung)'),
                #QgsProcessingParameterNumber.Double,
                #10
            #)
        #)
        
        self.addParameter(
            QgsProcessingParameterNumber(
                self.S1,
                self.tr('Pauschalgefälle in Grad - Stufe 1'),
                QgsProcessingParameterNumber.Double,
                30
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.S2,
                self.tr('Pauschalgefälle in Grad - Stufe 2'),
                QgsProcessingParameterNumber.Double,
                28
            )
        )
        
        self.addParameter(
            QgsProcessingParameterNumber(
                self.S3,
                self.tr('Pauschalgefälle in Grad - Stufe 3'),
                QgsProcessingParameterNumber.Double,
                26
            )
        )

        self.addParameter(
            QgsProcessingParameterFile(
                self.INPUT_ALS,
                self.tr('Geländemodell (DGM)'),
                defaultValue='DGM.tif',
            )
        )

        self.addParameter(
            QgsProcessingParameterEnum(
                self.RESAMPLE,
                self.tr('Soll die Auflösung des DGM verändert werden?'),
                self.selections2,
                defaultValue=1
            )
        )

        self.addParameter(
            QgsProcessingParameterFolderDestination(
                self.TEMP, 
                self.tr('Output-Ordner'), 
                defaultValue='C:/temp/GT'
            )
        )
#-------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------#
    def processAlgorithm(self, parameters, context, feedback):
        feedback = QgsProcessingMultiStepFeedback(1, feedback)
        results = {}
        outputs = {}
        
        laenge = 200

        # Define paths + timestamp + logfile
        timestamp = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
        temp_path = str(parameters[self.TEMP])+'/temp_'+str(timestamp)
        final_path = str(parameters[self.TEMP])+'/final_'+str(timestamp)
        
        if not os.path.exists(temp_path):
            os.makedirs(temp_path)
        if not os.path.exists(final_path):
            os.makedirs(final_path)
        
        logfile = open(final_path + '/' + 'logfile.txt','w')
        logfile.write("Gleitschnee-Tool gestartet" + "\n")
        logfile.write(str(timestamp)+"\n")

        # Get EPSG Code
        project = QgsProject.instance()
        crs_project = project.crs()
        epsg_project = crs_project.authid()
        epsg_project = epsg_project.split(':')[1]
        feedback.pushInfo(str(epsg_project))
        epsg = 'EPSG:'+str(epsg_project)        
        
        #-------------------------------------------------------------------------------------------------------------------------------------------#
        ## 1. Save Input-Shapefile
        input_shapefile = temp_path + '/' + 'input_shapefile.shp'
        alg_params = {
            'INPUT':str(parameters[self.INPUT_Shape]),
            'OUTPUT':input_shapefile,
            'LAYER_NAME':input_shapefile,
            'DATASOURCE_OPTIONS':'',
            'LAYER_OPTIONS':''
        }
        processing.run("native:savefeatures", alg_params, context=context, feedback=feedback, is_child_algorithm=True)

        #-------------------------------------------------------------------------------------------------------------------------------------------#
        ## 2. Height cutting in the area of the scribed surface
        # DGM        
        # Line or polygon
        sel = self.SELECTIONS[self.parameterAsEnum(parameters, self.SELECTION, context)]
        ## Line 
        if str(sel) == str('Linie'):
            Line_calc = temp_path + '/' + 'line_calc.shp'
            out_hl = temp_path + '/' + 'hsline.shp'
            line_calc_path = input_shapefile
            feedback.pushInfo(str(line_calc_path))
            gdf = gpd.read_file(line_calc_path)
            gdf['NR'] = 1
            gdf.to_file(Line_calc)
        
            # Resolve lines
            gdf['geometry'] = gdf.geometry.apply(lambda x: unary_union(x))
            gdf_dissolved = gdf.dissolve(by='NR')
            gdf_dissolved.to_file(out_hl)
            
            del gdf

        ## Polygon        
        if str(sel) == str('Polygon'):
            out_dgm_hl = temp_path + '/' + 'dgm_hl.tif'
            alg_params = {
                'ALPHA_BAND': False,
                'CROP_TO_CUTLINE': True,
                'DATA_TYPE': 0,
                'EXTRA': '',
                'INPUT': str(parameters[self.INPUT_ALS]),
                'KEEP_RESOLUTION': False,
                'MASK': input_shapefile,
                'MULTITHREADING': False,
                'NODATA': None,
                'OPTIONS': '',
                'SET_RESOLUTION': False,
                'SOURCE_CRS': QgsCoordinateReferenceSystem(epsg),
                'TARGET_CRS': QgsCoordinateReferenceSystem(epsg),
                'X_RESOLUTION': None,
                'Y_RESOLUTION': None,
                'OUTPUT': out_dgm_hl
            }
            processing.run('gdal:cliprasterbymasklayer', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

            # Contour
            hl = temp_path + '/' + 'hl.shp'
            alg_params = {
                'BAND': 1,
                'CREATE_3D': False,
                'EXTRA': '',
                'FIELD_NAME': 'ELEV',
                'IGNORE_NODATA': False,
                'INPUT': out_dgm_hl,
                'INTERVAL': 3,
                'NODATA': None,
                'OFFSET': 0,
                'OUTPUT': hl
            }
            outputs['Contour'] = processing.run('gdal:contour', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

            # Linien erstellen
            Line_calc = temp_path + '/' + 'line_calc.shp'
            out_hl = temp_path + '/' + 'hsline.shp'
            gdf = gpd.read_file(hl)
            gdf['NR'] = 1
            gdf.to_file(Line_calc)

            # Resolve lines
            gdf['geometry'] = gdf.geometry.apply(lambda x: unary_union(x))
            gdf_dissolved = gdf.dissolve(by='NR')
            gdf_dissolved.to_file(out_hl)
            
            del gdf

        #-------------------------------------------------------------------------------------------------------------------------------------------#
        ## 3. Define study area
        buffer = out = temp_path + '/' + 'buffer.shp'
        alg_params = {
            'DISSOLVE': False,
            'DISTANCE': 350,
            'END_CAP_STYLE': 2,
            'INPUT': out_hl,
            'JOIN_STYLE': 1,
            'MITER_LIMIT': 2,
            'SEGMENTS': 1,
            'OUTPUT': buffer
        }
        processing.run('native:buffer', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        logfile.write("Study area defined"+"\n")

        #-------------------------------------------------------------------------------------------------------------------------------------------#
        ## 4. Crop grid to layer mask
        out = temp_path + '/' + 'dgm_orig.tif'
        res = temp_path + '/' + 'dgm_res.tif'

        fileInfo = QFileInfo(str(parameters[self.INPUT_ALS]))
        absolutePath = fileInfo.absolutePath()
        baseName = fileInfo.baseName()
        asc = str(absolutePath)+'/'+str(baseName)+'.asc'
        asc = asc.replace("/","(\)")
        asc = asc.replace("(","")
        asc = asc.replace(")","")
        if str(parameters[self.INPUT_ALS]) == str(asc):
            alg_params = {
                'ALPHA_BAND': False,
                'CROP_TO_CUTLINE': True,
                'DATA_TYPE': 0,
                'EXTRA': '',
                'INPUT': str(parameters[self.INPUT_ALS]),
                'KEEP_RESOLUTION': False,
                'MASK': str(buffer),
                'MULTITHREADING': False,
                'NODATA': None,
                'OPTIONS': '',
                'SET_RESOLUTION': False,
                'SOURCE_CRS': None,
                'TARGET_CRS': QgsCoordinateReferenceSystem(epsg),
                'X_RESOLUTION': None,
                'Y_RESOLUTION': None,
                'OUTPUT': res
            }
            processing.run('gdal:cliprasterbymasklayer', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

            # Reproject Raster
            alg_params = {
                'DATA_TYPE': 0,
                'EXTRA': '',
                'INPUT': res,
                'MULTITHREADING': False,
                'NODATA': None,
                'OPTIONS': '',
                'RESAMPLING': 0,
                'SOURCE_CRS': QgsCoordinateReferenceSystem(epsg),
                'TARGET_CRS': QgsCoordinateReferenceSystem(epsg),
                'TARGET_EXTENT': None,
                'TARGET_EXTENT_CRS': None,
                'TARGET_RESOLUTION': None,
                'OUTPUT': out
            }
            outputs['TransformierenReprojizieren'] = processing.run('gdal:warpreproject', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
                
        if str(parameters[self.INPUT_ALS]) is not str(asc):
            alg_params = {
                'ALPHA_BAND': False,
                'CROP_TO_CUTLINE': True,
                'DATA_TYPE': 0,
                'EXTRA': '',
                'INPUT': str(parameters[self.INPUT_ALS]),
                'KEEP_RESOLUTION': False,
                'MASK': str(buffer),
                'MULTITHREADING': False,
                'NODATA': None,
                'OPTIONS': '',
                'SET_RESOLUTION': False,
                'SOURCE_CRS': QgsCoordinateReferenceSystem(epsg),
                'TARGET_CRS': QgsCoordinateReferenceSystem(epsg),
                'X_RESOLUTION': None,
                'Y_RESOLUTION': None,
                'OUTPUT': out
            }
            processing.run('gdal:cliprasterbymasklayer', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        del out

        #-------------------------------------------------------------------------------------------------------------------------------------------#
        ## 5. Get Cell Size and EPSG
        out = temp_path + '/' + 'dgm_orig.tif'
        raster_cs = gdal.Open(out)
        gt_cs =raster_cs.GetGeoTransform()
        proj = osr.SpatialReference(wkt=raster_cs.GetProjection()) 
        cs = gt_cs[1]
        epsg2 = proj.GetAttrValue('AUTHORITY',1)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(int(epsg2))
        logfile.write("Cellsize: " + str(cs)+"\n")
        logfile.write("EPSG: " + str(epsg2)+"\n")
        
        del raster_cs
        del gt_cs

        # 5.1. Resolution
        sel2 = self.SELECTIONS2[self.parameterAsEnum(parameters, self.RESAMPLE, context)]
            
        ## Resolution of Input-DTM
        if str(sel2) == str('1m'):
            leng_f = 1
            iter = 200
            out = temp_path + '/' + 'dgm_orig.tif'
            feedback.pushInfo(str(sel2))
            GPHmax = 80
            GPHmin = 20
            S1max = str(parameters[self.S1])
            S1max = S1max.replace(",",".")
            S1max = float(S1max)
            S2max = str(parameters[self.S2])
            S2max = S2max.replace(",",".")
            S2max = float(S2max)
            S3max = str(parameters[self.S3])
            S3max = S3max.replace(",",".")
            S3max = float(S3max)
            #imax30 = S1max + 1
            imax28 = S1max + 0.2   
            imax26 = S2max + 0.2
            imax20 = S3max + 0.2
            #imin30 = S3max - 2.5
            imin28 = S1max - 0.2     
            imin26 = S2max - 0.2
            imin20 = S3max - 0.2
                
        ## Resolution: 2 m
        if str(sel2) == str('2m'):
            leng_f = 2
            iter = 150
            out_orig = temp_path + '/' + 'dgm_orig.tif'
            out = temp_path + '/' + 'dgm_resample.tif'
            alg_params = {
                'DATA_TYPE': 0,
                'EXTRA': '',
                'INPUT': out_orig,
                'MULTITHREADING': False,
                'NODATA': None,
                'OPTIONS': '',
                'RESAMPLING': 1,
                'SOURCE_CRS': QgsCoordinateReferenceSystem('EPSG:'+str(epsg2)), #31287
                'TARGET_CRS': QgsCoordinateReferenceSystem('EPSG:'+str(epsg2)),
                'TARGET_EXTENT': None,
                'TARGET_EXTENT_CRS': None,
                'TARGET_RESOLUTION': 2,
                'OUTPUT': out
            }
            outputs['TransformierenReprojizieren'] = processing.run('gdal:warpreproject', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
            feedback.pushInfo(str(sel2))

            # Get Cellsize and EPSG code
            del proj
            del cs
            del epsg
            del srs
                
            raster_cs = gdal.Open(str(out))
            gt_cs =raster_cs.GetGeoTransform()
            proj = osr.SpatialReference(wkt=raster_cs.GetProjection()) 
            cs = gt_cs[1]
            epsg = proj.GetAttrValue('AUTHORITY',1)
            srs = osr.SpatialReference()
            srs.ImportFromEPSG(int(epsg))
            logfile.write("Cellsize: " + str(cs)+"\n")
            logfile.write("EPSG: " + str(epsg)+"\n")
        
            del raster_cs
            del gt_cs
                
            GPHmax = 80
            GPHmin = 20
            S1max = str(parameters[self.S1])
            S1max = S1max.replace(",",".")
            S1max = float(S1max)
            S2max = str(parameters[self.S2])
            S2max = S2max.replace(",",".")
            S2max = float(S2max)
            S3max = str(parameters[self.S3])
            S3max = S3max.replace(",",".")
            S3max = float(S3max)
            #imax30 = S1max + 1
            imax28 = S1max + 0.3   
            imax26 = S2max + 0.3
            imax20 = S3max + 0.3
            #imin30 = S3max - 2.5
            imin28 = S1max - 0.3     
            imin26 = S2max - 0.3
            imin20 = S3max - 0.3

        ## Resolution: 3 m
        if str(sel2) == str('3m'):
            leng_f = 3
            iter = 100
            out_orig = temp_path + '/' + 'dgm_orig.tif'
            out = temp_path + '/' + 'dgm_resample.tif'
            alg_params = {
                'DATA_TYPE': 0,
                'EXTRA': '',
                'INPUT': out_orig,
                'MULTITHREADING': False,
                'NODATA': None,
                'OPTIONS': '',
                'RESAMPLING': 1,
                'SOURCE_CRS': QgsCoordinateReferenceSystem('EPSG:'+str(epsg2)), #31287
                'TARGET_CRS': QgsCoordinateReferenceSystem('EPSG:'+str(epsg2)),
                'TARGET_EXTENT': None,
                'TARGET_EXTENT_CRS': None,
                'TARGET_RESOLUTION': 3,
                'OUTPUT': out
            }
            outputs['TransformierenReprojizieren'] = processing.run('gdal:warpreproject', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
            feedback.pushInfo(str(sel2))

            ## Cell Size and EPSG code
            del proj
            del cs
            del epsg
            del srs
                
            raster_cs = gdal.Open(str(out))
            gt_cs =raster_cs.GetGeoTransform()
            proj = osr.SpatialReference(wkt=raster_cs.GetProjection()) 
            cs = gt_cs[1]
            epsg = proj.GetAttrValue('AUTHORITY',1)
            srs = osr.SpatialReference()
            srs.ImportFromEPSG(int(epsg))
            logfile.write("Cellsize: " + str(cs)+"\n")
            logfile.write("EPSG: " + str(epsg)+"\n")
        
            del raster_cs
            del gt_cs
                
            GPHmax = 80
            GPHmin = 20
            S1max = str(parameters[self.S1])
            S1max = S1max.replace(",",".")
            S1max = float(S1max)
            S2max = str(parameters[self.S2])
            S2max = S2max.replace(",",".")
            S2max = float(S2max)
            S3max = str(parameters[self.S3])
            S3max = S3max.replace(",",".")
            S3max = float(S3max)
            #imax30 = S1max + 1
            imax28 = S1max + 0.5   
            imax26 = S2max + 0.5
            imax20 = S3max + 0.5
            #imin30 = S3max - 2.5
            imin28 = S1max - 0.5     
            imin26 = S2max - 0.5
            imin20 = S3max - 0.5
            
        ## Resolution: 5 m
        if str(sel2) == str('5m'):
            leng_f = 5
            iter = 60
            out_orig = temp_path + '/' + 'dgm_orig.tif'
            out = temp_path + '/' + 'dgm_resample.tif'
            alg_params = {
                'DATA_TYPE': 0,
                'EXTRA': '',
                'INPUT': out_orig,
                'MULTITHREADING': False,
                'NODATA': None,
                'OPTIONS': '',
                'RESAMPLING': 1,
                'SOURCE_CRS': QgsCoordinateReferenceSystem('EPSG:'+str(epsg)), #31287
                'TARGET_CRS': QgsCoordinateReferenceSystem('EPSG:'+str(epsg)),
                'TARGET_EXTENT': None,
                'TARGET_EXTENT_CRS': None,
                'TARGET_RESOLUTION': 5,
                'OUTPUT': out
            }
            outputs['TransformierenReprojizieren'] = processing.run('gdal:warpreproject', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
            feedback.pushInfo(str(sel2))

            # Cellsize and EPSG code
            del proj
            del cs
            del epsg
            del srs
                
            raster_cs = gdal.Open(str(out))
            gt_cs =raster_cs.GetGeoTransform()
            proj = osr.SpatialReference(wkt=raster_cs.GetProjection()) 
            cs = gt_cs[1]
            epsg = proj.GetAttrValue('AUTHORITY',1)
            srs = osr.SpatialReference()
            srs.ImportFromEPSG(int(epsg))
            logfile.write("Cellsize: " + str(cs)+"\n")
            logfile.write("EPSG: " + str(epsg)+"\n")
        
            del raster_cs
            del gt_cs
                
            GPHmax = 80
            GPHmin = 20
            S1max = str(parameters[self.S1])
            S1max = S1max.replace(",",".")
            S1max = float(S1max)
            S2max = str(parameters[self.S2])
            S2max = S2max.replace(",",".")
            S2max = float(S2max)
            S3max = str(parameters[self.S3])
            S3max = S3max.replace(",",".")
            S3max = float(S3max)
            #imax30 = S1max + 1
            imax28 = S1max + 0.5   
            imax26 = S2max + 0.5
            imax20 = S3max + 0.5
            #imin30 = S3max - 2.5
            imin28 = S1max - 0.5     
            imin26 = S2max - 0.5
            imin20 = S3max - 0.5
            
        logfile.write("Untersuchtes Gebiet herausgeklippt und resampled"+"\n")

        #-------------------------------------------------------------------------------------------------------------------------------------------# 
        ## 6. Read DTM as Array
        dtm_rds = gdal.Open(out)
        format = "GTiff"
        driver = gdal.GetDriverByName( format )
        band1_dtm = dtm_rds.GetRasterBand(1)
        gt = dtm_rds.GetGeoTransform()
        width = dtm_rds.RasterXSize
        height = dtm_rds.RasterYSize
        minx = gt[0]
        miny = gt[3] + width*gt[4] + height*gt[5] 
        maxx = gt[0] + width*gt[1] + height*gt[2]
        maxy = gt[3]
        extent=str(minx)+","+str(maxx)+","+str(miny)+","+str(maxy)
        extent2=(minx,maxx,miny,maxy)
        feedback.pushInfo(str(extent))
        #pix_dtm = band1_dtm.ReadAsArray() ##Error
        im_dtm = Image.open(out)
        im_flip = im_dtm.transpose(Image.FLIP_TOP_BOTTOM)
        pix_dtm = im_flip.load()
            
        #-------------------------------------------------------------------------------------------------------------------------------------------#
        ## 7. Generate profile points from breakline
        out11 = temp_path + '/' + 'line_vertices.shp'
        lines = gpd.read_file(out_hl)

        profile_points = []
        for line in lines.geometry:
            num_points = int(line.length / 1)  # Beispiel für einen Punkt alle Meter
            for distance in np.linspace(0, line.length, num_points):
                point = line.interpolate(distance)
                elevation = get_elevation(point.x, point.y, band1_dtm, gt)
                profile_points.append({'geometry': Point(point.x, point.y), 'elevation': elevation})
        gdf_points = gpd.GeoDataFrame(profile_points)
        gdf_points.to_file(out11)

        logfile.write("Profile points from breakline generated"+"\n")
            
        #-------------------------------------------------------------------------------------------------------------------------------------------#
        ## 8. Determining the direction of flow of sliding snow
        #grid_flow = np.zeros(shape=(width,height), dtype=np.float32)
            
        # 8.1. Extract coordinates
        count_sp = 0
        cell_list =[]
        coord = open(final_path + '/' + 'coord.txt', "w")
        coord.write("x;y;z;slope;len;id" +"\n")
        coord_pts = open(final_path + '/' + 'coord_pts.txt', "w")
        coord_pts.write("x;y;z;slope;len;id" +"\n")
        info = open(final_path + '/' + 'info.txt', "w")
                
        vLayer = QgsVectorLayer(out11, "layer2", "ogr")
        features = vLayer.getFeatures()
        for feature in features:
            geometry = feature.geometry()
            x, y = geometry.asPoint()
            x_save = x
            y_save = y            
         
            # Get starting pixel
            for row in range(0, width): #gdal - height
                x_neu = minx+(cs*row)-(cs/2)
                for col in range(0, height): #gdal - width
                    y_neu= miny+(cs*col)-(cs/2)
                    val_dtm = pix_dtm[row, col]
                    if x_neu > x - (cs/2) and x_neu < x + (cs/2) and y_neu > y - (cs/2) and y_neu < y + (cs/2):
                        x_start = row
                        y_start = col
                        lst = [row, col]
                        cell_list.append(lst)
                        count_sp = count_sp+1
            
        del row, col
        del x, y
        
        
        ## 8.2. Create Shapefile
        line = final_path + '/' + 'line.shp'
        driver = ogr.GetDriverByName('Esri Shapefile')
        ds = driver.CreateDataSource(line)
        layer = ds.CreateLayer('', srs, ogr.wkbLineString)
        layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
        layer.CreateField(ogr.FieldDefn('grad', ogr.OFTInteger))
        layer.CreateField(ogr.FieldDefn('len', ogr.OFTInteger))
        layer.CreateField(ogr.FieldDefn('grad_class', ogr.OFTReal))
        defn = layer.GetLayerDefn()
        feat = ogr.Feature(defn)
        id_zelle = 1
  
  
        ## 8.3. Loop through starting list
        for cell in cell_list:
            len = 0
            row = cell[0]
            col = cell[1]

            flow=(row, col) #col-1
            val_dtm = pix_dtm[row-1, col-1]
                
            x = minx+(cs*row)-(cs/2)
            y = miny+(cs*col)-(cs/2)
                
            feat.SetField('id', id_zelle)
            feat.SetField('grad', 0)
            feat.SetField('len', 0)
            geom = ogr.Geometry(ogr.wkbLineString)
            geom.AddPoint(x, y)
                
            coord.write(str(x) + ";" + str(y) + ";" + str(val_dtm) + ";" + str(0) + ";" + str(len) + ";" + str(id_zelle)+"\n")
            coord_pts.write(str(x) + ";" + str(y)+ ";" + str(val_dtm) + ";" + str(0) + ";" + str(len) + ";" + str(id_zelle) + "\n")
                               
                               
            # 8.4. 3x3 neighbourhood analysis
            start = -1 
            end = 2 

            # counter
            count = 0
            count2 = 0
                
            # loop
            while count <= iter: 
                row_l = flow[0]
                col_l = flow[1]                
                list_n =[]
                for dx2 in range(start,end):
                    for dy2 in range(start,end):
                        if (not (dx2 == dy2 == 0)):
                            ## D8- Ansatz
                            z2 = pix_dtm[row_l-1, col_l-1]
                            z_n2 = pix_dtm[row_l-1+dx2, col_l-1+dy2]
                            if dx2 == -1 and dy2 == 1 or dx2 == 1 and dy2 == 1 or dx2 == 1 and dy2 == -1 or dx2 == -1 and dy2 == -1:
                                cross2 = (math.sqrt(2)*cs)
                                diff2 = z2-z_n2
                                value_n2 = diff2/cross2
                            if dx2 == 0 and dy2 == 1 or dx2 == 1 and dy2 == 0 or dx2 == 0 and dy2 == -1 or dx2 == -1 and dy2 == 0:
                                diff2 = z2-z_n2
                                value_n2 = diff2/cs
                            list_n.append(value_n2)
                max_n = max(list_n)
                for dx3 in range(start,end):
                    for dy3 in range(start,end):
                        if (not (dx3 == dy3 == 0)):
                            z3 = pix_dtm[row_l-1, col_l-1]
                            z_n3 = pix_dtm[row_l-1+dx3, col_l-1+dy3]
                            if dx3 == -1 and dy3 == 1 or dx3 == 1 and dy3 == 1 or dx3 == 1 and dy3 == -1 or dx3 == -1 and dy3 == -1:
                                cross3 = (math.sqrt(2)*cs)
                                diff3 = z3-z_n3
                                value_n3 = diff3/cross3
                            if dx3 == 0 and dy3 == 1 or dx3 == 1 and dy3 == 0 or dx3 == 0 and dy3 == -1 or dx3 == -1 and dy3 == 0:
                                diff3 = z3-z_n3
                                value_n3 = diff3/cs
                                    
                            if max_n < 0: # Stops the loop if the pixel is a sink
                                count = iter
                                del flow
                                flow =(row_l+dx3, col_l+dy3)
                                x_flow2 = minx+(cs*(row_l+dx3))-(cs/2)
                                y_flow2= miny+(cs*(col_l+dy3))-(cs/2) 
                                z_flow2 = pix_dtm[row_l-1+dx3, col_l-1+dy3]
                                coord_pts.write(str(x_flow2) + ";" + str(y_flow2)+ ";" + str(z_n3) + ";" + str(round(grad2,2)) + ";" + str(l2) + ";" + str(id_zelle) + "\n")
                                coord.write(str(x_flow2) + ";" + str(y_flow2)+ ";" + str(z_n3) + ";" + str(2) + ";" + str(0) + ";" + str(id_zelle) + "\n")
                                
                            if max_n > 0:
                                if max_n == value_n3:
                                    del flow
                                    flow =(row_l+dx3, col_l+dy3)
                                    x_flow2 = minx+(cs*(row_l+dx3))-(cs/2)
                                    y_flow2= miny+(cs*(col_l+dy3))-(cs/2) 
                                    z_flow2 = pix_dtm[row_l-1+dx3, col_l-1+dy3]
                                        
                                    if z_n3 > 0:                                        
                                        h_diff2 = round((val_dtm - z_n3),2)
                                        laenge = round((math.sqrt((((dx3*cs)-0)**2)+(((dy3*cs)-0)**2))),2)
                                        len = len + laenge
                                        l2 = len
                                        rad2 = math.atan(h_diff2/l2)
                                        grad2 = int(abs(math.degrees(rad2)))
                                            
 
                                        # 8.5. 30 degree points
                                        if h_diff2 >= GPHmax:
                                            dx_save = dx3
                                            dy_save = dy3
                                            coord.write(str(x_flow2) + ";" + str(y_flow2)+ ";" + str(z_n3) + ";" + str(grad2) + ";" + str(l2) + ";" + str(id_zelle) + "\n")
                                            geom.AddPoint(x_flow2, y_flow2)
                                                        
                                        if h_diff2 <= GPHmax and h_diff2 >= GPHmin:
                                            dx_save = dx3
                                            dy_save = dy3
                                            coord_pts.write(str(x_flow2) + ";" + str(y_flow2)+ ";" + str(z_n3) + ";" + str(round(grad2,2)) + ";" + str(l2) + ";" + str(id_zelle) + "\n")
                                            coord.write(str(x_flow2) + ";" + str(y_flow2)+ ";" + str(z_n3) + ";" + str(grad2) + ";" + str(l2) + ";" + str(id_zelle) + "\n")
                                            count = iter
                                            x_end = x_flow2
                                            y_end = y_flow2
                                            row_dir = row_l+dx_save
                                            col_dir = col_l+dy_save
                                            len_dir = len
                                            geom.AddPoint(x_end, y_end)
                                            feat.SetGeometry(geom)
                                            layer.CreateFeature(feat)
                                            
                                            row_dir = row_dir+dx_save
                                            col_dir = col_dir+dy_save
                                            count2 = count2+1
                                                            
                count = count+1
            id_zelle = id_zelle + 1 
 
        feat = geom = None 
        ds = layer = feat = geom = None
        logfile.write("Possible snow mass direction lines created"+"\n")

        #-------------------------------------------------------------------------------------------------------------------------------------------#
        ## 9. Extend lines
        RW_final = final_path + '/' + 'RW_final.shp'
        alg_params = {
            'END_DISTANCE': 200,
            'INPUT': line,
            'START_DISTANCE': 0,
            'OUTPUT': RW_final
        }
        outputs['LinienErweitern'] = processing.run('native:extendlines', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
     
        #-------------------------------------------------------------------------------------------------------------------------------------------#
        ## 10. Points along a geometry
        pts_temp = temp_path + '/' + 'pts_temp.shp'
        alg_params = {
            'DISTANCE': cs,
            'END_OFFSET': 0,
            'INPUT': RW_final,
            'START_OFFSET': 0,
            'OUTPUT': pts_temp
        }
        outputs['PunkteEntlangEinerGeometrie'] = processing.run('native:pointsalonglines', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        
        #-------------------------------------------------------------------------------------------------------------------------------------------#
        # 11. Scan grid values
        pts_Z = temp_path + '/' + 'pts_Z.shp'
        alg_params = {
            'COLUMN_PREFIX': 'rvalue',
            'INPUT': pts_temp,
            'RASTERCOPY': out,
            'OUTPUT': pts_Z
        }
        outputs['RasterwerteAbtasten'] = processing.run('qgis:rastersampling', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        
        #-------------------------------------------------------------------------------------------------------------------------------------------#
        ## 12. Link attributes by field value 1
        line_out_Z = temp_path + '/' + 'line_out_Z.shp'
        alg_params = {
            'DISCARD_NONMATCHING': False,
            'FIELD': 'id',
            'FIELDS_TO_COPY': None,
            'FIELD_2': 'id',
            'INPUT': RW_final,
            'INPUT_2': pts_Z,
            'METHOD': 1,
            'PREFIX': '',
            'OUTPUT': line_out_Z#['Save']
        }
        outputs['AttributeNachFeldwertVerknpfen'] = processing.run('native:joinattributestable', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        
        #-------------------------------------------------------------------------------------------------------------------------------------------#
        ## 13. Link attributes by field value 2
        pts_out_Z = temp_path + '/' + 'pts_out_Z.shp'
        alg_params = {
            'DISCARD_NONMATCHING': False,
            'FIELD': 'id',
            'FIELDS_TO_COPY': None,
            'FIELD_2': 'id',
            'INPUT': pts_Z,
            'INPUT_2': line_out_Z,
            'METHOD': 1,
            'PREFIX': '',
            'OUTPUT': pts_out_Z#['Save']
        }
        outputs['AttributeNachFeldwertVerknpfen'] = processing.run('native:joinattributestable', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        
        #-------------------------------------------------------------------------------------------------------------------------------------------#
        ## 14. compute height difference to start cell
        pts_out_calc = temp_path + '/' + 'pts_out_calc.shp'
        alg_params = {
            'FIELD_LENGTH': 10,
            'FIELD_NAME': 'HDIFF',
            'FIELD_PRECISION': 5,
            'FIELD_TYPE': 0,
            'FORMULA': '\"rvalue1_2\" - \"rvalue1\" ',
            'INPUT': pts_out_Z,
            'NEW_FIELD': True,
            'OUTPUT': pts_out_calc
        }
        outputs['Feldrechner'] = processing.run('qgis:fieldcalculator', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        
        #-------------------------------------------------------------------------------------------------------------------------------------------#
        ## 15. compute atan of height difference and distance
        pts_out_calc_atan = temp_path + '/' + 'pts_out_calc_atan.shp'
        alg_params = {
            'FIELD_LENGTH': 10,
            'FIELD_NAME': 'atan',
            'FIELD_PRECISION': 5,
            'FIELD_TYPE': 0,
            'FORMULA': ' atan(\"HDIFF\" / \"distance\") ',
            'INPUT': pts_out_calc,
            'NEW_FIELD': True,
            'OUTPUT': pts_out_calc_atan
        }
        outputs['Feldrechner'] = processing.run('qgis:fieldcalculator', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        
        #-------------------------------------------------------------------------------------------------------------------------------------------#
        ## 16. compute degree from atan
        pts_out_calc_final = temp_path + '/' + 'pts_out_calc_final.shp'
        alg_params = {
            'FIELD_LENGTH': 10,
            'FIELD_NAME': 'GRD',
            'FIELD_PRECISION': 2,
            'FIELD_TYPE': 0,
            'FORMULA': ' degrees(\"atan\") ',
            'INPUT': pts_out_calc_atan,
            'NEW_FIELD': True,
            'OUTPUT': pts_out_calc_final
        }
        outputs['Feldrechner'] = processing.run('qgis:fieldcalculator', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        

        # 16.1. Search Cursor
        vLayer = QgsVectorLayer(pts_out_calc_final, "layer1", "ogr")
        features = vLayer.getFeatures()
        attr=[]
        fields = vLayer.fields()
        
        vLayer.startEditing()
        for feature in features:
            ## Degree
            grad = feature['GRD']
            if grad <= imax28 and grad >= imin28:
                #feature.changeAttribute(feature[20], 28)
                feature['grad_class'] = S1max
                vLayer.updateFeature(feature)
            if grad <= imax26 and grad >= imin26:
                feature['grad_class'] = S2max
                vLayer.updateFeature(feature)
            if grad <= imax20 and grad >= imin20:
                feature['grad_class'] = S3max
                vLayer.updateFeature(feature)
           
        vLayer.commitChanges()
        

        # 16.2. Search Cursor
        vLayer = QgsVectorLayer(pts_out_calc_final, "layer1", "ogr")
        features = vLayer.getFeatures()
        attr=[]
        fields = vLayer.fields()
        
        vLayer.startEditing()
        for feature in features:
            ## Degree
            grad = feature['distance']
            if grad == 0:
                feature['grad_class'] = 0
                vLayer.updateFeature(feature)
       
           
        vLayer.commitChanges()

        
        # 16.3. Delete duplicates by attribute
        pts = final_path + '/' + 'pts.shp'
        alg_params = {
            'FIELDS': ['id','grad_class'],
            'INPUT': pts_out_calc_final,
            'OUTPUT': pts
        }
        outputs['DuplikateNachAttributLschen'] = processing.run('native:removeduplicatesbyattribute', alg_params, context=context, feedback=feedback, is_child_algorithm=True)


        # 16.4. Points along a geometry
        out11_n = temp_path + '/' + 'line_vertices_red.shp'
        alg_params = {
            'DISTANCE': 10,
            'END_OFFSET': 0,
            'INPUT': out_hl,
            'START_OFFSET': 0,
            'OUTPUT': out11_n
        }
        outputs['PunkteEntlangEinerGeometrie'] = processing.run('native:pointsalonglines', alg_params, context=context, feedback=feedback, is_child_algorithm=True)


        #-------------------------------------------------------------------------------------------------------------------------------------------#
        ## 17. Determine the entire flow path
        # 17.1. Extract coordinates
        count_sp = 0
        coord_start = open(final_path + '/' + 'coord_start.txt', "w")
        coord_start.write("x;y" +"\n")
        cell_start = open(final_path + '/' + 'cell_start.txt', "w")
        cell_list =[]
        coord = open(final_path + '/' + 'coord_path.txt', "w")
        coord.write("x;y;z;slope;id" +"\n")
                
        vLayer = QgsVectorLayer(out11_n, "layer2", "ogr")
        features = vLayer.getFeatures()
        for feature in features:
            geometry = feature.geometry()
            x, y = geometry.asPoint()
            x_save = x
            y_save = y
                

            # 17.2. Extract starting pixel
            for row in range(0, width): #gdal - height
                x_neu = minx+(cs*row)-(cs/2)
                for col in range(0, height): #gdal - width
                    y_neu= miny+(cs*col)-(cs/2)
                    val_dtm = pix_dtm[row, col]
                    if x_neu > x - (cs/2) and x_neu < x + (cs/2) and y_neu > y - (cs/2) and y_neu < y + (cs/2):
                        x_start = row
                        y_start = col
                        coord_start.write(str(x_neu)+";"+str(y_neu)+"\n")
                        cell_start.write(str(row)+","+str(col)+"\n")
                        lst = [row, col]
                        cell_list.append(lst)
                        count_sp = count_sp+1
            
        del cell_start
        del row, col
        del x, y
            

        # 17.3. Create shapefile
        cell_start = open(final_path + '/' + 'cell_start.txt', "r")
            
        line_path = final_path + '/' + 'steepest_path.shp'
        driver_path = ogr.GetDriverByName('Esri Shapefile')
        ds_path = driver_path.CreateDataSource(line_path)
        layer_path = ds_path.CreateLayer('', srs, ogr.wkbLineString)
        layer_path.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
        defn_path = layer_path.GetLayerDefn()
        feat_path = ogr.Feature(defn_path)

        id_zelle = 1
            

        # 17.4. Loop through starting list
        for cell in cell_list:
            row = cell[0]
            col = cell[1]

            flow=(row, col) #col-1
            val_dtm = pix_dtm[row-1, col-1]
            z = val_dtm
                
            x = minx+(cs*row)-(cs/2)
            y = miny+(cs*col)-(cs/2)
                
            feat_path.SetField('id', id_zelle)
            feat_path.SetField('grad', 0)
            geom_path = ogr.Geometry(ogr.wkbLineString)
            geom_path.AddPoint(x, y)
                
            coord.write(str(x) + ";" + str(y) + ";" + str(z) + ";" + str(0) + ";" + str(id_zelle)+"\n")
                

            # 17.5. 3x3 neighborhooud analysis
            start = -1
            end = 2

            # Counter
            count = 0
                
            # Loop
            while count <= iter: 
                row_l = flow[0]
                col_l = flow[1]                
                list_n =[]
                for dx2 in range(start,end):
                    for dy2 in range(start,end):
                        if (not (dx2 == dy2 == 0)):
                            ## D8-approach
                            z2 = pix_dtm[row_l-1, col_l-1]
                            z_n2 = pix_dtm[row_l-1+dx2, col_l-1+dy2]
                            if dx2 == -1 and dy2 == 1 or dx2 == 1 and dy2 == 1 or dx2 == 1 and dy2 == -1 or dx2 == -1 and dy2 == -1:
                                cross2 = (math.sqrt(2)*cs)
                                diff2 = z2-z_n2
                                value_n2 = diff2/cross2
                            if dx2 == 0 and dy2 == 1 or dx2 == 1 and dy2 == 0 or dx2 == 0 and dy2 == -1 or dx2 == -1 and dy2 == 0:
                                diff2 = z2-z_n2
                                value_n2 = diff2/cs 
                            list_n.append(value_n2)
                max_n = max(list_n)
                for dx3 in range(start,end):
                    for dy3 in range(start,end):
                        if (not (dx3 == dy3 == 0)):
                            z3 = pix_dtm[row_l-1, col_l-1]
                            z_n3 = pix_dtm[row_l-1+dx3, col_l-1+dy3]
                            if dx3 == -1 and dy3 == 1 or dx3 == 1 and dy3 == 1 or dx3 == 1 and dy3 == -1 or dx3 == -1 and dy3 == -1:
                                cross3 = (math.sqrt(2)*cs)
                                diff3 = z3-z_n3
                                value_n3 = diff3/cross3
                            if dx3 == 0 and dy3 == 1 or dx3 == 1 and dy3 == 0 or dx3 == 0 and dy3 == -1 or dx3 == -1 and dy3 == 0:
                                diff3 = z3-z_n3
                                value_n3 = diff3/cs
                                
                            if max_n < 0: # Stops the loop if the pixel is a sink
                                count = iter
                                del flow
                                flow =(row_l+dx3, col_l+dy3)
                                x_flow2 = minx+(cs*(row_l+dx3))-(cs/2)
                                y_flow2= miny+(cs*(col_l+dy3))-(cs/2) 
                                z_flow2 = pix_dtm[row_l-1+dx3, col_l-1+dy3]
                                coord.write(str(x_flow2) + ";" + str(y_flow2) + ";" + str(z_flow2) + ";" + str(0) + ";" + str(id_zelle)+"\n")
                                
                            if max_n > 0:
                                if max_n == value_n3:
                                    if z_n3 > 0:
                                        x_flow2 = minx+(cs*(row_l+dx3))-(cs/2)
                                        y_flow2= miny+(cs*(col_l+dy3))-(cs/2)
                                        z_flow2 = pix_dtm[row_l-1+dx3, col_l-1+dy3]
                                        coord.write(str(x_flow2) + ";" + str(y_flow2) + ";" + str(z_flow2) + ";" + str(0) + ";" + str(id_zelle)+"\n")
                                                        
                                        geom_path.AddPoint(x_flow2, y_flow2)
                                        feat_path.SetGeometry(geom_path)
                                        layer_path.CreateFeature(feat_path)
                                            
                                        del flow
                                        flow =(row_l+dx3, col_l+dy3)
                count = count+1
            id_zelle = id_zelle + 1    

        feat_path = geom_path = None 
        ds_path = layer_path = feat_path = geom_path = None
        logfile.write("Complete flow path was determined"+"\n") 
        

        #-------------------------------------------------------------------------------------------------------------------------------------------#
        #-------------------------------------------------------------------------------------------------------------------------------------------#
        ## Visualization in QGIS
        # Add Layers to Map
        root = QgsProject.instance().layerTreeRoot()
        mygroup_name = "Snowslide-Tool - Ergebnisse"
        mygroup = root.findGroup(mygroup_name)

        # If the group doesn't exist, create it
        if mygroup is None:
            mygroup = root.insertGroup(0, "Snowslide-Tool - Ergebnisse")
        else:
            clone = mygroup.clone()
            parent = mygroup.parent() if mygroup.parent() else root
            parent.removeChildNode(mygroup)
            newGroup = parent.insertChildNode(0, clone)
            mygroup = newGroup or clone
        
        # Add Impact-Points
        layer1 = QgsVectorLayer(pts, "Impact-Points", "ogr")

        # Symbology
        symbol = QgsSymbol.defaultSymbol(layer1.geometryType())    
        
        cat0_wert = 0
        cat0_symbol = QgsMarkerSymbol.createSimple({'color': 'black', 'size': '1', 'outline_color': 'black'})
        cat0 = QgsRendererCategory(cat0_wert, cat0_symbol, 'Start')

        cat1_wert = S1max
        cat1_symbol = QgsMarkerSymbol.createSimple({'color': '#0550f2', 'size': '1.5', 'outline_color': 'black'})
        cat1 = QgsRendererCategory(cat1_wert, cat1_symbol, str(S1max) + ' Grad')
        
        cat2_wert = S2max
        cat2_symbol = QgsMarkerSymbol.createSimple({'color': '#63bef6', 'size': '1.5', 'outline_color': 'black'})
        cat2 = QgsRendererCategory(cat2_wert, cat2_symbol, str(S2max) + ' Grad')
        
        cat3_wert = S3max
        cat3_symbol = QgsMarkerSymbol.createSimple({'color': '#ffed1f', 'size': '1.5', 'outline_color': 'black'})
        cat3 = QgsRendererCategory(cat3_wert, cat3_symbol, str(S3max) + ' Grad')
        
        categories =(cat0,cat1,cat2,cat3)
        renderer = QgsCategorizedSymbolRenderer('grad_class', categories)
        if renderer is not None:
            layer1.setRenderer(renderer)
        layer1.triggerRepaint()
        
        QgsProject.instance().addMapLayer(layer1, False)
        mygroup.addLayer(layer1)

        ## Add flow paths
        layer2 = QgsVectorLayer(line_path, "Fliesspfade", "ogr")
        label = QgsPalLayerSettings()
        labeler2 = QgsVectorLayerSimpleLabeling(label)
        layer2.setLabelsEnabled(True)
        layer2.setLabeling(labeler2)
        layer2.renderer().symbol().setColor(QColor("#33b5f5"))
        layer2.triggerRepaint()
        QgsProject.instance().addMapLayer(layer2, False)
        mygroup.addLayer(layer2)
        
        ## Add directional routes
        layer3 = QgsVectorLayer(line, "Richtungswege", "ogr")
        labeler3 = QgsVectorLayerSimpleLabeling(label)
        layer3.setLabelsEnabled(True)
        layer3.setLabeling(labeler3)
        symbol3 = QgsLineSymbol.createSimple({'line_style': 'dash', 'color': 'black'})
        layer3.renderer().setSymbol(symbol3)
        #layer3.renderer().symbol().setColor(QColor("Black"))
        layer3.triggerRepaint()
        QgsProject.instance().addMapLayer(layer3, False)
        mygroup.addLayer(layer3)



        #-------------------------------------------------------------------------------------------------------------------------------------------#
        # THE END        
        logfile.write("Layers successfully added to map"+"\n") 
        logfile.write("SnowslideTool finished"+"\n") 
        timestamp = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
        logfile.write(str(timestamp)+"\n")
        
        
        outputs['LastStep'] = line_path
        results['Work finished'] = outputs['LastStep']
        return results
            


