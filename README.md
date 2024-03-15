<h1>pyAlpineRisk<em> - natural hazard management of alpine torrential catchments</em></h1>

<p>At a time when severe weather events are becoming more frequent, the availability of detailed and high-resolution digital terrain models is key to the management of natural hazards.  For this reason, the author, in close collaboration with die.wildbach, developed a set of tools (pyAlpineRisk) that can be applied in the field of natural hazard management and adapted to other scientific enquiries. By using simple terrain models (DTM) and various spatial analysis methods, they were able to collect important data that can be used to address landform issues and support natural hazard management.
All pyAlpineRisk tools are under constant development.
</p>
 
<h2>Snowslide Tool</h2>
<u>QGIS-SnowslideTool: GIS-supported range determination of snowslide processes (flat-rate slope method)</u>
</p>

<p>The QGIS-SnowslideTool uses high-resolution terrain models and GIS techniques to automatically calculate and display precise slope and range estimates in the terrain. It improves the assessment of potential natural hazards and reduces the amount of manual work involved. The tool is particularly suitable for mapping areas of snowslides in areas below 1500 m above sea level. It automatically determines flow paths and slope direction lines and visualizes them within the QGIS project. The automated terrain analysis with the QGIS-SnowslideTool is a promising approach for precise and efficient data analysis and interpretation from remote sensing data.</p>

<h2>Installation/Application</h2>
<p>The latest release is written for PyQGIS 3.28 and can be used by installing QGIS 3.28 or above.

To install QGIS tools developed for QGIS 3.x, copy them into
~/AppData/Roaming/QGIS/QGIS3/profiles/default/processing/scripts or in the upper part of the toolbox dialog you can add the scripts with ![mIconPythonFile](https://user-images.githubusercontent.com/52344347/136413201-b4a1f7d3-4053-4aa6-b11c-9433ae617057.png) Scripts - Add Script to Toolbox ...

After that the tools can be found in the QGIS "Processing Toolbox" - Scripts</p>

<h2>Case Examples</h2>
<p>
<i lang="id">Input Parameters:</i>

![B1](https://github.com/pyAlpineRisk/SnowslideTool/assets/52344347/aef70ba7-f783-42a4-8eb2-e7fbc70e2be8)
<p>The described input form is a conventional QGIS tool interface that provides a user-friendly way to edit and analyze geospatial data. In addition to the input parameters, the form also includes a help window that provides the user with a brief description of the steps to be performed. The input parameters are divided into different categories. The first category is the input geometry, which can be either a line or an area. In the second category, there are three options for slope values that can be used to approximate the range. The third category is the terrain model, which forms the basis for the calculation. The resolution is another important parameter that affects the accuracy of the results. The output folder for temporary and final records allows the user to save the calculation results. Upon completion of the calculation, the results are immediately loaded into the current QGIS interface and displayed according to a predefined symbology. The results include the hit points in various slope levels, the starting position as well as direction paths and flow paths.</p>

<i lang="id">Results:</i>

![B2](https://github.com/pyAlpineRisk/SnowslideTool/assets/52344347/a2b977a9-00b7-4368-bbc2-5abfd08bbce6)

![B3](https://github.com/pyAlpineRisk/SnowslideTool/assets/52344347/7166cfeb-bcd3-4e06-b5d1-c2a3f1296ace)
