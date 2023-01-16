-----------------------------------
# CaRPA QUICK TUTORIAL
-----------------------------------

Welcome to CaRPA (Calcium Recordings Processing and Analysis). Here is a quick tutorial and some tips to get you started:

## 1. Starting CaRPA:
	
To use CaRPA, first make sure the CARPA folder is in MATLAB's path.
 
 <img src="CARPA/_readmePics/path.png" alt="Path" width="500"/>
 
Next, create an instance of the CaRPA object:

```
carp = carpa;
```
A menu will show up, listing several options:

<img src="CARPA/_readmePics/options.png" alt="Options" width="500"/>


### 1.1. Work on an existing animal:
		
To pick up the processing of an existing folder. The file structure that CaRPA works with is a root folder on top, then the folders for individual days/animals/experiment sessions, with the files for these sessions (e.g., the h5 files) inside.


| 	      | ROOT FOLDER | ANIMAL/DAY/EXPERIMENT FOLDER | SESSION FILES |
| ----------- | ----------- | ----------- | ----------- |
|E.g. 	   | Mice 50xx	   | Mouse-5012-20150228-eightPorts 	        | concat_recording_20150228_114834.h5|
|	   |		   |  					        | concat_recording_20150228_120021.h5|
|     	   |               | Mouse-5012-20150228-linearTrack		| concat_recording_20150228_192756.h5|
|	   |		   | Mouse-5013-20150301-eightPorts 	        | concat_recording_20150301_114211.h5|
|	   |	           |					        | concat_recording_20150301_115923.h5|

		
CaRPA will ask for the root folder. It is only possible to use an instance of CaRPA on one root folder at a time. 

**IMPORTANT:**

> When selecting the root folder, make sure the name of the folder you have selected appears in the window as show below. 

 <img src="CARPA/_readmePics/make_sure.png" alt="Make Sure" width="500"/>

**IMPORTANT:**

>This exact syntax of the folder name is important for CaRPA to understand the different fields:  Mouse-XXXX-YYYYMMDD-EXPERIMENT. e.g.     Mouse-5012-20150228-eightPorts 

**IMPORTANT:**

> If a given day has two different experiments, they should go in different folders (see example above).    To change the experiment names in a convenient name, use the command:
```
carp.setExperimentNames;
```


### 1.2. Create file structure for a new animal

Creates the file structure described in the previous point from a folder with several calcium_recording files. For this to work properly, the first number sequence in the file name should correspond to the day of the recording, e.g.  concat_recording_20150228_192756.h5

**IMPORTANT:**

> In CaRPA, spatially downsampled non-processed calcium files are called concat_recording_*YYYYMMDD*_*hhmmss*.h5. 


### 1.3. Download files from server

Allows to download one or more sessions of an animal directly from the cluster. The path of the server is defined by the archiveRawCaPath property. To change it use the command:

```
carp.archiveRawCaPath = ‘‘C://new_path’’;
```

It is assumed that in the path folder there are folders for the different animals, and inside each animal folder there are folders for the different days, in a similar manner as the folder structure specified in 1.1. 

This function will also decompress .raw files using the Inscopix Decompressor. For this to work the nVistaHD decompressor must be installed on the system, the file raw2xxx.py must exist in the nVistaHD directory, and the file InscopixDecompress.bat must be visible by Matlab. The InscopixDecompress.bat file contains the installation directory of nVistaHD, which is by default: C:\Program Files\Inscopix\nVistaHD\.

Finally this function will spatially downsample calcium files which do not begin by concat_recording*. It will also spatially downsample newly decompressed files. The amount of spatial downsample is defined by the propriety spatialDS, and can be changed by using the command:

```
carp.spatialDS = 4;
```

**IMPORTANT**

>The custom interface to connect to the server is mysftp, which uses the SSH2 library. Therefore, this library and the custom interface must be in the matlab path for this function to work. Besides, a username and password must be provided before connecting to the server.

## 2. The folder structure

CaRPA organizes your files in a folder structure. This is a list of the files CaRPA understands:

| File | Description | Syntax |
|--------|--------|--------|
| concat | Unprocessed, spatially downsampled ca files. | concat_recording_*YYYYMMDD*_*hhmmss*.h5 |
| dfof | Processed ca files. | Mouse-*XXXX*-*YYYYMMDD*_*hhmmss*&*hhmmss*&...*-EXPERIMENT*-dfof.h5 |
| downsample | Temporally downsampled processed ca files. | Mouse-*XXXX*-*YYYYMMDD*_*hhmmss*&*hhmmss*&...*-EXPERIMENT*-dfof-downsampled.h5 |
| analysis | .mat files with extracted cells and traces. | Mouse-*XXXX*-*YYYYMMDD*_*hhmmss*&*hhmmss*&...*-EXPERIMENT*-emAnalysis.mat OR Mouse-*XXXX*-*YYYYMMDD*_*hhmmss*&*hhmmss*&...*-EXPERIMENT*-pcaicaAnalysis.mat |
| decisions | .mat files with cell decisions.| Mouse-*XXXX*-*YYYYMMDD*_*hhmmss*&*hhmmss*&...*-EXPERIMENT*-emAnalysisSorted.mat OR Mouse-*XXXX*-*YYYYMMDD*_*hhmmss*&*hhmmss*&...*EXPERIMENT*-ICDecisions.mat |
| behavior | .avi files with animal behavior. | .avi |
| logs | .log files with recording data. | .txt OR  .html |
| tracesEvents | .mat files with the final processed traces and position data. | Mouse-*XXXX*-*YYYYMMDD*_*hhmmss*&*hhmmss*&...*-EXPERIMENT*-TracesAndEvents.mat | 
| other | Other files | * |


The folder structure is accessed using the property folderStruct, and it is updated every time it is modified. However, you can manually update the folder structure by using the command:

```
carp.buildFolderStructure;
```

## 3. Processing the data

By default CaRPA tracks the processing stage of the different sessions and only shows and applies the selected processing stages to files which have not been yet processed on those stages. To change this behavior and view/process all available sessions, regardless of preprocessing stage, change the propierty showProcessed to 1:

```
carp.showProcessed = 1; (Not recommended)
```

To start the processing, type the following command into the console:

```
carp.menu;
```

This will bring up a menu with different processing options. Select one or more options and click ok.

 <img src="CARPA/_readmePics/processing_options.png" alt="Processing Options" width="500"/>

This will bring up a menu with different processing stages. Select one or more options and click ok.

If a single processing stage is selected a menu will pop up next, which allows to process all available sessions or choose specific ones. 

 <img src="CARPA/_readmePics/single_multiple.png" alt="Single/Multiple" width="500"/>

If more than one stage is selected all the corresponding sessions will be processed.

## 3.1 Preprocess

Lorem ipsum dolor sit amet, consectetur adipiscing elit. In eget purus urna. Integer vel viverra neque. Nullam felis massa, gravida tincidunt erat quis, porta elementum tellus. Nulla sit amet bibendum purus, vitae ultrices ipsum. In condimentum pellentesque nunc, a semper augue maximus sit amet. Suspendisse mauris leo, porta a tellus quis, mattis rutrum justo. Aenean vel vehicula mauris. Etiam pretium in felis sed porttitor. Suspendisse tincidunt turpis imperdiet faucibus laoreet. Praesent elementum erat vitae augue ultrices condimentum. 

## 3.2 Extract
Proin fringilla consectetur scelerisque. In hac habitasse platea dictumst. Class aptent taciti sociosqu ad litora torquent per conubia nostra, per inceptos himenaeos. Nam vulputate orci eget ex iaculis, id volutpat lectus consequat. Etiam sit amet nunc consequat, gravida velit sit amet, gravida felis. Etiam viverra nisi ut lobortis condimentum. Duis varius justo quam, et ultrices nisi viverra at. Curabitur sit amet pretium erat. Pellentesque ut arcu purus. 

## 3.3 Filter

Donec tristique condimentum metus. Proin volutpat consectetur dignissim. Vestibulum consequat pharetra dolor in suscipit. Nam accumsan nisl et sem eleifend, eget dapibus arcu iaculis. Donec semper rhoncus massa, cursus condimentum ipsum mattis vel. Integer imperdiet convallis lectus vel laoreet. Aliquam dignissim, mi nec vestibulum dapibus, tellus elit egestas purus, eget malesuada ipsum sapien vitae erat. Vestibulum at nunc porttitor, vehicula lorem vestibulum, aliquam lorem. Interdum et malesuada fames ac ante ipsum primis in faucibus. Nam non bibendum augue, auctor scelerisque nisi. Vivamus nibh ligula, tincidunt quis vestibulum ac, semper vel turpis. Integer malesuada massa eu purus dapibus, id pretium ligula elementum. Fusce interdum, quam ac facilisis imperdiet, dui ligula tristique lorem, posuere vestibulum velit est vel ex. 

 <img src="CARPA/_readmePics/filter.png" alt="Filter" width="500"/>

## 3.4 Postprocess

Interdum et malesuada fames ac ante ipsum primis in faucibus. Aliquam mi nunc, laoreet at sem id, ultricies porta tellus. Sed vitae suscipit ante. Nunc placerat porta ante. Etiam vulputate volutpat turpis, quis convallis est tempor eget. Sed vel finibus est. Sed ut ipsum tortor. Aliquam ac lacus sem. Morbi vitae maximus dui. In hac habitasse platea dictumst. Pellentesque accumsan ipsum erat, tristique fringilla justo tincidunt ut. Donec congue hendrerit nisl, ac lobortis leo tincidunt nec. Quisque feugiat diam et mollis tempor. Sed a interdum orci, mattis venenatis tortor. Fusce ultricies dui non nunc malesuada, vitae convallis erat viverra. Fusce aliquam mauris nulla, sed vulputate mi euismod non. 

## 4. Carpa Outputs

| Variable Name | Dimensions | Description|
|--------------|--------------|--------------|
| rawProb | time x neurons  | Maximum likelihood dfof of each neuron, given by the CELLMAX algorithm|
| rawTraces| time x neurons | DFOF value of each neuron |
| position | time x coordinates | Detected tracker position. |
| velocity | time | Absolute velocity |
| tresholdEvents | time x neurons | Whether there is a spike at each time given by the thresholding algorithm method developed by Biafra Bahanonu: \miniscope_analysis-bahanonu-miniscopeAnalysisClass_18_08_17\signal_processing\computeSignalPeaks.m | 
|spikeDeconv | time x neurons | Whether there is a spike at each time given by the spike deconvolution algorithm from Pnevmatikakis et al 2013. Bayesian spike inference from calcium imaging data https://github.com/zhoupc/OASIS_matlab |
| spikeML | time x neurons | Whether there is a spike at each time given by the spikeML algorithm from : Deneux T, Kaszas A, Szalay G, Katona G, Lakner T, Grinvald A, et al. Accurate spike estimation from noisy calcium signals for ultrafast three-dimensional imaging of large neuronal populations in vivo. 2016 https://github.com/MLspike/spikes | 
| spikeDeconvTrace | time x neurons | Denoised df/f traces given by the spike deconvolution algorithm |
| cellAnatomicLocat | neurons x coordinate | Location of each neuron centroid | 


