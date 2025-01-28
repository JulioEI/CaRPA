# miniscope_analysis
Scripts to process and analyze miniature microscope data.

##License
*	Don't distribute code outside the Schnitzer Lab without asking first.

##Usage

*	Run __loadRepoFunctions.m__ before using functions in the directory. This adds all directories and sub-directories to Matlab path.

##Top-level directories
Below are a list of the top-level directories and what types of functions or files should be placed there.

*	__cellmax__ - Contains code for running CELLMax algorithm for extracting traces from movies.
*	__file\_exchange__ - Contains any outside code from Matlab's File Exchange that are dependencies in repository functions.
*	__hdf5__ - Functions concerned with HDF5 input/output.
*	__image__ - Functions concerned with processing images (or [x y] matrices).
*	__io__ - Contains functions concerned with file or function input-output.
*	__motion\_correction__ - Functions concerned with motion correction.
*	__movie\_processing__ - Functions concerned with manipulating movies, e.g. removing noise.
*	__settings__ - Functions concerned with settings for other functions.
*	__view__ - Functions concerned with displaying data or information to the user, no real processing of data.