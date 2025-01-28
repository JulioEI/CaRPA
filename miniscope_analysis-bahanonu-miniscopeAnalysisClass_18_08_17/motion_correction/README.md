# motion correction
Scripts to register movies to remove motion.

##Turboreg

*	Below is an example usage of `turboregMovie`.

```Matlab
ioptions.inputDatasetName = '/1';
ioptions.turboregRotation = 1;
ioptions.RegisType = 1;
ioptions.parallel = 1;
ioptions.meanSubtract = 1;
ioptions.normalizeType = 'divideByLowpass';
ioptions.registrationFxn = 'transfturboreg';
ioptions.normalizeBeforeRegister = 'imagejFFT';
ioptions.imagejFFTLarge = 10000;
ioptions.imagejFFTSmall = 80;
ioptions.saveNormalizeBeforeRegister = [];
ioptions.cropCoords = [];
ioptions.closeMatlabPool = 0;
ioptions.refFrame = 1;
ioptions.refFrameMatrix = [];
regMovie = turboregMovie('pathToDir\filename.h5','options',ioptions);
% OR
regMovie = turboregMovie(inputMovieMatrix,'options',ioptions);
```