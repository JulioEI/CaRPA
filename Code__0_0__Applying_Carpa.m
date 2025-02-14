%% Code_0_0_  Applying Carpa for pre - to postprocessing stages
close all; clear; clc;
disp('  @ @ @ There is ');
disp(repmat(' ',[1,100]));
disp([repmat('*',[1,35]) , '        Applying carpa       ' , repmat('*',[1,35])]);
disp(repmat(' ',[1,100]));
disp(' This is the code for calling Carpa and applying the pre to post processing stages of data analyses ')
disp(' * Each time you apply the POSTPROCESS make sure that the threshold of the tracker is set correctly  * ')
disp(' *  It is defined in the function "getMouseTrajectory" --> options.backgroundTreshold * ')
disp(' *  So far, the best practice is:')
disp('    Background threshold = -0.8 for Global Remaping , and');
disp('    Background t= -0.9 for Novel Object Location ');
disp(' ** When applying POSTPROCESS make sure that the name of xml files start with the Animnals Numbers, and end with "-log"');
disp(' *** Also, the timing indicated in the behavior video names should mathc exactly the log files name');
disp(' **** When applying POSTPROCESS make sure that the behavior videos are in the same directories');
disp(repmat(' ',[1,100]));
disp(' If you want, just to have it printed at the end for your own information, not to get confused when running many animals and sessions together');
A = input(' Enter the Animal number and the session you are going to work on : ' , 's');
Pro = input(' And the processing Stage of Carapa you are going to apply: ' ,'s');
input(' Now press ENTER to choose the corresponding folder and start analysing');
addpath(genpath('C:\Users\usuario\Documents\GitHub'));
carp = carpa;
carp.menu; 
disp(repmat('~',[1,100]));

if ~isempty(A)
    disp([repmat('~',[1,21]) , ['   ',Pro,' was succesfully applied on Animal ' , A , '   '] , repmat('~',[1,21])]);
    disp(repmat('~',[1,100]));
end
disp([repmat('~',[1,45]) , '    End   ' , repmat('~',[1,45])]);
disp(repmat('~',[1,100]));

disp('started at 10:45')
disp('ended:')
datetime