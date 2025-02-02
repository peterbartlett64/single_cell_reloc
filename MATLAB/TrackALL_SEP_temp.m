%%
%function t = TrackALL(pos, date, posFingerprint, length, path_seg, path_raw, par, quant)
%Convert to the right class types
format_projN = "%s_%s";
%posFingerprint = char(posFingerprint);
%length = double(length);
%path_seg = char(path_seg);
%path_raw  = char(path_raw);

date = 'Sep13';
pos = 'd0214r2p090300'
posFingerprint = char('position090300');
length = double(97);
path_seg = char('E:\Microfluidics\RESULTS\d0214r2p090300\GFP_mKO_mKa');
path_raw  = char('E:\Microfluidics\Analyze\20180214_hr2000_MMS\2018-02-14_20-25-52');
%projectName = char('testing_stuff');
par = false;
quant = false;

%Add path to the tracX
addpath(genpath("C:\Users\pcnba\OneDrive\Desktop\Test\tracx\src\"));
%addpath(genpath("C:\Users\Nikon\Desktop\tracx\src\"));

%Set the TracX variables from input
projectName = char(sprintf(format_projN, pos, date)); %char(sprintf(format_projN, pos, date));
fileIdentifierFingerprintImages = 'Sub'; % Image idenfier for Brigthfield images;
fileIdentifierWellPositionFingerprint = char(posFingerprint); % Well position identifier if multiwell experiment.
fileIdentifierCellLineage = []; % Image identifier for the Cell Lineage reconstruction (i.e bud neck marker).
imageCropCoordinateArray = []; % Empty if no crop has been applied in CellX, add CellX cropRegionBoundaries coordinates otherwise (from CellX_SCerevisiae_Parameter.xml).
movieLength = length; % Number of timepoints to track
cellsFilePath = char(path_seg); % Path to segmentation results (CellX Style).
imagesFilePath = char(path_raw); % Path to raw images.
cellDivisionType = 'asym'; % Cell division type.


%Import Tracker parameters
Tracker = TracX.Tracker();
Tracker.createNewTrackingProject(projectName, imagesFilePath, ...
        cellsFilePath, fileIdentifierFingerprintImages, ...
        fileIdentifierWellPositionFingerprint, fileIdentifierCellLineage, ...
        imageCropCoordinateArray, movieLength, cellDivisionType);

Tracker.revertSegmentationImageCrop();
Tracker.configuration.ParameterConfiguration.setMaxExpectedMovement(70.0000);
Tracker.configuration.ParameterConfiguration.setMeanCellDiameter(28.0000);
Tracker.configuration.ParameterConfiguration.setNeighbourhoodSearchRadius(50.0000);

%Run the tracker
Tracker.runTracker();
%%Save tracker results
Tracker.saveCurrentTrackerState(); % Saves the tracker state as mat file (to continue work anytime later)
Tracker.saveTrackingProject(); % Saves the tracking project.

% Saves the tracking results as one column seperated table for further analysis. This is dependent on if the Quantification was captured for
% segmentation
%%
Tracker.saveTrackerResultsAsTable();
Tracker.saveTrackerProjectControlImages('isParallel', false, 'maxWorkers', 1);

%if quant
 %   Tracker.saveTrackerResultsAsTable();
%else
 %   Tracker.saveTrackerResultsAsTable_nQ();
%end 

%Save the control images
%if par
 %   Tracker.saveTrackerProjectControlImages('isParallel', false, 'maxWorkers', 1);
%else
 %   Tracker.saveTrackerProjectControlImages('isParallel', true, 'maxWorkers', 4); % This is the number of logical 
%end

%Tracker.data.QuantificationData = []
%rmfield(Tracker.data, QuantificationData)

%end
