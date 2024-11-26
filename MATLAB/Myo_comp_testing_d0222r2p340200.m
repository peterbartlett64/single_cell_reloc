%function t = tracking(pos, date, posFingerprint, length, path_seg, path_raw)
%Convert to the right class types
%format_projN = "%s_%s";
posFingerprint = char('position340200');
length = 69;
path_seg = char('G:\testingSUPPORT\New_training\Just the one\GFP_mKO_mKa');
path_raw  = char('G:\testingSUPPORT\New_training\Just the one\d0222r2p340200');
%%
%Add path to the tracX
addpath(genpath("C:\Users\Nikon\Desktop\tracx\src"));
%TrackX_GUI = TracX.GUI.TracXGUI();
%%

%Set the TracX variables from input
projectName = char('test_composite_d0222r2p340200');
fileIdentifierFingerprintImages = 'Sub'; % Image idenfier for Brigthfield images;
fileIdentifierWellPositionFingerprint = char(posFingerprint); % Well position identifier if multiwell experiment.
fileIdentifierCellLineage = 'Add'; % Image identifier for the Cell Lineage reconstruction (i.e bud neck marker).
imageCropCoordinateArray = []; % Empty if no crop has been applied in CellX, add CellX cropRegionBoundaries coordinates otherwise (from CellX_SCerevisiae_Parameter.xml).
movieLength = length; % Number of timepoints to track
cellsFilePath = path_seg; % Path to segmentation results (CellX Style).
imagesFilePath = path_raw; % Path to raw images.
cellDivisionType = 'asym'; % Cell division type.

%%

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
%% Show a few plots with the default values before tuning
Tracker.testBudNeckSegmentation(60);
%%
%Tracker.testBudNeckSegmentation(62);
%%
%Tracker.testBudNeckSegmentation(64);
%%
%Tracker.testBudNeckSegmentation(66);
%%
%Test first with the default params
%Tracker.configuration.ParameterConfiguration.setDivisionMarkerMaxExpectedMovement(20);
%Tracker.configuration.ParameterConfiguration.setDivisionMarkerMaxObjCenterDisplacement(20);
%Tracker.configuration.ParameterConfiguration.setDivisionMarkerMaxObjSizeDecrease(0.7);
%Tracker.configuration.ParameterConfiguration.setDivisionMarkerMaxObjSizeIncrease(4);
%Tracker.configuration.ParameterConfiguration.setDivisionMarkerMaxTrackFrameSkipping(9);
%Tracker.configuration.ParameterConfiguration.setDivisionMarkerEdgeSensitivityThresh(0.1)

%% Here are all the division marker variables. Starting with the default values
Tracker.configuration.ParameterConfiguration.setDivisionMarkerEdgeSensitivityThresh(0.30); % This will likely need to be changed to something like 0.15
Tracker.configuration.ParameterConfiguration.setDivisionMarkerDenoising(25);
Tracker.configuration.ParameterConfiguration.setDivisionMarkerConvexAreaUpperThresh(100);
Tracker.configuration.ParameterConfiguration.setDivisionMarkerConvexAreaLowerThresh(5);
Tracker.configuration.ParameterConfiguration.setDivisionMarkerMaxObjCenterDisplacement(25);
Tracker.configuration.ParameterConfiguration.setDivisionMarkerMaxExpectedMovement(15);
Tracker.configuration.ParameterConfiguration.setDivisionMarkerAverageObjSizeGrowth(0.0300);
Tracker.configuration.ParameterConfiguration.setDivisionMarkerMaxObjSizeDecrease(0.75);
Tracker.configuration.ParameterConfiguration.setDivisionMarkerMaxObjSizeIncrease(2.5000);
Tracker.configuration.ParameterConfiguration.setDivisionMarkerMaxMajorAxisRotation(0); % This could be causing issue
Tracker.configuration.ParameterConfiguration.setDivisionMarkerIndividualFunctionPenalty(500);
Tracker.configuration.ParameterConfiguration.setDivisionMarkerMaxTrackFrameSkipping(2); %This might need to be 9
Tracker.configuration.ParameterConfiguration.setDivisionMarkerUsedFunctionsForCostMatrix([1101]);
Tracker.configuration.ParameterConfiguration.setDivisionMarkerMeanObjDiameter(30);
Tracker.configuration.ParameterConfiguration.setDivisionMarkerMeanObjDiameterScalingFactor(2.5000);
Tracker.configuration.ParameterConfiguration.setDivisionMarkerFingerprintHalfWindowSideLength(25);
Tracker.configuration.ParameterConfiguration.setDivisionMarkerFingerprintResizeFactor(32);
Tracker.configuration.ParameterConfiguration.setDivisionMarkerFingerprintMaxConsideredFrequencies(8);
Tracker.configuration.ParameterConfiguration.setDivisionMarkerNeighbourhoodSearchRadius(40);
Tracker.configuration.ParameterConfiguration.setDivisionMarkerRadiusVectorFilterX(100);
Tracker.configuration.ParameterConfiguration.setDivisionMarkerRadiusVectorFilterY(100);
Tracker.configuration.ParameterConfiguration.setDivisionMarkerRadiusVectorFilterZ(100);
Tracker.configuration.ParameterConfiguration.setDivisionMarkerLineProfileLength(100);
Tracker.configuration.ParameterConfiguration.setDivisionMarkerLineProfilePeakLowerBound(10);
Tracker.configuration.ParameterConfiguration.setDivisionMarkerLineProfilePeakUpperBound(90);

Tracker.testMarkerTrackingParameters([30, 60], 'divisionMarkerMaxTrackFrameSkipping', 9);
%%

%%
Tracker.configuration.ParameterConfiguration.setDivisionMarkerMeanObjDiameter(9)
Tracker.testBudNeckSegmentation(27, 'divisionMarkerEdgeSensitivityThresh', 0.2)
%%
Tracker.configuration.ParameterConfiguration.setDivisionMarkerMaxExpectedMovement(10)
Tracker.configuration.ParameterConfiguration.setDivisionMarkerMaxObjCenterDisplacement(10)
Tracker.configuration.ParameterConfiguration.setDivisionMarkerMaxObjSizeDecrease(1);
Tracker.configuration.ParameterConfiguration.setDivisionMarkerMaxObjSizeIncrease(5);
%%
Tracker.testMarkerTrackingParameters([30, 60], 'divisionMarkerMaxTrackFrameSkipping', 4)
%%
%Tracker.configuration.ParameterConfiguration.setDivisionMarkerEdgeSensitivityThresh(0.2);% This is very low but was needed to find the mKa. This could be an issue down the line

%All of the below are dependent on multiple frames. Should test first with
%the test keys

Tracker.configuration.ParameterConfiguration.setDivisionMarkerMeanObjDiameterScalingFactor(1);
%Tracker.configuration.ParameterConfiguration.setDivisionMarkerMaxTrackFrameSkipping(1)
%Tracker.configuration.ParameterConfiguration.setDivisionMarkerEdgeSensitivityThresh(0.18)
Tracker.configuration.ParameterConfiguration.setDivisionMarkerMaxObjSizeDecrease(0.5);

Tracker.testMarkerTrackingParameters([1, 30], 'divisionMarkerMaxTrackFrameSkipping', 1)
%%
%Tracker.testBudNeckSegmentation(42, 'divisionMarkerEdgeSensitivityThresh', 0.2)
%% Save tracker results
Tracker.saveCurrentTrackerState(); % Saves the tracker state as mat file (to continue work anytime later)
Tracker.saveTrackingProject(); % Saves the tracking project.
Tracker.saveTrackerResultsAsTable(); % Saves the tracking results as one column seperated table for further analysis.
Tracker.saveTrackerProjectControlImages('isParallel', true) %This is parrallel bc it is being run as a test individually

%% Run the lineage calculations
Tracker.runLineageReconstruction('symmetricalDivision', false,'WriteControlImages', false); %At this time the control images cannot be outputed
%%
% Plot linaege for track 1 to console
trackToDisplay = 9;
Tracker.imageVisualization.plotLineageTree(trackToDisplay);
% Plot linaege for track 1 to figure
fh = Tracker.imageVisualization.plotLineageTree(trackToDisplay, 'plotToFigure', true);
fh{trackToDisplay}.Visible = true;

%% Print the cell cycle phase information as table and save it
Tracker.lineage.cellCyclePhaseTable;
Tracker.saveTrackerCellCycleResultsAsTable()
%%
Tracker.imageProcessing.generateLineageMovie(2, 5/60, 'Track21LineageMovie')
%% Generate an animated movie of a lineage of interest
Tracker.imageVisualization.plotLineageTree('symmetrical', 2 , 'plotToFigure', false);

%% Save  marker time series control images for lineage roots
% This creates a mask overlay with the budneck marker channel and nuclear
% channel signal
Tracker.imageProcessing.generateMarkerTimeSeriesImage('cycleMarker', ...
    {'Add'}, 'isPlotMarkerSignal', true)
