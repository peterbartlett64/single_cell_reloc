function t = tracking(pos, date, posFingerprint, length, path_seg, path_raw, TracX_dir)
%Convert to the right class types
format_projN = "%s_%s";
posFingerprint = char(posFingerprint);
% length = length; %* This may be unnecessary
path_seg = char(path_seg);
path_raw  = char(path_raw);

%Add path to the tracX
addpath(genpath(TracX_dir));
%TrackX_GUI = TracX.GUI.TracXGUI();

%Set the TracX variables from input
projectName = char(sprintf(format_projN, pos, date));
fileIdentifierFingerprintImages = 'Sub'; % Image idenfier for Brigthfield images;
fileIdentifierWellPositionFingerprint = char(posFingerprint); % Well position identifier if multiwell experiment.
fileIdentifierCellLineage = 'Add'; % Image identifier for the Cell Lineage reconstruction (i.e bud neck marker).
imageCropCoordinateArray = []; % Empty if no crop has been applied in CellX, add CellX cropRegionBoundaries coordinates otherwise (from CellX_SCerevisiae_Parameter.xml).
movieLength = length; % Number of timepoints to track
cellsFilePath = path_seg; % Path to segmentation results (CellX Style).
imagesFilePath = path_raw; % Path to raw images.
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
Tracker.configuration.ParameterConfiguration.setDivisionMarkerEdgeSensitivityThresh(0.05) % This is very low but was needed to find the mKa. This could be an issue down the line

%Run the tracker
Tracker.runTracker();

%% Save tracker results
Tracker.saveCurrentTrackerState(); % Saves the tracker state as mat file (to continue work anytime later)
Tracker.saveTrackingProject(); % Saves the tracking project.
Tracker.saveTrackerResultsAsTable(); % Saves the tracking results as one column seperated table for further analysis.
Tracker.saveTrackerProjectControlImages('isParallel', false) %This is parrallel bc it is being run as a test individually

%% Run the lineage calculations
Tracker.runLineageReconstruction('symmetricalDivision', false,'WriteControlImages', true);
%%
% Plot linaege for track 1 to console
trackToDisplay = 3;
Tracker.imageVisualization.plotLineageTree(trackToDisplay);
% Plot linaege for track 1 to figure
fh = Tracker.imageVisualization.plotLineageTree(trackToDisplay, 'plotToFigure', true);
fh{trackToDisplay}.Visible = true;

%% Print the cell cycle phase information as table and save it
Tracker.lineage.cellCyclePhaseTable;
Tracker.saveTrackerCellCycleResultsAsTable()

%% Generate an animated movie of a lineage of interest
Tracker.imageVisualization.plotLineageTree('symmetrical', 2 , 'plotToFigure', false);

%% Save  marker time series control images for lineage roots
% This creates a mask overlay with the budneck marker channel and nuclear
% channel signal
Tracker.imageProcessing.generateMarkerTimeSeriesImage('cycleMarker', ...
    {'mKOk', 'mKate'}, 'isPlotMarkerSignal', true)

end