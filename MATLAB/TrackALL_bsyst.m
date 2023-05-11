
function t = tracking(pos, date, posFingerprint, length, path_seg, path_raw, prog_path)
%Convert to the right class types
format_projN = "%s_%s";
posFingerprint = char(posFingerprint)
length = double(length)
path_seg = char(path_seg)
path_raw  = char(path_raw)
prog_path = char(prog_path)

%Add path to the tracX
addpath(genpath(prog_path));
Tracker = TracX.Tracker();

%Set the TracX variables from input
projectName = char(sprintf(format_projN, pos, date));
fileIdentifierFingerprintImages = 'Sub'; % Image idenfier for Brigthfield images;
fileIdentifierWellPositionFingerprint = char(posFingerprint); % Well position identifier if multiwell experiment.
fileIdentifierCellLineage = []; % Image identifier for the Cell Lineage reconstruction (i.e bud neck marker).
imageCropCoordinateArray = []; % Empty if no crop has been applied in CellX, add CellX cropRegionBoundaries coordinates otherwise (from CellX_SCerevisiae_Parameter.xml).
movieLength = length; % Number of timepoints to track
cellsFilePath = char(path_seg); % Path to segmentation results (CellX Style).
imagesFilePath = char(path_raw); % Path to raw images.
cellDivisionType = 'asym'; % Cell division type.

Tracker.createNewTrackingProject(projectName, imagesFilePath, cellsFilePath, fileIdentifierFingerprintImages, fileIdentifierWellPositionFingerprint, fileIdentifierCellLineage, imageCropCoordinateArray, movieLength, cellDivisionType);

Tracker.runTracker();

%Import Tracker parameters
Tracker = TracX.Tracker();
Tracker.createNewTrackingProject(projectName, imagesFilePath, ...
        cellsFilePath, fileIdentifierFingerprintImages, ...
        fileIdentifierWellPositionFingerprint, fileIdentifierCellLineage, ...
        imageCropCoordinateArray, movieLength, cellDivisionType);
Tracker.revertSegmentationImageCrop()

Tracker.configuration.ParameterConfiguration.setMaxExpectedMovement(70.0000);
Tracker.configuration.ParameterConfiguration.setMeanCellDiameter(28.0000);
Tracker.configuration.ParameterConfiguration.setNeighbourhoodSearchRadius(50.0000);

%Run the tracker
Tracker.runTracker();

%Save tracker results
Tracker.saveCurrentTrackerState(); % Saves the tracker state as mat file (to continue work anytime later)
Tracker.saveTrackingProject(); % Saves the tracking project.
Tracker.saveTrackerResultsAsTable(); % Saves the tracking results as one column seperated table for further analysis.
Tracker.saveTrackerProjectControlImages('isParallel', false, 'maxWorkers', 1)
end
