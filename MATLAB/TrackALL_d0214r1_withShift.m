
%function t = tracking(pos, date, posFingerprint, length, path_seg, path_raw)
%Convert to the right class types
%format_projN = "%s_%s";
posFingerprint = char('position090200');
length = 83;
path_seg = char("C:\Users\pcnba\Grant Brown's Lab Dropbox\Peter Bartlett\Peter Bartlett Data\Testing_shift_fix\d0214r1p090200\GFP_mKO_mKa");
path_raw  = char("C:\Users\pcnba\Grant Brown's Lab Dropbox\Peter Bartlett\Peter Bartlett Data\Testing_shift_fix\Images");


%Add path to the tracX

addpath(genpath("C:\Users\pcnba\Grant Brown's Lab Dropbox\Peter Bartlett\Peter Bartlett Data\Code\single_cell_reloc\tracx\src\"));

%Set the TracX variables from input
projectName = char("testing_w_big_shift");
fileIdentifierFingerprintImages = 'Sub'; % Image idenfier for Brigthfield images;
fileIdentifierWellPositionFingerprint = char(posFingerprint); % Well position identifier if multiwell experiment.
fileIdentifierCellLineage = []; % Image identifier for the Cell Lineage reconstruction (i.e bud neck marker).
imageCropCoordinateArray = []; % Empty if no crop has been applied in CellX, add CellX cropRegionBoundaries coordinates otherwise (from CellX_SCerevisiae_Parameter.xml).
movieLength = length; % Number of timepoints to track
cellsFilePath = char(path_seg); % Path to segmentation results (CellX Style).
imagesFilePath = char(path_raw); % Path to raw images.3
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
Tracker.saveTrackerResultsAsTable(); % Saves the tracking results as one column seperated table for further analysis.
Tracker.saveTrackerProjectControlImages('isParallel', false, 'maxWorkers', 1)

%end
