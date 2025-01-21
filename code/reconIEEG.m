function subject = reconIEEG(subject_n)

allsubject = readtable('file_paths.csv', Delimiter=',');

subject = ieeg_recon;
subject.preImplantMRI = allsubject.t1w{subject_n};
subject.postImplantCT = allsubject.ct{subject_n};
subject.postImplantCT_electrodes = allsubject.electrodes{subject_n};

subject.output = ['derivatives/'  allsubject.sub{subject_n}];
subject.fslLoc = '/usr/local/fsl/bin';
subject.itksnap = '/Applications/ITK-SNAP.app/Contents/bin';
subject.freeSurfer = '/Applications/freesurfer/7.3.2/bin';

%% Run Module 1

subject.module1;

%% Run Module 2

ignoreFlag = true;
fileLocations = subject.module2('gc',ignoreFlag);

% % % Quality Assurance Visualization
% imageviewer ='freeview';
% subject.module2_QualityAssurance(fileLocations,imageviewer);

end
