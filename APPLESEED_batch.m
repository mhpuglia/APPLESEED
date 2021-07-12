% APPLESEED_batch - Loops across multiple EEG datasets to run APPLESEED: The Automated
% 		Preprocessing Pipe-Line for the Estimation of Scale-wise Entropy from EEG Data on
%		each dataset
%
% 		Author   : Meghan H. Puglia (meghan.puglia@virginia.edu) | 2021
%		Download : https://github.com/mhpuglia/APPLESEED
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Usage: 	
%
%       The user must set the following variables under "User-defined variables," below.
%
%		parentdir 	 - 	[string] full path to the parent directory that contains the input 
%					   	datasets
%		chanfile	 - 	[string] full path to the channel location file, required for
%					   	channel interpolation. To learn how to generate a channel file,  
%					   	see https://eeglab.org/tutorials/04_Import/Channel_Locations.html
%		binfile		 - 	[string] full path to the bin file defining event codes in the 
%					  	dataset (see https://github.com/lucklab/erplab/wiki/Assigning-Events-to-Bins-with-BINLISTER:-Tutorial).
%					   	required for task-based data; for resting-state data, this variable 
%					   	should be deleted
%		filetype     -  [string] file extension for EEG data type (e.g. '.edf' (European 
%						Data Format), '.vhdr' (BrainVision Core Data Format header file),
%						'.set' (MATLAB toolbox EEGLAB), '.bdf' (Biosemi), etc.)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% About:
% 
% 		APPLESEED_batch uses APPLESEED_setup() to prepare a dataset for APPLESEED, and 
%		APPLESEED() to preprocess and compute scale-wise entropy on the dataset.
%
%		As written, APPLESEED_setup() will accept any of the four Brain Imaging Data 
%		Structure (BIDS)-approved formats (i.e., '.edf', '.vhdr', '.set','.bdf'). Users may 
%		edit Step 3 of APPLESEED_setup() to import raw EEG data of a different format via 
%		an EEGLAB data inport plugin (see https://eeglab.org/tutorials/04_Import/Importing_Continuous_and_Epoched_Data.html).
%
%		This script runs APPLESEED() with its default parameters. The user may customize
%		preprocessing and scale-wise entropy computation parameters via the optional input
%		arguments. See APPLESEED() for more details.
%
%
% 		An example dataset is available for download from https://openneuro.org/datasets/ds003710.
%		The user should place all APPLESEED scripts within the MATLAB path, download the 
%		example dataset to an "APPLESEED_Example_Dataset" directory, and change into this 
%		directory in MATLAB. The example code below will then run APPLESEED_setup() and 
%		APPLESEED() from the MATLAB command line on the example dataset which includes: EEG 
%		data from 48 recording sessions from 13 infants and a channel location file and a 
%		bin file (located in the "APPLESEED_Example_Dataset > code" directory). 
%
%		See Puglia, Slobin, & Williams (2021) for more information on the example dataset 
%		and the development and optimization of APPLESEED.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Citations:
%
%		In publications, please reference: Puglia, M.H., Slobin, J.S., Williams, C.L., 2021.
%		The Automated Preprocessing Pipe-Line for the Estimation of Scale-wise Entropy from 
%		EEG Data (APPLESEED): Development and validation for use in pediatric populations.
%		bioRxiv. https://doi.org/10.1101/2021.07.10.450198.
%	
% 		Additionally, APPLESEED makes use of the following toolboxes/plugins that should be 
%		installed prior to use and cited in any resulting manuscripts (see also the output 
%		logfile for required citations):
%
% 		EEGLAB: https://sccn.ucsd.edu/eeglab/download.php
%				Citation: Delorme, A., Makeig, S., 2004. EEGLAB: An open source toolbox for 
%							analysis of single-trial EEG dynamics including independent 
%							component analysis. J. Neurosci. Methods 134, 9–21.
%
%		ERPLAB: will be automatically installed as an EEGLAB plugin if necessary
%			    Citation: Lopez-Calderon, J., Luck, S.J., 2014. ERPLAB: an open-source 
%							toolbox for the analysis of event-related potentials. Front. 
%							Hum. Neurosci. 8, 213.
% 
%		MADE:   https://github.com/ChildDevLab/MADE-EEG-preprocessing-pipeline
%			    Citation: Debnath, R., Buzzell, G.A., Morales, S., Bowers, M.E., Leach, S.C., 
%							Fox, N.A., 2020. The Maryland analysis of developmental 
%							EEG (MADE) pipeline. Psychophysiology 57, e13580.
%
%		ADJUST: will be automatically installed as an EEGLAB plugin if necessary
%				Citation: Mognon, A., Jovicich, J., Bruzzone, L., Buiatti, M., 2011. ADJUST: 
%							An automatic EEG artifact detector based on the joint use of 
%							spatial and temporal features. Psychophysiology 48, 229–240. 
%
%		FASTER:	will be automatically installed as an EEGLAB plugin if necessary
%				Citation: Nolan, H., Whelan, R., Reilly, R.B., 2010. FASTER: Fully 
%							Automated Statistical Thresholding for EEG artifact Rejection. 
%							J. Neurosci. Methods 192, 152–162.	
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% License:
%
%    Copyright 2021 Meghan H. Puglia, University of Virginia (meghan.puglia@virignia.edu)
%    
%    This file is part of APPLESEED.
%
%    APPLESEED is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    APPLESEED is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with APPLESEED.  If not, see <https://www.gnu.org/licenses/>.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                           %
%                                           %
%          User-defined variables:          %
%                                           %
%                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parentdir = fullfile(pwd); %Full path to the parent directory that contains the input dataset(s); the default present working directory (pwd) assumes you have changed into the parent directory housing your dataset
chanfile  = fullfile(parentdir,'code','brainvision_32chan_locs.ced'); %Full path to the channel location file
binfile   = fullfile(parentdir,'code','bins_viewing_condition.txt'); %Full path to the bin file, required for task-based analyses; remove this line for resting-state analysis
filetype  = '.vhdr'; %File extension of the eeg data type. If you are not using a BIDS-approved data format (i.e. European Data Format (.edf), BrainVision Core Data Format (.vhdr), MATLAB toolbox EEGLAB (.set), or Biosemi (.bdf)), you must edit Step 3 of APPLESEED_setup() to import your data via the appropriate EEGLAB import plugin.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                             %
%                                                             %
%   EXAMPLE 1: Run APPLESEED for all datasets in parentdir    %
%                                                             %
%                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Find all eeg datasets with the specified file extension within parentdir
	Files=dir(fullfile(parentdir, '**',['*',filetype])); 

%Confirm which files will be run in the batch 
	clc;
	allfiles={Files.name,};
	fprintf('%s\n\n', 'APPLESEED will be run for the following files:');
	fprintf('   %s\n',allfiles{:});

	yn=input('\nDo you wish to proceed (y/n)?: ', 's');
	while ((~ strcmpi(yn, 'y') && ~strcmpi(yn, 'n')))
		yn=input('\nInvalid input. Do you wish to proceed (y/n)?: ', 's');
		disp(yn)
	end

	if strcmpi(yn, 'y') 
		disp('Proceeding...');

	%Looping across all files, run APPLESEED_setup(), then APPLESEED()
		for f = 1:size(Files,1)
	
			%Extract filename base from filename
			[~, filenamebase, ~]=fileparts(Files(f).name);
	
			%Prepare dataset for APPLESEED
			APPLESEED_setup(filenamebase, parentdir, chanfile, 'filetype', filetype); 
	
			if ~exist('binfile', 'var') || isempty(binfile) 
		
				%Run APPLESEED if resting-state
				APPLESEED(filenamebase, parentdir); 
	
			else
	
				%Run APPLESEED if task-based
				APPLESEED(filenamebase, parentdir, binfile); 
	
			end	%if binfile
	
		end %for files	

	else 

		disp('Process cancelled.');
		return;

	end %if y/n	




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                             %
%                                                             %
% EXAMPLE 2: Run APPLESEED for specific datasets in parentdir %
%                                                             %
%                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%List the names of the files in a cell array
Files={'sub-01_ses-1_task-appleseedexample_eeg', 'sub-01_ses-2_task-appleseedexample_eeg', 'sub-01_ses-3_task-appleseedexample_eeg', 'sub-01_ses-4_task-appleseedexample_eeg'}; 

for f = 1:size(Files,2)
	
	%Set filenamebase
	filenamebase=Files{f};   
	
	%Prepare dataset for APPLESEED
	APPLESEED_setup(filenamebase, parentdir, chanfile); 
	
	if ~exist('binfile', 'var') || isempty(binfile) 
		
		%Run APPLESEED if resting-state
		APPLESEED(filenamebase, parentdir); 
	
	else
		
		%Run APPLESEED if task-based
		APPLESEED(filenamebase, parentdir, binfile); 
	
	end	%if binfile
	
end %for files
