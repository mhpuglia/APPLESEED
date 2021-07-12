% APPLESEED_setup() - Prepares an EEG dataset for use in APPLESEED: The Automated
% 		Preprocessing Pipe-Line for the Estimation of Scale-wise Entropy from EEG Data 
%
% 		Author   : Meghan H. Puglia (meghan.puglia@virginia.edu) | 2021
%		Download : https://github.com/mhpuglia/APPLESEED
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Usage: 	
%
%	APPLESEED_setup(filenamebase, parentdir, chanfile); 
%
%
% 	Inputs:
%		filenamebase - 	[string] name of the input dataset without the file extension
%		parentdir 	 - 	[string] full path to the parent directory that contains the input 
%						dataset
%		chanfile	 - 	[string] full path to the channel location file, required for
%						channel interpolation; to learn how to generate a channel file,  
%						see https://eeglab.org/tutorials/04_Import/Channel_Locations.html
%
% 	Optional inputs:
%		'filetype'   -  [string] file extension of the EEG data {default: will attempt to 
%						load a BIDS-approved file type i.e., '.edf' (European Data Format), 
%						'.vhdr' (BrainVision Core Data Format header file), '.set' (MATLAB 
%						toolbox EEGLAB), or '.bdf' (Biosemi)}
%		'biosmeiref' -  [string] name of the channel to use as a reference for BIOSEMI type 
%						data. {default: channel closest to the central zenith}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% About:
%
%		This script prepares EEG data for use in APPLESEED, which involves loading the raw
%		data, adding a channel location file, and saving the dataset in EEGLAB format. 
%		
%		As written, this script will load EEG data in any of the four Brain Imaging Data 
%		Structure (BIDS)-approved formats (i.e., '.edf', '.vhdr', '.set','.bdf'). Users may 
%		edit Step 3 of this script to import raw EEG data of a different format via 
%		an EEGLAB data inport plugin (see https://eeglab.org/tutorials/04_Import/Importing_Continuous_and_Epoched_Data.html).
%
%		See Puglia, Slobin, & Williams (2021) for more information on the development and
%		optimization of APPLESEED. 
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
%		installed prior to use and cited in any resulting manuscripts (see also output 
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



function outdir = APPLESEED_setup(filenamebase, parentdir, chanfile, varargin)
	
	p = [];
	if nargin < 3
	
		help APPLESEED_setup; return;
	
	elseif nargin > 3
		
		%Read optional input parameters
 		for index = 1:2:length(varargin)
			p = setfield(p, varargin{index}, varargin{index+1});
 		end
	
	end %if nargin
	
		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   %
% Step 1: Determine input file type %
%                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		

	%Determine file type if not provided
	if ~isfield(p, 'filetype')
		locs=dir(fullfile(parentdir, '**', [filenamebase, '*']));
		if ~isempty(locs)
			[~, ~, p.filetype]=cellfun(@fileparts, {locs.name}, 'UniformOutput', false);
		else
			error(['Oh no! APPLESEED_setup has FAILED! A file named ''', [filenamebase],''' cannot be found within ',parentdir]);
		end
		
		if sum(strcmpi(p.filetype, '.vhdr')) > 0 && sum(strcmpi(p.filetype, '.eeg')) > 0; p.filetype='.vhdr'; end
		if sum(strcmpi(p.filetype, '.set')) > 0; p.filetype= '.set'; end
		if sum(strcmpi(p.filetype, '.bdf')) > 0; p.filetype = '.bdf'; end
		if sum(strcmpi(p.filetype, '.edf')) > 0; p.filetype = '.edf'; end		
	
	end
	
		


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                      %
% Step 2: Check that input files exist %
%                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%Dataset
	datasetdir=fullfile(parentdir,filenamebase);

	if exist(fullfile(datasetdir, [filenamebase, p.filetype]), 'file')==0 
		if exist(fullfile(parentdir, [filenamebase, p.filetype]), 'file')==0  
			try
				locs=dir(fullfile(parentdir, '**', [filenamebase, p.filetype]));
				if size(locs,2)==1
					inputdir=locs.folder;
				end %if subdir contains dataset
			catch
				error(['Oh no! APPLESEED_setup has FAILED! The EEG dataset ''', [filenamebase, p.filetype],''' cannot be found within ',parentdir]);
			end %try dir
		else
			inputdir=parentdir;
		end %if parentdir contains dataset
	else
		inputdir=datasetdir;
	end %if datasetdir contains dataset


	%Channel location file
	if exist(chanfile, 'file')==0
		error(['Oh no! APPLESEED_setup has FAILED! The channel location file does not exist: ', chanfile]);
	end

	%Initialize EEGLAB 
	eeglab;	


	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            %
% Step 3: Import EEG dataset %
%                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% NOTE!: If your EEG data is not in a BIDS-approved format, edit this section to import your data using an alternate EEGLAB data inport plugin (see https://eeglab.org/tutorials/04_Import/Importing_Continuous_and_Epoched_Data.html)
	
	switch p.filetype

		case '.vhdr'
			
			%Make sure bva-io plugin is installed 
			Test=which('pop_loadbv');
			if isempty(Test)
				warning('The bva-io plugin is required to import Brain Vision Analyzer EEG data files. APPLESEED_setup is installing this plugin now.');
				plugin_askinstall('bva-io', [], true);
			end

			EEG = pop_loadbv(inputdir, [filenamebase,p.filetype], [], []);
		
		case '.set'
			
			EEG = pop_loadset('filename',[filenamebase,p.filetype],'filepath',fullfile(inputdir));
		
		case '.bdf'
		
			%Set reference channel for import (required for BIOSEMI)
			if isfield(p, 'biosemiref')
				refindex=find(strcmpi(p.biosemiref, {A.labels}));
			else
				A=readlocs(chanfile);
				[~,refindex]=sortrows(horzcat([A.radius]',abs([A.theta])'), [1,2]);
				refindex=refindex(1);
				warning(['You did not specify a channel for importing BIOSEMI data. Channel', A(refindex).label,' was selected for you.']);
			end
			
			EEG = pop_biosig(fullfile(inputdir, [filenamebase,p.filetype]), 'ref',refindex);
			
		case '.edf'
		
			EEG = pop_biosig(fullfile(inputdir, [filenamebase,p.filetype]));
		
		otherwise
		
			error('Oh no! APPLESEED_setup has FAILED! Your data is not in a BIDS-approved format. You must edit Step 3 of this script to import your data type.');

	end %switch filetype


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   %
% Step 4: Add channel location file %
%                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	EEG = pop_chanedit(EEG, 'lookup',chanfile);
	

%%%%%%%%%%%%%%%%%%%%%%%%
%                      %
% Step 5: Save dataset %
%                      %
%%%%%%%%%%%%%%%%%%%%%%%%

	%Create new 'derivatives' directory in parentdir with sub/[ses] level directory
	outdir=	fullfile(parentdir, 'derivatives', strrep(inputdir, parentdir, ''));
	mkdir(outdir);
	EEG = pop_saveset(EEG, 'filename', [filenamebase,'.set'], 'filepath', outdir);

		