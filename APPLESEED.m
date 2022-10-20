% APPLESEED() - The Automated Preprocessing Pipe-Line for the Estimation of Scale-wise Entropy from EEG Data 
%
%		Version  : 1.0
% 		Author   : Meghan H. Puglia (meghan.puglia@virginia.edu) | 2021
%		Download : https://github.com/mhpuglia/APPLESEED
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Usage:
% 	
%		APPLESEED(filenamebase, parentdir, binfile); %for task-based data
%		APPLESEED(filenamebase, parentdir);	%for resting-state data
% 
%
% 	Inputs:
%		filenamebase - 	[string] name of the input dataset without the file extension
%		parentdir 	 - 	[string] full path to the parent directory that contains the input 
%					   	dataset
%		binfile		 - 	[string] required for task-based analyses; full path to the bin 
%					   	file defining event codes in the dataset (see 
%					   	https://github.com/lucklab/erplab/wiki/Assigning-Events-to-Bins-with-BINLISTER:-Tutorial)
%
%
% 	Optional inputs:
%       'outputdir'     -   [string] specify an alternate output directory
%		'saverobust'	-	['on'|'off'] save interim datasets enabling the inspection of 
%						  	artifacts/components marked for rejection {default: 'on'}
%		'resamp'		- 	[in Hz] resampling rate {default: 250}
%		'hp'			- 	[in Hz] high-pass filter cutoff  {default: 0.3}
%		'lp'			- 	[in Hz] low-pass filter cutoff {default: 50}
%		'eplen'			- 	[in ms] epoch length {default: 1000}
%		'arxt'			- 	[in µV] threshold for identification of extreme voltage 
%						  	artifacts {default: 500}
%		'runica'		- 	['on'|'off'] decompose data into independent components for the 
%						  	identification of data artifacts {default: 'on'}
%		'icarejmethd' 	- 	['adjusted_ADJUST'|'ADJUST'|'manual'] the method to identify 
%						  	components contaminated with artifacts {default: 'adjusted_ADJUST'}
%		'arsd'			- 	[in µV] threshold for identification of artifacts that exceed 
%						  	this maximum standard deviation change within a 200 ms moving 
%						  	window {default: 80}
%       'arauto'        -   ['on'|'off'] automatically detect and reject artifactual epochs. If 
%                           'off' is selected, users will manually reject artifacts, and must 
%                           create a text file which contains the index of each epoch they wish 
%                           to remove, and supply the full path to this text file with the 'arfile' 
%                           input argument. If the 'arfile' input argument is specified, this parameter
%                           will automatically be set to 'off' {default: 'on'}
%       'arfile'        -   A string specifying the full path to the file containing manually identified 
%                           artifactual epochs. Epochs marked as '0' will be retained, epochs marked as '1' 
%                           will be removed. The file must be a plain text file with the number of rows 
%                           corresponding to the number of epochs. {defualt: This argument is only required 
%                           when 'arauto' is set to 'off' 
%		'chaninterp'	- 	['on'|'off'] interpolate channels identified as problematic via 
%						  	the FASTER algorithm {default: 'on'}
%		'chanfile'		- 	[string] full path to the channel location file, required for
%						  	channel interpolation; to learn how to generate a channel file,  
%						  	see https://eeglab.org/tutorials/04_Import/Channel_Locations.html
%						  	{default: this file should be saved to the input dataset following 
%						  	data import}
%		'fastref'		- 	[string] name of the channel to serve as the reference for the 
%						  	FASTER algorithm {default: 'Cz'}
%		'reref'			- 	[string] name of the channel(s) for referencing; multiple channels
%							should be space separated {default: 'Average'}
%		'trselcnt'		- 	[integer] number of trials to retain for scale-wise entropy 
%						  	computation; this number should be identical for all datasets in 
%						  	a study {default: 10}
%		'trselmethd'	- 	['gfp'|'first'|'middle'|'last'] method to identify which trials 
%						  	should be retained for scale-wise entropy computation 
%						  	{default: 'gfp' i.e., global field power}
%		'm'				- 	[integer] pattern length for scale-wise entropy computation 
%						  	{default: 2}
%		'r'				- 	[numeric] similarity criterion for scale-wise entropy computation 
%						  	{default: .5}		
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% About:
%
% 		APPLESEED (Puglia, Slobin, & Williams, 2022) preprocesses task-based or resting-state 
%		EEG data and then computes scale-wise entropy on the data. Scale-wise entropy is an 
%		adaptation of multiscale entropy (Costa et al., 2002) in which sample entropy 
%		(Richman & Moorman, 2000) is calculated across discontinuous data epochs (Grandy 
%		et al., 2016) on the residuals (i.e., after subtracting the within-person average 
%		response across trials/epochs) of the preprocessed EEG data, with the entropy 
%		similarity criterion parameter, r, recalculated at each scale. 
%		
%		Costa, M., Goldberger, A.L., Peng, C.-K., 2002. Multiscale Entropy Analysis of 
%			Complex Physiologic Time Series. Phys. Rev. Lett. 89, 068102. 
%		Grandy, T.H., Garrett, D.D., Schmiedek, F., Werkle-Bergner, M., 2016. On the 
%			estimation of brain signal entropy from sparse neuroimaging data. Sci. Rep. 
%			6, 23073. 
%		Puglia, M.H., Slobin, J.S., Williams, C.L., 2022. The Automated Preprocessing 
%			Pipe-Line for the Estimation of Scale-wise Entropy from EEG Data (APPLESEED): 
%			Development and validation for use in pediatric populations. Developmental Cognitive 
%           Neuroscience, 101163.
%		Richman, J.S., Moorman, J.R., 2000. Physiological time-series analysis using 
%			approximate entropy and sample entropy. Am. J. Physiol. Circ. Physiol. 278, 
%			H2039ÐH2049. 
%
%
% 		An example dataset is available for download from https://openneuro.org/datasets/ds003710.
%		The user should place all APPLESEED scripts within the MATLAB path, download the 
%		example dataset to an "APPLESEED_Example_Dataset" directory, and change into this 
%		directory in MATLAB. The example code in APPLESEED_batch will then run  APPLESEED() from 
%       the MATLAB command line on the example dataset which includes:  EEG data from 48 recording 
%       sessions from 13 infants and a channel location file and a bin file (located in the 
%       "APPLESEED_Example_Dataset > code" directory). 
%
%		See Puglia, Slobin, & Williams (2022) for more information on the example dataset 
%		and the development and optimization of APPLESEED.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Citations:
%
%		In publications, please reference: Puglia, M.H., Slobin, J.S., Williams, C.L., 2022.
%		The Automated Preprocessing Pipe-Line for the Estimation of Scale-wise Entropy from 
%		EEG Data (APPLESEED): Development and validation for use in pediatric populations.
%		Developmental Cognitive Neuroscience, 101163.
%
% 		APPLESEED makes use of the following toolboxes/plugins that should be installed prior 
%		to use and cited in any resulting manuscripts (see also output logfile for citations):
%
% 		EEGLAB: https://sccn.ucsd.edu/eeglab/download.php
%				Citation: Delorme, A., Makeig, S., 2004. EEGLAB: An open source toolbox for 
%							analysis of single-trial EEG dynamics including independent 
%							component analysis. J. Neurosci. Methods 134, 9Ð21.
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
%							spatial and temporal features. Psychophysiology 48, 229Ð240. 
%
%		FASTER:	will be automatically installed as an EEGLAB plugin if necessary
%				Citation: Nolan, H., Whelan, R., Reilly, R.B., 2010. FASTER: Fully 
%							Automated Statistical Thresholding for EEG artifact Rejection. 
%							J. Neurosci. Methods 192, 152Ð162.	
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


function APPLESEED_new20221014(filenamebase, parentdir, varargin)

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            %
% Set parameters/directories %
%                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Display help file if not enough inputs
	if nargin < 2
		help APPLESEED; return;%!quit;
	end


	%If an odd number of inputs, take 3rd input to be binfile
	if rem(length(varargin),2) == 1 
		p.binfile=varargin{1};
		varargin=varargin(2:size(varargin,2));
	end

	
	%Read optional input parameters
 	if ~exist('p', 'var'); 	p = []; end
	for index = 1:2:length(varargin)
 	   p.(varargin{index}) = varargin{index+1}; %setfield(p, varargin{index}, varargin{index+1});
	end

	%If optional parameters not specified, set to default
    if ~isfield(p,'saverobust');	p.saverobust	= 'on'; 			 end
	if ~isfield(p,'resamp');		p.resamp		= 250; 				 else; if ~isnumeric(p.resamp); p.resamp=str2double(p.resamp); end; end
	if ~isfield(p,'hp');			p.hp			= .3; 				 else; if ~isnumeric(p.hp);  p.hp=str2double(p.hp); end; end
	if ~isfield(p,'lp');			p.lp			= 50; 				 else; if ~isnumeric(p.lp); p.lp=str2double(p.lp); end;end
	if ~isfield(p,'eplen');			p.eplen			= 1000; 			 else; if ~isnumeric(p.eplen); p.eplen=str2double(p.eplen); end;end
	if ~isfield(p,'arxt');			p.arxt			= 500; 				 else; if ~isnumeric(p.arxt); p.arxt=str2double(p.arxt); end;end
	if ~isfield(p,'runica');		p.runica		= 'on'; 			 end
	if ~isfield(p,'icarejmethd');	p.icarejmethd	= 'adjusted_ADJUST'; end 
	if ~isfield(p,'arsd');			p.arsd			= 80; 				 else; if ~isnumeric(p.arsd); p.arsd=str2double(p.arsd); end;end
	if ~isfield(p,'arauto');		p.arauto		= 'on';              end
    if ~isfield(p,'chaninterp');	p.chaninterp 	= 'on';				 end 
	if ~isfield(p,'fastref');		p.fastref		= 'Cz'; 			 end
	if ~isfield(p,'reref');			p.reref			= 'Average';		 else; p.reref=strsplit(p.reref, ' '); end 
	if ~isfield(p,'trselcnt');		p.trselcnt		= 10; 				 else; if ~isnumeric(p.trselcnt); p.trselcnt=str2double(p.trselcnt); end;end
	if ~isfield(p,'trselmethd');	p.trselmethd	= 'gfp'; 			 end
	if ~isfield(p,'m');				p.m				= 2; 				 else; if ~isnumeric(p.m); p.m=str2double(p.m); end;end
	if ~isfield(p,'r');				p.r				= .5; 				 else;if ~isnumeric(p.r); p.r=str2double(p.r); end;end
    if isfield(p,'arfile');         p.arauto        = 'off';             end
    
	%Set filename & path
	if isnumeric(filenamebase); filenamebase=num2str(filenamebase); else; filenamebase=strrep(filenamebase,'.set',''); end
	datasetdir=fullfile(parentdir,'derivatives',filenamebase);

	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Check that all specified files exist %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Input dataset
	if exist(fullfile(datasetdir, [filenamebase, '.set']), 'file')==0 || exist(fullfile(datasetdir, [filenamebase, '.fdt']), 'file')==0 
		if exist(fullfile(parentdir, [filenamebase, '.set']), 'file')==0 || exist(fullfile(parentdir, [filenamebase, '.fdt']), 'file')==0  
			try
				locs=dir(fullfile(parentdir, '**', [filenamebase, '.set']));
				if size(locs,2)==1
					inputdir=locs.folder;
				end %if subdir contains dataset
			catch
				error(['Oh no! APPLESEED has FAILED! The EEG dataset ''', [filenamebase, '.set'],''' cannot be found within ',parentdir]);
			end %try dur
		else
			inputdir=parentdir;
		end %if parentdir contains dataset
	else
		inputdir=datasetdir;
	end %if datasetdir contains dataset
	
	
	%Outdir
	if ~isfield(p,'outputdir')     
        outdir = strrep(fullfile(parentdir, 'appleseed', strrep(inputdir, parentdir, '')), [filesep, 'derivatives'], ''); 
	else
        outdir = fullfile(p.outputdir, strrep(inputdir, parentdir, '')); %[EDIT 20221014]
	end
	if ~exist(outdir, 'dir'); mkdir(outdir); end
	
	
	%Binfile
	if isfield(p,'binfile') && exist(p.binfile, 'file')==0
		error(['Oh no! APPLESEED has FAILED! The bin file does not exist: ', p.binfile]);
	end 

    %ARfile [EDIT 20221014]
	if isfield(p,'arfile') && exist(p.arfile, 'file')==0
		error(['Oh no! APPLESEED has FAILED! The artifact rejection file does not exist: ', p.arfile]);
	end 
	
	%If binfile specified, extract binfile name from fullpath
	if isfield(p, 'binfile'); [~,p.binfilename,~]=fileparts(p.binfile); end


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Check that required packages are installed & in path %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%EEGLAB
	try 
		eeglab;
	catch
		error('Oh no! APPLESEED has FAILED! APPLESEED requires the EEGLAB toolbox, but it is not in your MATLAB path.\n%s', 'Please download EEGLAB and add it to your path: https://sccn.ucsd.edu/eeglab/download.php'); 
	end

	
	%ERPLAB
	Test=which('pop_basicfilter');
	if isempty(Test)
		warning('The ERPLAP plugin is required for APPLESEED. APPLESEED is installing this plugin now.')
		plugin_askinstall('ERPLAB', [], true);
	end
	clear Test
	
	
	%MADE (if using adjusted_ADJUST)
	if strcmp(p.runica,'on') == 1 && strcmp(p.icarejmethd,'adjusted_ADJUST')== 1
		Test=which('ajusted_ADJUST');
		if isempty(Test)
			P=what('MADE-EEG-preprocessing-pipeline-master');
			if size(P,1)==1
				addpath(genpath(P.path));
			elseif size(P,1)==0
				error('Oh no! APPLESEED has FAILED! You have selected adjusted_ADJUST to remove components, but the MADE-EEG-preprocessing-pipeline is not in your MATLAB path.\n%s','Please download MADE and add it to your path: https://github.com/ChildDevLab/MADE-EEG-preprocessing-pipeline.');
			else
			   addpath(genpath(P(end).path));
			end %if no MADE
			
			%Edit MARA script to fix error w/ spectopo - [EDIT ADDED 20220321]
			[filepath,~,~]=fileparts(which('MARA_extract_time_freq_features'));
			if ~exist(fullfile(filepath, 'appleseed_edit'), 'dir') || ~contains(filepath, 'appleseed_edit')
				A = regexp(fileread(which('MARA_extract_time_freq_features')),'\n','split');
				whichline = find(contains(A,'[pxx,freq,speccomp,contrib,specstd] = spectopo(icacomps(ic,:), fs, fs);'));

				if ~isempty(whichline)
					A{whichline}=strrep(A{whichline}, 'fs, fs', '0, fs');
					mkdir(fullfile(filepath, 'appleseed_edit'));
					fid=fopen(fullfile(filepath, 'appleseed_edit', 'MARA_extract_time_freq_features.m'), 'w');
					for i = 1:numel(A)
						fprintf(fid,'%s\n', A{i});
					end
					fclose(fid);
					addpath(fullfile(filepath, 'appleseed_edit'));
				end
			else
				addpath(fullfile(filepath, 'appleseed_edit'));
			end %if edited MARA does not exist		
			
		end %if no adjusted_ADJUST
	end %if runica w/ adjusted_ADJUST
	clear Test

	
	%ADJUST (if using ADJUST or adjusted_ADJUST)
	if strcmp(p.runica,'on') == 1 && strcmp(p.icarejmethd,'adjusted_ADJUST')== 1 || strcmp(p.runica,'on') == 1 && strcmp(p.icarejmethd,'ADJUST')== 1
		Test=which('compute_GD_feat');
		if isempty(Test)
			warning(['The ADJUST plugin is required for ICA rejection using ', p.icarejmethd, '. APPLESEED is installing this plugin now.']);
			plugin_askinstall('ADJUST', [], true);
		%else %edited 20220321
			%!fid=fopen(which('MARA_extract_time_freq_features'));
		end
	end
	clear Test

	
	%FASTER (if running chaninterp)
	if strcmp(p.chaninterp, 'on') == 1
		Test=which('channel_properties');
		if isempty(Test)
			warning('The FASTER plugin is required for channel interpolation. APPLESEED is installing this plugin now.')
			plugin_askinstall('FASTER', [], true);
		end
	end
	clear Test



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            %
%    Initialize Analysis     %
%                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	eeg_global;
	
	%Load prepared dataset
	CurrentFileName=filenamebase;
	EEG = pop_loadset('filename',[CurrentFileName, '.set'],'filepath',inputdir); 


	%Check for channel file - required only for FASTER
	if strcmp(p.chaninterp, 'on')
		if ~isfield(EEG.chaninfo, 'filename') 
			if ~isfield(p, 'chanfile')
				error('Oh no! APPLESEED has FAILED! A channel location file is required for channel interpolation via FASTER.\n%s','See https://eeglab.org/tutorials/04_Import/Channel_Locations.html'); 
			else
				if exist(p.chanfile, 'file')==0
					error(['Oh no! APPLESEED has FAILED! The channel location file does not exist: ', p.chanfile]);
				else
					EEG=pop_chanedit(EEG, 'lookup',p.chanfile);
				end %if file doesn't exist
			end %if no defined file
		else
			p.chanfile=EEG.chaninfo.filename;
		end %if no set file
        
        %Check that faster re-ref channel exists
        if isempty(find(strcmpi({EEG.chanlocs.labels},p.fastref),1))
           error(['Oh no! APPLESEED has FAILED! ''',p.fastref, ''' was selected as the reference channel for channel interpolation, but this channel does not exist in the location file specified: ', p.chanfile])
        else
            p.fastrefn=find(strcmpi({EEG.chanlocs.labels},p.fastref));
        end
            
	end %if chaninterp


	%Check that sampling rate is not lower than resampling rate
	if EEG.srate <= p.resamp
		p.resampi=p.resamp;
		p.resamp=EEG.srate;
	end	


	%Establish final file name [Edits made 20221014] 	
	FinalFileName=[filenamebase, '_srate-', num2str(p.resamp), '_filt-', num2str(p.hp), '-', num2str(p.lp)];
	if isfield(p,'binfilename')
		FinalFileName=[FinalFileName, '_', strrep(p.binfilename, '.txt', ''), '-be0-', num2str(p.eplen)];
	end
	if strcmp(p.runica,'on') 
		FinalFileName=[FinalFileName,'_arxt-',num2str(p.arxt),'_ica-', char(strrep(p.icarejmethd, "_", ""))];
	end
    if strcmp(p.arauto, 'on')
        FinalFileName=[FinalFileName,'_arsd-', num2str(p.arsd)];
    else
        FileNameARAuto=FinalFileName;
        FinalFileName=[FinalFileName,'_arauto'];
    end
	if strcmp(p.chaninterp, 'on')
		FinalFileName=[FinalFileName,'_chaninterp'];
	end
    if strcmp(p.reref, 'Average')
        FinalFileName=[FinalFileName,'_reref-Avg'];
    else
        FinalFileName=[FinalFileName,'_reref-',char(strrep(p.reref, " ", "-"))];
    end
        
    FinalFileName=[FinalFileName,'_trialsel-', p.trselmethd, num2str(p.trselcnt)];
	logfilename=fullfile(outdir, [FinalFileName, '_logfile.txt']);
	
	%Check if task/rest and if analysis already run on the input dataset
	if ~isfield(p,'binfile') %rest
		%Check rest data for stimulus codes
		if contains([EEG.event.code], 'Stimulus')
			p.taskwarn=length(strfind([EEG.event.code],'Stimulus'));
		end

		%Check if rest analysis already run on the input dataset
		if size(dir(fullfile(outdir,[FinalFileName,'_scale-wise_entropy.txt'])),1) > 0
			disp(repmat('%',1,100))
			disp(['% APPLESEED already completed. See analysis log file: ', logfilename]);
			disp(repmat('%',1,100))
			close all;
			return;%!quit;
		end

	else %task

		%Check if task analysis already run on the input dataset
		if size(dir(fullfile(outdir,[FinalFileName,'*scale-wise_entropy.txt'])),1) == length(strfind(fileread(p.binfile),'bin'))
			disp(repmat('%',1,100))
			disp(['% APPLESEED already completed. See analysis log file: ', logfilename]);
			disp(repmat('%',1,100))
			close all;
			return;%!quit;
		end
	end %if task/rest
    
    %Check for manual AR [EDIT: 20221014]
    if strcmp(p.arauto, 'off') %run armanual
        if exist(fullfile(outdir,[FileNameARAuto,'.set']),'file') %already ran up to arrej
            if isfield(p, 'arfile') %if file specified
                RunStep=0; %prepare to coninue
                logfile=fopen(logfilename,'a');
            else
                error('ERROR - ready to procede with manual AR rejection but no AR file specified. Please specify before continuing')
                %return;
            end
        else
            RunStep=1;
            logfile=fopen(logfilename,'w');
        end
    else %run arauto
        RunStep=1;
        logfile=fopen(logfilename,'w');
    end


	%Initiate logfile
    if RunStep == 1
        fprintf(logfile,'%s\n\n  %s\n               %s\n               %s\n','Automated Preprocessing Pipe-Line for the Estimation of Scale-wise Entropy from EEG Data (APPLESEED)', ...
            '-Citations : Puglia, M.H., Slobin, J.S., Williams, C.L., 2022. The Automated Preprocessing Pipe-Line for the Estimation of Scale-wise Entropy from EEG Data (APPLESEED): Development and validation for use in pediatric populations. Developmental Cognitive Neuroscience, 101163.',...
            'Delorme, A., Makeig, S., 2004. EEGLAB: An open source toolbox for analysis of single-trial EEG dynamics including independent component analysis. J. Neurosci. Methods 134, 9Ð21.', ...
            'Lopez-Calderon, J., Luck, S.J., 2014. ERPLAB: an open-source toolbox for the analysis of event-related potentials. Front. Hum. Neurosci. 8, 213.');

        if strcmp(p.runica, 'on') && strcmp(p.icarejmethd, 'adjusted_ADJUST') 
            fprintf(logfile,'               %s\n               %s\n','Debnath, R., Buzzell, G.A., Morales, S., Bowers, M.E., Leach, S.C., Fox, N.A., 2020. The Maryland analysis of developmental EEG (MADE) pipeline. Psychophysiology 57, e13580.', ... 
                'Mognon, A., Jovicich, J., Bruzzone, L., Buiatti, M., 2011. ADJUST: An automatic EEG artifact detector based on the joint use of spatial and temporal features. Psychophysiology 48, 229Ð240.'); 
        end
        if strcmp(p.runica, 'on') && strcmp(p.icarejmethd, 'ADJUST'); fprintf(logfile,'               %s\n','Mognon, A., Jovicich, J., Bruzzone, L., Buiatti, M., 2011. ADJUST: An automatic EEG artifact detector based on the joint use of spatial and temporal features. Psychophysiology 48, 229Ð240.'); end
        if strcmp(p.chaninterp, 'on'); fprintf(logfile,'               %s\n','Nolan, H., Whelan, R., Reilly, R.B., 2010. FASTER: Fully Automated Statistical Thresholding for EEG artifact Rejection. J. Neurosci. Methods 192, 152Ð162.'); end
        fprintf(logfile,'\n%s\n\nAnalysis initiated %s\n\nInput dataset    : %s\n',repmat('-',1,200),datetime,fullfile(inputdir,[filenamebase, '.set']));
        if isfield(p,'chanfile'); fprintf(logfile,'Channel file     : %s\n',fullfile(p.chanfile)); end
        if ~isfield(p,'binfile') && isfield(p, 'taskwarn')
            fprintf(logfile,'\nWARNING: No bin file was specified, but %g stimulus event(s) were detected in the dataset.\n         %s\n\n', p.taskwarn,'Proceeding as if resting-state data...'); 
            warning(['No bin file has been specified, but there are ', num2str(p.taskwarn),' stimulus events in the dataset. Proceeding as if resting state data.']);
        elseif isfield(p,'binfile')
        fprintf(logfile,'Bin file         : %s\n\n',p.binfile);
        end
        fprintf(logfile,'Output directory : %s\n\n%s\n\n', fullfile(outdir),repmat('-',1,200));
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            %
%   Complete Preprocessing   %
%                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if RunStep == 1
        fprintf(logfile,'%s\n\n', 'Preprocessing: ');

        %%%%%%%%%%%%
        % Resample %
        %%%%%%%%%%%%
        FileNameAppend=['_srate-', num2str(p.resamp)];
        CurrentFileName=[CurrentFileName,FileNameAppend];

        if isfield(p, 'resampi')
            fprintf(logfile,'  -Resampling not completed.\n     Original sampling rate (%g Hz) is less than user input (%g Hz)\n\n',EEG.srate,p.resampi);
        else
            EEG = pop_resample( EEG, p.resamp);
            [ALLEEG,EEG,~] = pop_newset(ALLEEG, EEG, 0,'setname', CurrentFileName,'gui','off'); %#ok ALLEEG
            fprintf(logfile,'  -Resampling to %g Hz.\n\n',p.resamp);
        end


        %%%%%%%%%%
        % Filter %
        %%%%%%%%%%
        FileNameAppend=['_filt-', num2str(p.hp), '-', num2str(p.lp)];
        CurrentFileName=[CurrentFileName,FileNameAppend];

        EEG  = pop_basicfilter( EEG,  1:EEG.nbchan , 'Boundary', -99, 'Cutoff', [p.hp p.lp], 'Design', 'butter', 'Filter', 'bandpass', 'Order',  2, 'RemoveDC', 'on' ); 
        [ALLEEG,EEG,~] = pop_newset(ALLEEG, EEG, 0,'setname', CurrentFileName,'gui','off');

        fprintf(logfile,'  -Band-pass filtering %g-%g Hz.\n\n',p.hp,p.lp);


        %%%%%%%%%
        % Epoch %
        %%%%%%%%%
        if isfield(p,'binfilename') ~= 0 %task-based analysis		

            FileNameAppend=['_', strrep(p.binfilename, '.txt', ''), '-be0-', num2str(p.eplen)];
            CurrentFileName=[CurrentFileName,FileNameAppend];
            fprintf(logfile,'  -Segmenting task-based data into 0-%g ms epochs time-locked to stimulus onset.\n',p.eplen);  

            %Create elist
            EEG = pop_creabasiceventlist( EEG , 'AlphanumericCleaning', 'on', 'BoundaryNumeric', { -99 }, 'BoundaryString', { 'boundary' } );

            %Assign bins & epoch: 
            [EEG, EVENTLIST, binOfBins, ~] = neobinlister2( EEG , p.binfile, 'none', 'none', [], [], 0);	
            EEG.EVENTLIST=EVENTLIST;

            if ~isempty(find(binOfBins==0,1))
                fprintf(logfile,'     %s\n     %s\n','FAILURE: not all bins were identified in the dataset','Successful epochs per bin:'); 


                for b = 1:length(binOfBins)
                    fprintf(logfile,'         bin %g:  %3g\n', b, binOfBins(b));
                end 
                fprintf(logfile,'\n%s\n\nAnalysis terminated %s\n\n%s',repmat('-',1,200),datetime,repmat('-',1,200));
                fclose(logfile);
                close all;

                fprintf('%s\n%s%s\n%s\n',repmat('%',1,100),'% APPLESEED completed. See analysis log file: ',logfilename,repmat('%',1,100));

                return;%!quit;
            end		


            EEG = pop_epochbin( EEG , [0  p.eplen],  'none');
            EEG = pop_rmbase( EEG, [],[]);
            [ALLEEG,EEG,~] = pop_newset(ALLEEG, EEG, 0,'setname', CurrentFileName,'gui','off');

            %Check if enough trials to continue
            a=hist([EEG.event.bini],1:length(find(unique([EEG.event.bini])>0))); %#ok hist


            %Write to logfile
            if isempty(find(unique([EEG.event.bini])>0,1)) || sum(a(1:length(find(unique([EEG.event.bini])>0)))>=p.trselcnt)~=length(find(unique([EEG.event.bini])>0))		
                fprintf(logfile,'     %s\n\n     %s\n          %s\n','FAILURE: insufficient epochs to continue.','Epochs in original dataset:','Bin      #');  

                for b = 1:length(find(unique([EEG.event.bini])>0))
                    fprintf(logfile, '         %3g  %6g\n', b, a(b));
                end
                fprintf(logfile, '         %s\n         Total %5g\n\n',repmat('-',1,12), sum(a));
                fprintf(logfile,'%s\n\nAnalysis terminated %s\n\n%s', repmat('-',1,200),datetime,repmat('-',1,200));
                fclose(logfile);
                close all;
                fprintf('%s\n%s%s\n%s\n',repmat('%',1,100),'% APPLESEED completed. See analysis log file: ',logfilename,repmat('%',1,100));
                return;%!quit;

            else		
                fprintf(logfile,'     %s\n          %s\n','Epochs in original dataset:','Bin      #');

                for b = 1:length(find(unique([EEG.event.bini])>0))
                    fprintf(logfile, '         %3g  %6g\n', b, a(b));
                end

                fprintf(logfile, '         %s\n         Total %5g\n\n',repmat('-',1,12), sum(a));

            end %if enough trials	


        else %rest-based analysis

            %Determine latencies of evenly-spaced epochs
            EEG = pop_creabasiceventlist( EEG , 'AlphanumericCleaning', 'on', 'BoundaryNumeric', { -99 }, 'BoundaryString', { 'boundary' } );
            NumNewEvents = floor((EEG.xmax+1/EEG.srate)/(p.eplen/1000));
            if NumNewEvents > 0
                %Create unique event label for new epochs
                Events=unique([EEG.event.type]);
                EventLabel=Events(end)+1;

                %Remove old events and reate new evenly-spaced epochs
                EEG=rmfield(EEG, 'event');

                for n = 1:NumNewEvents
                    EEG.event(n)=struct('bepoch', 0, 'bini', 1, 'binlabel', '""', 'bvmknum',n,'bvtime',[],'channel',0,'codelabel', '""', 'duration', 0, 'enable',1,'flag',0, 'item', n, 'latency', (p.eplen/1000)*(n-1)*EEG.srate+1, 'type', EventLabel);
                end
                [EEG.data,~,indices,~]= epoch(EEG.data, [EEG.event.latency], [0 p.eplen/1000]*EEG.srate);

                for n = 1:size(indices,1)
                    EEG.event(n).epoch=n;
                end
                EEG = eeg_checkset(EEG, 'eventconsistency');

                EEG = pop_rmbase( EEG, [],[]);		
                [ALLEEG,EEG,~] = pop_newset(ALLEEG, EEG, 0,'setname', 'Test','gui','off');

                %Check if enough trials to continue
                a=hist([EEG.event.bini],1:length(find(unique([EEG.event.bini])>0))); %#ok hist

            else

                a=0;

            end %if NumNewEvents > 0

            %Write to logfile
            fprintf(logfile,'  -Segmenting resting-state data into evenly-spaced %g ms epochs.\n',p.eplen);  
            if NumNewEvents == 0 || sum(a(1:length(find(unique([EEG.event.bini])>0)))>=p.trselcnt)~=length(find(unique([EEG.event.bini])>0))		
                fprintf(logfile,'     %s\n\n     %s\n          %s\n','FAILURE: insufficient epochs to continue.','Epochs in original dataset:','Bin      #');  
                for b = 1:length(find(unique([EEG.event.bini])>0))
                    fprintf(logfile, '         %3g  %6g\n', b, a(b));
                end
                fprintf(logfile, '         %s\n         Total %5g\n\n',repmat('-',1,12), sum(a));
                fprintf(logfile,'%s\n\nAnalysis terminated %s\n\n%s',repmat('-',1,200),datetime,repmat('-',1,200));
                fclose(logfile);
                fprintf('%s\n%s%s\n%s\n',repmat('%',1,100),'% APPLESEED completed. See analysis log file: ',logfilename,repmat('%',1,100));
                close all;
                return;%!quit;

            else		
                fprintf(logfile,'     %s\n          %s\n','Epochs in original dataset:','Bin      #');

                for b = 1:length(find(unique([EEG.event.bini])>0))
                    fprintf(logfile, '         %3g  %6g\n', b, a(b));
                end

                fprintf(logfile, '         %s\n         Total %5g\n\n', repmat('-',1,12),sum(a));
            end %if enough trials

        end %if task/rest


        if strcmp(p.runica,'on') == 1 	

            %%%%%%%%%%%%%%%%
            % Step 5: ARXT %
            %%%%%%%%%%%%%%%%

            FileNameAppend=['_arxt-', num2str(p.arxt)];
            CurrentFileName=[CurrentFileName,FileNameAppend];

            %Count epochs
            prea=groupcounts([EEG.event.bini]')';
            %![prea,~]=hist([EEG.event.bini],1:length(find(unique([EEG.event.bini])>0)));
            EEG  = pop_artextval( EEG , 'Channel',  1:EEG.nbchan, 'Flag',  1, 'Threshold', [ -p.arxt p.arxt], 'Twindow', [ 0 p.eplen] );
            fprintf(logfile,'  %s\n     %s %g uV\n','-Rejecting artifacts via extreme voltage.','Threshold: >= +/-',p.arxt);

            %Check if enough trials to continue
            if length(find(EEG.reject.rejmanual == 1)) == EEG.trials
                fprintf(logfile,'     FAILURE: 100%% of epochs rejected. Insufficient epochs to continue\n\n%s\n\nAnalysis terminated %s\n\n%s',repmat('-',1,200),datetime,repmat('-',1,200));
                fclose(logfile);		
                fprintf('%s\n%s%s\n%s\n',repmat('%',1,100),'% APPLESEED completed. See analysis log file: ',logfilename,repmat('%',1,100));
                close all;
                return;%!quit;
            end %if enough trials

            [ALLEEG,EEG,~] = pop_newset(ALLEEG, EEG, 0,'setname', CurrentFileName,'gui','off');
            if strcmp(p.saverobust, 'on') == 1
                EEG = pop_saveset(EEG, 'filename', [CurrentFileName,'.set'], 'filepath', outdir);
            end %if save
            EEG = pop_rejepoch( EEG, find(EEG.reject.rejmanual == 1) ,0);


            %Count epochs 
            a=hist([EEG.event.bini],1:length(find(unique([EEG.event.bini])>0))); %#ok hist

            %Check if enough trials to continue
            if sum(a(1:length(find(unique([EEG.event.bini])>0)))>=p.trselcnt)~=length(find(unique([EEG.event.bini])>0))		
                fprintf(logfile,'     %s\n\n     %s\n          %s\n','FAILURE: insufficient epochs to continue.','Epochs after artifact rejection:','Bin   #(%) accepted   #(%) rejected');  
                for b = 1:length(find(unique([EEG.event.bini])>0))
                    fprintf(logfile, '         %3g  %6g(%5s) %8g(%5s)\n', b, a(b), sprintf('%.1f', a(b)/prea(b)*100), prea(b)-a(b), sprintf('%.1f',100-a(b)/prea(b)*100));			
                end							
                fprintf(logfile, '         %s\n         Total %5g(%5s) %8g(%5s)\n\n',repmat('-',1,35), sum(a), sprintf('%.1f', sum(a(b))/sum(prea(b))*100), sum(prea)-sum(a), sprintf('%.1f', (sum(prea)-sum(a))/sum(prea)*100));

                fprintf(logfile,'%s\n\nAnalysis terminated %s\n\n%s', repmat('-',1,200),datetime,repmat('-',1,200));
                fclose(logfile);
                close all;
                fprintf('%s\n%s%s\n%s\n',repmat('%',1,100),'% APPLESEED completed. See analysis log file: ',logfilename,repmat('%',1,100));
                return;%!quit;

            else
                %Write logfile
                if strcmp(p.saverobust, 'on') == 1
                    fprintf(logfile,'     Output file: %s.set\n',CurrentFileName);
                end %if save
                fprintf(logfile,'     %s\n          %s\n','Epochs after artifact rejection:','Bin   #(%) accepted   #(%) rejected');
                for b = 1:length(find(unique([EEG.event.bini])>0))
                    fprintf(logfile, '         %3g  %6g(%5s) %8g(%5s)\n', b, a(b), sprintf('%.1f', a(b)/prea(b)*100), prea(b)-a(b), sprintf('%.1f',100-a(b)/prea(b)*100));			
                end %for b			

                fprintf(logfile, '         %s\n         Total %5g(%5s) %8g(%5s)\n\n',repmat('-',1,35), sum(a), sprintf('%.1f', sum(a(b))/sum(prea(b))*100), sum(prea)-sum(a), sprintf('%.1f', (sum(prea)-sum(a))/sum(prea)*100));

            end %if enough trials


            %%%%%%%%%%%
            % Run ICA %
            %%%%%%%%%%%
            FileNameAppend='_ica';
            CurrentFileName=[CurrentFileName,FileNameAppend];
            fprintf(logfile,'  %s\n','-Running ICA.');

            EEG = pop_rmbase( EEG, []); %baseline zero (necessary for discontinuous epochs)
            EEG = pop_runica(EEG, 'extended',1, 'pca', EEG.nbchan);


            [ALLEEG,EEG,~] = pop_newset(ALLEEG, EEG, 0,'setname', CurrentFileName,'gui','off');
            if strcmp(p.saverobust, 'on') == 1
                EEG = pop_saveset(EEG, 'filename', [CurrentFileName,'.set'], 'filepath', outdir);
            end %if save

            if strcmp(p.saverobust, 'on') == 1
                fprintf(logfile,'     Output file: %s.set\n\n',CurrentFileName);			
            end %if save


            %%%%%%%%%%%%%%%%%%%%%%%%%
            % Reject Bad Components %
            %%%%%%%%%%%%%%%%%%%%%%%%%
            FileNameAppend=['-', char(strrep(p.icarejmethd, "_", ""))];
            CurrentFileName=[CurrentFileName,FileNameAppend];
            %Determine what components are bad based on rejection method


            if strcmp(p.icarejmethd,'manual') == 1
                pop_selectcomps(EEG,1:size(EEG.icawinv,2));
                uiwait;
                BadComps= find(EEG.reject.gcompreject == 1);
                p.icarejmethd='manual inspection';
            else	
                BadComps=eval([p.icarejmethd,'(EEG, ''',fullfile(outdir, [CurrentFileName, '_report.txt''']),')']);

                %Rename/move adjusted_ADJUST output image file
                par_id=strsplit(EEG.setname,'_');
                if exist([par_id{1}, '.jpg'], 'file')
                    movefile([par_id{1}, '.jpg'], fullfile(outdir, [CurrentFileName, '.jpg']))
                end
            end %if icarejmethd


            %Remove bad components 
            EEG = pop_subcomp( EEG, BadComps, 0); 

            %Write to log
            [ALLEEG,EEG,~] = pop_newset(ALLEEG, EEG, 0,'setname', CurrentFileName,'gui','off');
            if strcmp(p.saverobust, 'on') || strcmp(p.arauto, 'on')
                EEG = pop_saveset(EEG, 'filename', [CurrentFileName,'.set'], 'filepath', outdir);
            end %if save
            fprintf(logfile,'  -Removing bad components identified via %s.\n     %g component(s) removed\n\n',p.icarejmethd,length(BadComps));

        end %if runica	
	end %if runstep = 1


	%%%%%%%%
	% ARSD %
	%%%%%%%%
	if strcmp(p.arauto, 'off') %[Added 20221014]
        if ~isfield(p, 'arfile')
            fprintf('%s\n%s\n%s\n%s%s\n%s\n',repmat('%',1,100),'APPLESEED preprocessing is complete up to artifcat rejection.','You have selected manual artifact rejection. You must supply an artifact rejection file and re-run APPLESEED to continue.','Preprocessed dataset: ',fullfile(outdir,CurrentFileName),repmat('%',1,100));
            return;
        else
            %Load preprocessed dataset
            EEG = pop_loadset('filename',[FileNameARAuto, '.set'],'filepath',outdir); 
            CurrentFileName=FileNameARAuto;
            
            %Mark artifacts
            FileNameAppend='_arauto';
            CurrentFileName=[CurrentFileName,FileNameAppend];
            arfid=fopen(p.arfile, 'r');
            EEG.reject.rejmanual=fscanf(arfid,'%d')';
            if size(EEG.reject.rejmanual,2) ~= size(EEG.event,2)
                error(['Oh no! APPLESEED has FAILED! The EEG dataset contains ', num2str(size(EEG.event,2)), ' epochs, but the artifact rejection file ''', p.arfile, ''' contained rejection marks for ', num2str(size(EEG.reject.rejmanual,2)), ' epochs.']);
            else
                fprintf(logfile,'  %s\n     %s %s\n','-Rejecting artifacts via manual identification.','Artifact rejection file:',fullfile(outdir,p.arfile));
            end %if unexpected number of ar marks
        end %if no arfile specified
    else
            
        FileNameAppend=['_arsd-', num2str(p.arsd)];
        CurrentFileName=[CurrentFileName,FileNameAppend];
       
    	%Artifact detection SD: (Based on "Neural correlates of infants' sensitivity to vocal expressions of peers")	
		testwindow =  [0  p.eplen];
		winms      =  200;
		stepms     =  100;	
		fs       = EEG.srate;

		ntrial   = EEG.trials;
		interARcounter = zeros(1,ntrial);
		interARcounterE = zeros(EEG.nbchan,ntrial);
	
		winpnts  = floor(winms*fs/1000);
		stepnts  = floor(stepms*fs/1000);

		Ktime  = 1000;
		internum = testwindow/Ktime;

		toffsa = abs(round(EEG.xmin*fs))+1;

		p1 = round(internum(1)*fs) + toffsa;
		p2 = round(internum(2)*fs) + toffsa;


		chanArray= 1:EEG.nbchan;
		nch      = length(chanArray);

		for ch=1:nch
			for k=1:ntrial
				for j=p1:stepnts:p2-(winpnts-1)
					w1  = EEG.data(chanArray(ch), j:j+winpnts-1,k);
					vs = std(w1);
					if vs>p.arsd
							interARcounter(k) = 1 ;    
							interARcounterE(chanArray(ch),k) = 1;
					end
				end
			end
		end

		%Mark artifacts
		EEG.reject.rejmanual=interARcounter;
		EEG.reject.rejmanualE=interARcounterE;
		fprintf(logfile,'  %s\n     %s %g uV\n','-Rejecting artifacts via SD moving window.','Threshold:>= +/-',p.arsd);

	end %if autoar
		
    %Check if enough trials to continue
    if length(find(EEG.reject.rejmanual == 1)) == EEG.trials
        fprintf(logfile,'     FAILURE: 100%% of epochs rejected. Insufficient epochs to continue\n\n%s\n\nAnalysis terminated %s\n\n%s',repmat('-',1,200),datetime,repmat('-',1,200));
        fclose(logfile);
        fprintf('%s\n%s%s\n%s\n',repmat('%',1,100),'% APPLESEED completed. See analysis log file: ',logfilename,repmat('%',1,100));
        close all;
        return;%!quit;
    end %if enough trials


    %Count epochs
    [prea,~]=hist([EEG.event.bini],1:length(find(unique([EEG.event.bini])>0))); %#ok hist

    %Save file
    [ALLEEG,EEG,~] = pop_newset(ALLEEG, EEG, 0,'setname', CurrentFileName,'gui','off');
    if strcmp(p.saverobust, 'on') == 1
        EEG = pop_saveset(EEG, 'filename', [CurrentFileName,'.set'], 'filepath', outdir);
    end %if save

    %Reject Artifacts
    EEG = pop_rejepoch( EEG, find(EEG.reject.rejmanual == 1) ,0);

    %Count epochs
    a=hist([EEG.event.bini],1:length(find(unique([EEG.event.bini])>0))); %#ok hist

    %Check if enough trials to continue
    if sum(a(1:length(find(unique([EEG.event.bini])>0)))>=p.trselcnt)~=length(find(unique([EEG.event.bini])>0))		

        fprintf(logfile,'     %s\n\n     %s\n          %s\n','FAILURE: insufficient epochs to continue.','Epochs after artifact rejection:','Bin   #(%) accepted   #(%) rejected');  
        for b = 1:length(find(unique([EEG.event.bini])>0))
            fprintf(logfile, '         %3g  %6g(%5s) %8g(%5s)\n', b, a(b), sprintf('%.1f', a(b)/prea(b)*100), prea(b)-a(b), sprintf('%.1f',100-a(b)/prea(b)*100));			
        end	%for b		
        fprintf(logfile, '         %s\n         Total %5g(%5s) %8g(%5s)\n\n',repmat('-',1,35), sum(a), sprintf('%.1f', sum(a(b))/sum(prea(b))*100), sum(prea)-sum(a), sprintf('%.1f', (sum(prea)-sum(a))/sum(prea)*100));
        fprintf(logfile,'%s\n\nAnalysis terminated %s\n\n%s', repmat('-',1,200),datetime,repmat('-',1,200));
        fclose(logfile);
        fprintf('%s\n%s%s\n%s\n',repmat('%',1,100),'% APPLESEED completed. See analysis log file: ',logfilename,repmat('%',1,100));
        close all;
        return;%!quit;

    else

        if strcmp(p.saverobust, 'on') == 1
            fprintf(logfile,'     Output file: %s.set\n',CurrentFileName);
        end %if save
        fprintf(logfile,'     %s\n          %s\n','Epochs after artifact rejection:','Bin   #(%) accepted   #(%) rejected');
        for b = 1:length(find(unique([EEG.event.bini])>0))
            fprintf(logfile, '         %3g  %6g(%5s) %8g(%5s)\n', b, a(b), sprintf('%.1f', a(b)/prea(b)*100), prea(b)-a(b), sprintf('%.1f',100-a(b)/prea(b)*100));			
        end	%for b		
        fprintf(logfile, '         %s\n         Total %5g(%5s) %8g(%5s)\n\n',repmat('-',1,35), sum(a), sprintf('%.1f', sum(a(b))/sum(prea(b))*100), sum(prea)-sum(a), sprintf('%.1f', (sum(prea)-sum(a))/sum(prea)*100));

    end %if enough trials



	%%%%%%%%%%%%%%%%%%%%%%%%%
	% Channel Interpolation %
	%%%%%%%%%%%%%%%%%%%%%%%%%
	if strcmp(p.chaninterp, 'on') == 1

		FileNameAppend='_chaninterp';
		CurrentFileName=[CurrentFileName,FileNameAppend];

		list_properties = channel_properties(EEG, 1:EEG.nbchan, p.fastrefn); % run faster
		FASTbadIdx=min_z(list_properties);
		FASTbadChans=find(FASTbadIdx==1);
		FASTbadChans=FASTbadChans(FASTbadChans~=p.fastrefn);
	
		EEG = pop_interp(EEG, FASTbadChans, 'spherical');

		%Save dataset
		[ALLEEG,EEG,~] = pop_newset(ALLEEG, EEG, 0,'setname', CurrentFileName,'gui','off');

		fprintf(logfile,'  -Interpolating bad channels identified via FASTER with %s reference.\n     %g channel(s) interpolated\n',p.fastref,length(FASTbadChans));
		if length(FASTbadChans) >= 1
			BadChans={EEG.chanlocs(FASTbadChans).labels};
			for b = 1:length(BadChans)
				fprintf(logfile, '         %s\n', BadChans{b});
			end
		end
		fprintf(logfile,'\n');

	end %if chaninterp


	%%%%%%%%%%%%%%%%
	% Re-reference %
	%%%%%%%%%%%%%%%%	
	FileNameAppend='_avgref';
	CurrentFileName=[CurrentFileName,FileNameAppend];
	if exist(fullfile(outdir, [CurrentFileName, '.set']), 'file')== 0 %if file does not exist
		
		if strcmp(p.reref, 'Average')
			%Compute average ref
			EEG = pop_reref( EEG, []);	
		else
			rerefs=[];		
			for r=1:size(p.reref,2)
				rerefs=horzcat(rerefs, find(strcmpi({EEG.chanlocs.labels},p.reref{r}))); %#ok grow
			end			
			%Re-reference to specified channel(s)
			EEG = pop_reref( EEG, sort(rerefs)); 
		end %if avg reref
		
		%Save dataset
		[ALLEEG,EEG,~] = pop_newset(ALLEEG, EEG, 0,'setname', CurrentFileName,'gui','off');
		fprintf(logfile,'  -Re-referencing to %s.\n',p.reref);
        if strcmp(p.saverobust, 'on') == 1
            EEG = pop_saveset(EEG, 'filename', [CurrentFileName,'.set'], 'filepath', outdir);
            fprintf(logfile,'     Output file: %s.set\n',CurrentFileName);
        end %if save
	end %if file exists

	
	%%%%%%%%%%%%%%%%%%%
	% Trial Selection %
	%%%%%%%%%%%%%%%%%%%
	FileNameAppend=['_trialsel', p.trselmethd, num2str(p.trselcnt)];
	CurrentFileName=[CurrentFileName,FileNameAppend];

	
	if strcmp(p.trselmethd,'gfp')==1
	%Select trials: GFP 			

		%Calculate Total SS
		TotSS=[];
		for j = 1:size(EEG.data, 3)
			x=reshape(EEG.data(:,:,j), 1,numel(EEG.data(:,:,j)));
			P1=sum(x.^2);
			P2=(sum(x).^2)/size(x,2);
			TotSS=horzcat(TotSS,P1-P2); %#ok grow
		end %for trials


		TrialsToKeep=[];
		for b = 1:length(find(unique([EEG.event.bini])>0))
			TrialSelect{b}.Trials=find([EEG.event.bini]==b); %#ok grow
			if length(TrialSelect{b}.Trials)>p.trselcnt
				[~,TrialSelect{b}.SortedIndices]= sort(abs(TotSS(TrialSelect{b}.Trials)-median(TotSS))); %#ok grow
				TrialSelect{b}.SortedTrials=TrialSelect{b}.Trials(TrialSelect{b}.SortedIndices); %#ok grow
				TrialSelect{b}.KeepTrials=sort(TrialSelect{b}.SortedTrials(1:p.trselcnt)); %#ok grow
			else
				TrialSelect{b}.KeepTrials=find([EEG.event.bini]==b); %#ok grow
			end
			TrialsToKeep=horzcat(TrialsToKeep,TrialSelect{b}.KeepTrials); %#ok grow
		end %for bins


	elseif strcmp(p.trselmethd,'first')==1	
		
		TrialsToKeep=1:p.trselcnt;
	
	elseif strcmp(p.trselmethd,'middle')==1	
	
		TrialsToKeep=floor(size(EEG.data,3)/2-p.trselcnt/2+1:size(EEG.data,3)/2-p.trselcnt/2+p.trselcnt);
	
	elseif strcmp(p.trselmethd,'last')==1	
		
		TrialsToKeep=(size(EEG.data,3)-p.trselcnt+1:size(EEG.data,3));

	end %trialsel method


	%Discard excess trials
	EEG = pop_rejepoch( EEG, setdiff(1:size(EEG.data,3),TrialsToKeep), 0);
	[ALLEEG,EEG,~] = pop_newset(ALLEEG, EEG, 0,'setname', CurrentFileName,'gui','off');
	EEG = pop_saveset(EEG, 'filename', [CurrentFileName,'.set'], 'filepath', outdir);

	
	%Count remaining
	a=hist([EEG.event.bini],1:length(find(unique([EEG.event.bini])>0))); %#ok hist

	
	%Write to logfile
	fprintf(logfile,'\n  -Selecting trials: %s.\n     %s\n          %s\n',p.trselmethd,'Epochs in dataset after trial selection:','Bin      #');

	for b = 1:length(find(unique([EEG.event.bini])>0))
		fprintf(logfile, '         %3g  %6g\n', b, a(b));
	end

	fprintf(logfile, '         %s\n         Total %5g\n\n     Preprocessed dataset: %s.set\n\n',repmat('-',1,12), sum(a),CurrentFileName);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                %
%   Compute scale-wise entropy   %
%                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	fprintf(logfile,'%s\n\nCalculating scale-wise entropy:\n\n', repmat('-',1,200));

	if exist(fullfile(outdir, strrep([CurrentFileName, '_',EEG.EVENTLIST.bdf(length(find(unique([EEG.event.bini])>0))).description,'_scale-wise_entropy.txt'], '__', '_')),'file')== 0 %if file does not exist

		TotalBins=length(find(unique([EEG.event.bini])>0));
		for b = 1:length(find(unique([EEG.event.bini])>0))

			
			%Load data & remove trials not associated with current condition
			EEG = pop_loadset('filename',[CurrentFileName,'.set'],'filepath',outdir);
			[ALLEEG,EEG,~] = pop_newset(ALLEEG, EEG, 0,'setname', CurrentFileName,'gui','off');
			EEG = pop_rejepoch( EEG, setdiff(1:size(EEG.data,3),find([EEG.epoch.eventbini]==b)), 0);

			
			%Transform data to 1xnbchan cell array containing each trial as a column; calculate residuals
			for Chan = 1:size(EEG.data,1)
				Tmp=[];
				for Trial = 1:size(EEG.data,3)
					Tmp=horzcat(Tmp,EEG.data(Chan,:,Trial)'); %#ok grow
				end %for trial
				Raw{Chan}=Tmp; %#ok grow
				Mean{Chan}=mean(Raw{Chan},2); %#ok grow
				Res{Chan}=Raw{Chan}-Mean{Chan}; %#ok grow
			end %for chan		
			
			
			%Calculate entropy on Residuals	

			disp('Computing scale-wise entropy');		
			EntropyOut=[];
			for Chan = 1:EEG.nbchan
				
				%Create data structure for scale-wise entropy analysis
				for Trial = 1:size(EEG.data,3)
					Data{Trial}=Res{Chan}(:,Trial)'; %#ok grow
				end %for Trial

				Scales=floor(size(Data{1},2)/(p.m+1));
				SE = zeros(numel(Scales),1);

				for Scale = 1:Scales
				
					%Coarse-grain the time series
					DataCG=cellfun(@(x) mean(buffer(x, Scale), 1)',Data, 'UniformOutput',false);
							
					%Concatenate the data
					%This section of code computes entropy across discontinuous segments and has been adapted from Grandy et al. (2016).
					DataCat = [];
					for d = 1:length(DataCG)
						DataCat = [DataCat; DataCG{d}]; %#ok grow

						%Generate indices of points to match such that pattern-matching does not cross segment border
						if d == 1
							Points = (1:(length(DataCG{d}) - p.m))';
						elseif d > 1
							Points = [Points; (1:(length(DataCG{d}) - p.m))' + Points(end) + 2]; %#ok grow
						end
		
					end %for discontinuous segments
					
					%Generate indices of points to match across pattern length
					for MatchInds = 2:(p.m+1)
						Points(:,MatchInds) = Points(:,MatchInds-1)+1;
					end
								
					%Calculate r
					SD = std(DataCat); 
					r_new = p.r * SD;

					%Count pattern matches 
					Count = zeros(p.m+1,1);

					for j = 1:size(Points,1)
						for l = (j+1):size(Points,1)
							k = 1;
							while k <= p.m && abs(DataCat(Points(j,k)) - DataCat(Points(l,k))) <= r_new
								Count(k) = Count(k) + 1;
								k = k+1;
							end %while k
							if k == (p.m+1) && abs(DataCat(Points(j,k)) - DataCat(Points(l,k))) <= r_new
								Count(k) = Count(k)+1;
							end %if k
						end %for l
					end %for j

					%Compute sample entropy
					if Count(p.m+1) == 0 || Count(p.m) == 0
						SE(Scale) = -log(1/((size(Points,1))*(size(Points,1)-1)));
					else
						SE(Scale) = -log(Count(p.m+1)/Count(p.m));
					end %if no matches
				
				end %for scale

				%Save scale-wise entropy for each channel 
				EntropyOut=horzcat(EntropyOut,SE'); %#ok grow
			end %for channel
			

			%Write scale-wise entropy output
			entropyfile=fopen(fullfile(outdir,strrep([CurrentFileName, '_',strrep(deblank(char(EEG.EVENTLIST.bdf(b).description)), ' ', '-'),'_scale-wise_entropy.txt'], '__', '_')),'w');			
			Header=cat(2, 'Scale',{EEG.chanlocs.labels});
			fprintf(entropyfile, [repmat('%s\t',1,size({EEG.chanlocs.labels},2)+1), '\n'], Header{:});
			fprintf(entropyfile, ['%1d\t', repmat('%f\t',1,size({EEG.chanlocs.labels},2)),'\n'], horzcat((1:Scales)', EntropyOut)');
			fclose(entropyfile);
			
		end %for condition


		fprintf('%s\n%s%s\n%s\n',repmat('%',1,100),'% APPLESEED successfully completed! See analysis log file: ',logfilename,repmat('%',1,100));

		fprintf(logfile,'%s\n     Output file(s): %s\n','  -Computing scale-wise entropy with variance-normalization on data residuals.',strrep([CurrentFileName, '_',strrep(deblank(char(EEG.EVENTLIST.bdf(b).description)), ' ', '-'),'_scale-wise_entropy.txt'], '__', '_'));
		
		if TotalBins > 1
			for b = 1:length(find(unique([EEG.event.bini])>0))
				fprintf(logfile,'                     %s\n',strrep([CurrentFileName, '_',strrep(deblank(char(EEG.EVENTLIST.bdf(b).description)), ' ', '-'),'_scale-wise_entropy.txt'], '__', '_')); 
			end
		end
		fprintf(logfile,'\n     %s\n\n     %s\n          m=%g\n          r=%g\n          Scales: 1 (%g Hz) to %g (%g Hz)\n\n','Coarse graining procedure: original (moving window average)','Parameters:', p.m,p.r,EEG.srate,Scales,floor(EEG.srate/Scales));		

	else

		fprintf(logfile,'  -Scale-wise entropy with variance-normalization on data residuals already computed.\n     Output file(s): %s\n',fullfile(outdir, strrep([CurrentFileName, '_',strrep(deblank(char(EEG.EVENTLIST.bdf(b).description)), ' ', '-'),'_scale-wise_entropy.txt'], '__', '_')));		
		fprintf('%s\n%s%s\n%s\n',repmat('%',1,100),'% APPLESEED already completed. See analysis log file: ',logfilename,repmat('%',1,100));

	end %if file exists


	fprintf(logfile,'%s\n\nAnalysis terminated %s\n\n%s', repmat('-',1,200),datetime,repmat('-',1,200));
	fclose(logfile);
	close all;
	return;%!quit;

