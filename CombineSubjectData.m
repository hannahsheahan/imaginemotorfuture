%% FUNCTION: Combines subject data for each group into a single datafile
% INPUTS:  N/A
% OUTPUTS: N/A
% NOTES:    N/A
% ISSUES:   N/A
% REFS:     Wolpert analysis method
% AUTHOR:   Hannah Sheahan, sheahan.hannah@gmail.com
%--------------------------------------------------------------------------
%% settings
clear all; close all; clc;

%% specify subject data
ngroups   = 1;
nsubjects = 2; % per group
datName   = {
   
       
      % button-box + fixation + transfer test
      %{
        'OJ-EyeTrackingTransfer-12Feb2018-0945'; % balance
        'IC-EyeTrackingTransfer-balance-14Feb2018-1335'; % balance
        'FW-EyeTrackingTransfer-balance-16Feb2018-1335'; % balance
        'TL-EyeTrackingTransfer-balance-19Feb2018-1020'; % balance
        'OB-EyeTrackingTransfer-counter-12Feb2018-1330'; % counter
        'SL-EyeTrackingTransfer-counter-13Feb2018-0928'; % counter
        'TP-EyeTrackingTransfer-counter-16Feb2018-0940'; % counter
        'TM-EyeTrackingTranfer-counter-20Feb2018-1413'; % counter
        %}
        % Imagine follow-throughs
      %'EB-imagine-balance-16nov2016-1400'; % pilot subject, no button-box, instructed to fixate but no eye-tracking
      
      % button-box + free eye-movements (no eye-tracker)
      %{
      'MS-imagine-freeeye-21Nov2017-1325';          % balance
      'GG-imagine-balance-22Nov2017-1010';          % balance
      'HK-balance-dfreeeyemovement-29nov2017-1350'; % balance, storing data from misstrials (used slightly diff config)
      'GT-freeeyemovements-balance-21Dec2017-1318'; % balance
      'YQY-imagined-counterbalance-22Nov2017-1353'; % counter
      'ARH-imagined-counterbalance-24Nov2017-1331'; % counter       
      'SP-imagine-counterbalance-20Nov2017-0958';   % counter,  storing data from misstrials
      'AMF-freeyemovements-counterbalance-12Jan2018-0935'; % counter
      %}
      % button-box + fixation (with eye-tracker)
      %'EDDS-balance-28Nov2017-1415'; % balance  - lots of fixation errors but didn't record data from misstrials so cannot analyse if saccading to target
      
      % button-box + fixation (with eye-tracker) - NOTE that these 2
      % terminated the recording of eye-movements when a misstrial was
      % triggered
    %  'LA-balance-30Nov2017-0950'; % balance - recording misstrials
    %  'AC-fixation-counter-29nov2017-1704'; % counterbalance - recording misstrials
      
    % button-box + fixation (with eye-tracker and recording properly)
    %{
    'TM-balance-1Dec2017-1404'; % balance
    'YB-balance-5Dec2017-0939'; % balance
    'JK-balance-imageryfixation-17dec2017-1339'; % balance
    'BAB-imaginefollowthrough-balance-18Dec2017-1416'; % balance
    'RR-eyetracking-balance-4dec2017-1337'; % counterbalance
    'YO-eyetracking-counterbalance-14Dec2017-1019'; % counterbalance
    'GM-counter-imageryfixation-14dec2017-1717'; % counterbalance
    'OS-eyetracking-counterbalance-19Dec2017-1332'; % counterbalance
    %}
   % no imagery + fixation
   %{
   % 'FP-noimagerycontrol-09dec2017-1041'; % balance - exclude for imagining
    'JJM-noimagerycontrol-balance-10dec2017-1054'; % balance
    'EL-NoImageryControl-balance-13Dec2017-0944'; % balance
    'EY-noimagerycontrol-balance-15Dec2017-09210930'; % balance - found it hard to concentrate, error numbers are higher but not terrible. 
    'CH-NoImageryControl-balance-11Jan2017-0928'; % balance
    'ASM-noimagerycontrol-counter-09dec2017-231'; % counterbalance
    'FH-counter-noimagerycontrol-11dec2017-448'; % counterbalance
    'BM-NoImageryControl-counterbalance-12Dec2017-0929'; % counterbalance
    'RC-NoImageryControl-counterbalance-04Jan2018-0927';  % counterbalance
    %}
   
    % no imagery + fixation + transfer test

    'AB-NoImageryTransfer-balance-06Mar2018-1645'; % balance
    'ELG-NoImageryTransfer-counterbalance-07Mar2018-0931'; % counter

   
   };

MIQscore = [
      
    % button-box + fixation + transfer test
      % 4.71; 5.36; 4.86; 5.21; 5.14 ; 7; 4.57; 5.57;
       
      % button-box + free eye-movements (no eye-tracker)
      % 5.50; 5.36; 5.57; 4.43; 5.14; 3.93; 5.07; 6.64 
      
      % button-box + fixation (with eye-tracker)
      % 6.57; 6.07;   
      
      % button-box + fixation (with eye-tracker and recording properly)
     % 5.36; 6.21; 5.64; 4.71; 6.14; 6.43; 5; 5.71;
     
      % no imagery + fixation
     % 0;0;0;0;0;0;0;0;
     
     % no imagery + fixation + transfer
     0;0;
    
    ];

MIQvisual = [
    

    % button-box + fixation + transfer test
    % 5.14; 6.29; 4.57; 4.86; 5.57; 7; 5.71; 5.86;

    % button-box + free eye-mvements (no eye-tracker)
    %  5.71; 4.71; 6.00; 6.14; 6.14; 4.14; 5; 6.71;
      
    % button-box + fixation (with eye-tracker)
    % 6.71; 6.29;
    
     % button-box + fixation (with eye-tracker and recording properly)
    %  5.57; 6.43; 6.29; 5.14; 5.71;  6.28; 4.43; 5.43;
    
    % no imagery + fixation
    % 0;0;0;0;0;0;0;0;
  
    % no imagery + fixation + transfer
     0;0;
];

MIQkinesthetic = [
   
    % button-box + fixation + transfer test
    % 4.28; 4.43; 5.14; 5.57; 4.71; 7; 3.42; 5.29;

    % button-box + free eye-mvements (no eye-tracker)
    % 5.29; 6.00; 5.14; 2.71; 4.14; 3.71; 5.14; 6.57; 
    
    % button-box + fixation (with eye-tracker)
    % 6.43; 5.86;
    
     % button-box + fixation (with eye-tracker and recording properly)
    %  5.14; 6.00;  5.57; 4.29; 6.00; 6.57; 5.57; 6;
    
    % no imagery + fixation
    % 0;0;0;0;0;0;0;0;
    
    % no imagery + fixation + transfer
     0;0;
    
];

maintenanceMIscore = [
    
    % button-box + fixation + transfer test
    %{
      [2,1;5,2;7,3;5,2;6,2;4,2;2,2;4,2;6,2]; % OJ
      [6,2;6,2;6,2;6,2;6,2;6,2;6,2;6,2;6,2]; % IC
      [6,2;6,2;6,2;6,3;6,2;6,3;6,2;6,3;6,3] % FW
      [4,2;4,1;4,2;5,1;5,2;5,2;5,2;4,2;4,2]; % TL
      [5,2;4,1;5,1;5,2;6,2;6,2;6,2;6,2;6,2]; % OB
      [6,2;6,2;6,2;5,2;5,2;6,2;5,2;4,1;4,1]; % SL
      [4,2;5,2;5,2;5,3;5,3;5,3;5,3;6,3;6,2] % TP
      [3,1;5,2;5,2;3,1;5,3;6,2;6,2;3,1;6,2]; % TM
     %}
      
     % button-box + free eye-mvements (no eye-tracker)
     %{
      [4,2; 5,3; 5,3; 5,2; 5,3; 6,3; 5,3; 6,3]; % MS
      [4,2; 4,2; 5,2; 5,2; 4,2; 5,2; 5,2; 5,2]; % GG
      [6,2; 6,3; 5,3; 7,3; 6,2; 5,1; 6,2; 6,2]; % HK
      [5,2; 4,1; 5,2; 5,2; 5,2; 5,2; 5,2; 5,2]; % GT
      [3,2; 3,2; 4,2; 4,2; 3,2; 4,2; 5,2; 5,2]; % YQY
      [4,2; 5,2; 3,2; 4,2; 5,3; 5,3; 3,3; 4,2]; % ARH
      [6,2; 7,3; 6,2; 7,3; 7,3; 6,2; 7,3; 7,3]; % SP
      [6,2; 7,3; 6,2; 6,2; 6,2; 6,2; 7,3; 7,3]; % AMF
     %}
      % button-box + fixation (with eye-tracker)
      % [5,2; 5,1; 5,2; 6,3; 7,2; 7,2; 6,2; 6,2]; % LA
      % [3,1; 4,2; 3,2; 4,2; 3,2; 5,3; 4,2; 3,2]; % AC
      
       % button-box + fixation (with eye-tracker and recording properly)
      %{
       [3,2;3,1;3,1;3,1;1,1;2,1;1,1;3,2];
      [5,2;4,2;5,3;5,2;5,1;4,2;6,3;6,3];
      [3,2;3,2;4,2;4,2;4,2;4,2;4,2;5,2];
      [6,2;6,2;6,2;3,1;4,2;4,2;3,1;4,2];
      [3,1;4,2;4,2;4,2;5,2;5,2;4,2;3,2];
      [5,3;6,3;7,1;6,2;6,2;7,2;6,3;7,2];
      [3,1;4,1;4,2;4,2;5,2;5,2;6,2;5,3];
      [5.5,2;6,3;5.5,3;6,3;6.5,3;6,3;5,2;6.5,3];
      %}
    % no imagery + fixation
    %{
       [0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0];
       [0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0];
       [0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0];
       [0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0];
       [0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0];
       [0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0];
       [0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0];
       [0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0];
    %}
    
    % no imagery + fixation + transfer
       [0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0];
       [0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0];
    
    ];

%% specify parameters you want to save
grouplist = 1:ngroups;  % number subjects by group
group  = repmat(grouplist,nsubjects,1);
group  = group(:);
nfiles = size(datName,1);
%fdir = [1;1;1;1;-1;-1;-1;-1];
%recordfbtime = [0;0;1;1;0;0;1;1];  % data recording during feedback delay period
%fixation = 0;

fdir = [1;-1];
recordfbtime = [1;1];
fixation = 1;

D = {};

%loop over files to load data and remove elements you don't want
for k=1:nfiles
    sprintf('Subject file: %d...',k);
    ww=['../01 Import Data/data/' datName{k}]; %filenames with .dat/.mat stripped off
    tmp=DATAFILE_Load(ww);
    tmpmisstrials = tmp;

    % save misstrial data separately in datafile
    if isfield(tmp,'MissTrialFlag')
        selecttrials = find(~tmp.MissTrialFlag & tmp.FrameData.State(:,5)==4);
        selectmisstrials = find(tmp.MissTrialFlag);
    else
        selecttrials = find(tmp.FrameData.State(:,5)==4);  % exclude passive returns but include misstrials separately
        selectmisstrials = [];
    end
    
    tmp=DATAFILE_Select(tmp,selecttrials);
    tmpmisstrials = DATAFILE_Select(tmpmisstrials,selectmisstrials);
    
    % Remove the fields we don't want
    % Regular fields:
    fields = {'Files'...
        'MovementDurationTooSlowToVia'...
        'MovementDurationTooSlow'...
        'MovementDurationTooFast'...
        'ViaSpeedThreshold'...
        'ViaTimeOutTime'...
        'ViaToleranceTime'...
        'PMoveEndPosition'...
        'PMoveStartPosition'...
        };
    for j=1:length(fields)
        if isfield(tmp,fields(j))
            tmp = rmfield( tmp,fields(j));
            tmpmisstrials = rmfield( tmpmisstrials,fields(j));
        end
    end
    % Framedata fields:
    fields = {'GraphicsVerticalRetraceNextOffsetTime',...
        'GraphicsSwapBuffersToVerticalRetraceTime',...
        'GraphicsSwapBuffersCount',...
        'PMoveStatePosition',...
        'PMoveStateRampValue',...
        'CursorPosition',...
        'HandleTorques',...
        'HandleForces',...
        'ForcesFunctionPeriod',...
        'ForcesFunctionLatency',...
        'StateGraphics',...
        'GraphicsVerticalRetraceNextOnsetTime' ...
        'PMoveStateTime' ...
        'PMoveState' ...
        };
    for j=1:length(fields)
        tmp.FrameData = rmfield( tmp.FrameData,fields{j});
        tmpmisstrials.FrameData = rmfield( tmpmisstrials.FrameData,fields{j});
    end
    
    % add the group and subject info to the datafile
    tmp.group=group(k)*ones(size(tmp.TrialNumber))';
    tmp.subj=k*ones(size(tmp.TrialNumber))';
    tmp.recordfbtime=recordfbtime(k).*ones(size(tmp.TrialNumber))';
    tmpmisstrials.group=group(k)*ones(size(tmpmisstrials.TrialNumber))';
    tmpmisstrials.subj=k*ones(size(tmpmisstrials.TrialNumber))';
    D = DATAFILE_Append(D,tmp);
    if ~exist('MissTrialData','var')
        MissTrialData = [];
    end
    if isfield(tmp,'MissTrialFlag')
        MissTrialData = DATAFILE_Append(MissTrialData,tmpmisstrials);
    end
    clear tmp tmpmisstrials;
end

D.MissTrialData = MissTrialData;
D.datName=datName;
D.fdir=fdir;
D.MIscore = MIQscore;
D.MIQvisual = MIQvisual;
D.MIQkinesthetic = MIQkinesthetic;
D.maintenanceMIscore = maintenanceMIscore;
D.fixation = fixation;

% remove the z-dimension from the kinetic/kinematic data (keep files small)
D.FrameData.RobotForces(:,:,3)=[];
D.FrameData.RobotPosition(:,:,3)=[];
D.FrameData.RobotVelocity(:,:,3)=[];

%%
% save it all in a single datafile (much faster to load later)
%save -v7.3 Imagine_combined D
%save -v7.3 ImagineFixation_combined D
%save -v7.3 NoImagineFixation_combined D
%save -v7.3 ImagineFixationTransfer_combined D
save -v7.3 NoImageryFixationTransfer_combined D

