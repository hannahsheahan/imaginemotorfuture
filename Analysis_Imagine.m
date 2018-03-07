%% Motor Planning, Not Execution, Separates Motor Memories
% Sheahan, Franklin & Wolpert. (2016), Neuron 92(4)
% - This is an analysis script for assessing learning in 4 main groups of
%   subjects (n=6 per group), and a control group (n=4).

% Group 1: no imagery
% Group 2: imagery fixation
% Group 3: imagery no fixation
% Group 4: planning only
% Group 5: full follow-through

% Author: Hannah Sheahan, sheahan.hannah@gmail.com
% Date:   Jan 2018
%--------------------------------------------------------------------------

%% settings
clc; clear all; close all;
tic
FontSize = 14;
set(0, 'DefaultAxesFontSize',FontSize); clear FontSize;
set(0,'DefaultFigureWindowStyle','docked');
set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesTickDirMode', 'manual');
FLAGS.normMPE = 0;              % normalise MPE by group peak speed
FLAGS.plotextra = 0;            % plot additional figures (subject-by-subject analysis, hand paths, target appearance time distributions)

%% load subject data
ngroups = 5;
for group = 1:ngroups  % load in by group because 'planning only' and 'control' groups have different total no. trials
    
    switch group
        case 1
            load NoImagineFixation_combined;
        case 2
            load ImagineFixation_combined;  
        case 3
            load Imagine_combined;
        case 4
            load Planning8_combined;
        case 5
            load Fullfollowthrough8_combined;
    end
    %% define general variables
    
    % Frame, trial, subject, group numbers
    trial  = D.ResampTrialNumber';
    ntrial = max(trial);                        %number of trials for each subject
    N      = D.Trials;                          %total number of trials across all subjects
    ndata  = size(D.FrameData.RobotPosition,2); %number of samples (data points) in longest trial
    nsubj  = max(D.subj);
    
    % Field types
    S.FIELD_NULL     = 0 +1;
    S.FIELD_VISCOUS  = 1 +1;
    S.FIELD_CHANNEL  = 2 +1;
    S.FIELD_PMOVE    = 3 +1;
    
    % State names
    S.STATE_GO       = 5;
    S.STATE_MOVEWAIT = 6;
    S.STATE_MOVING0  = 7;
    S.STATE_MOVING1  = 8;
    S.STATE_FINISH   = 9;
    
    % Label groups
    G.NOIMAGERY     = 1;
    G.IMAGERY       = 2;
    G.IMAGERYNOFIX  = 3;
    G.PLANNING      = 4;
    G.FOLLOWTHROUGH = 5;
    
    
    % binary field-type matrix
    fieldlist = unique(D.FieldType);
    field = zeros(N,length(fieldlist)+1);
    for i=1:length(fieldlist)
        FN = fieldlist(i)+1;
        field(:,FN) = D.FieldType==fieldlist(i);
    end
    
    % Experiment phases
    if group==G.FOLLOWTHROUGH
        phaselist = cumsum([5 150 3]);
    else
        phaselist = cumsum([3 75 2]);
    end
    Baseline = D.PhaseIndex<=phaselist(1);  indB = find(Baseline==1);
    exposurephase = D.PhaseIndex>phaselist(1) & D.PhaseIndex<=phaselist(2);  indE = find(exposurephase==1);
    Post     = D.PhaseIndex>phaselist(2) & D.PhaseIndex<=phaselist(3);  indP = find(Post==1);
    
    clear i FN fieldlist phaselist Baseline Exposure Post
    usefulthings = ws2struct();             % save some general access variables
    frames = D.FrameData.Frames;
    state = D.FrameData.State;
    trialtime = D.FrameData.TrialTime;
    
    %----------------------------
    %% Kinematic variable processing
    
    % Reformat kinematic variables
    %   - pad kinematic variables with NaNs in the excess frames, ready for resampling
    posx   = squeeze(D.FrameData.RobotPosition(:,:,1));
    posy   = squeeze(D.FrameData.RobotPosition(:,:,2));
    velx   = squeeze(D.FrameData.RobotVelocity(:,:,1));
    vely   = squeeze(D.FrameData.RobotVelocity(:,:,2));
    forcex = squeeze(D.FrameData.RobotForces(:,:,1));
    forcey = squeeze(D.FrameData.RobotForces(:,:,2));
    
    ind1 = D.FrameData.Frames;
    ind2 = ndata-ind1;
    
    counter = repmat(D.fdir, [1,ntrial])';  % indicate which subjects are counterbalanced for field direction
    counter = counter(:);
    fdir = -sign(D.TargetAngle);            % use this to signal each force-field direction
    fdir = fdir.*counter;                   % multiply all *-1 if field association is counterbalanced
    usefulthings.fdir = fdir;
    
    %pad with Nans (since each trial lasts a different number of frames)
    padnan = @(x,ind1,ind2,N) cell2mat( arrayfun(@(k)([x(k,1:ind1(k)) NaN(1,ind2(k))]), 1:N, 'UniformOutput', false)');
    
    posx   = padnan(posx,ind1,ind2,N);
    posy   = padnan(posy,ind1,ind2,N);
    velx   = padnan(velx,ind1,ind2,N);
    vely   = padnan(vely,ind1,ind2,N);
    forcex = padnan(forcex,ind1,ind2,N);
    forcey = padnan(forcey,ind1,ind2,N);
    
    posfullx = posx;  % store the full trajectory position data for plotting later
    posfully = posy;  % store the full trajectory position data for plotting later
    
    % Remove data after viapoint
    ind2 = findfirst(D.FrameData.State >= S.STATE_MOVING1,2);
    ind1 = ones(size(ind2));
    pad  = ndata -(ind2-ind1) -1;
    
    chopnan = @(x,ind1,ind2,N) cell2mat(arrayfun(@(k)([x(k,ind1(k):ind2(k)) NaN(1,pad(k))]), 1:N, 'UniformOutput', false)');
    
    posx = chopnan(posx,ind1,ind2,N);
    posy = chopnan(posy,ind1,ind2,N);
    velx = chopnan(velx,ind1,ind2,N);
    vely = chopnan(vely,ind1,ind2,N);
    forcex = chopnan(forcex,ind1,ind2,N);
    forcey = chopnan(forcey,ind1,ind2,N);
    pathlength = nansum(sqrt(diff(posy,1,2).^2+diff(posx,1,2).^2),2);
    
    %----------------------------
    %% Resample kinematic variables (to 1000 samples per trial)
    nsamp = 1000;
    
    % keep some full trajectory data
    original.state = state;
    original.trialtime = trialtime;
    original.posy = posfully;
    original.posx = posfullx;
    original.lngth = size(state,2);
    
    % Resample data for first section of movement to via point
    ind1 = findfirst(state >= S.STATE_MOVING0,2);
    ind2 = findfirst(state >= S.STATE_MOVING1,2);
    
    resample = @(x,frames) interp1(1:length(x),x,linspace(1,length(x),frames));
    kinematic_resample = @(x,ind1,ind2,N) cell2mat( arrayfun(@(k) resample(x(k,ind1(k):ind2(k)),nsamp), 1:N, 'UniformOutput', false)');
    
    F.state     = kinematic_resample(state,ind1,ind2,N);
    F.trialtime = kinematic_resample(trialtime,ind1,ind2,N);
    posx      = kinematic_resample(posx,ind1,ind2,N);
    posy      = kinematic_resample(posy,ind1,ind2,N);
    velx      = kinematic_resample(velx,ind1,ind2,N);
    vely      = kinematic_resample(vely,ind1,ind2,N);
    forcex    = kinematic_resample(forcex,ind1,ind2,N);
    forcey    = kinematic_resample(forcey,ind1,ind2,N);
    
    % save a copy of variables before rotating (for later trajectory plotting)
    F.rpos = [reshape(posx',[1,nsamp,N]);reshape(posy',[1,nsamp,N])];
    F.rvel = [reshape(velx',[1,nsamp,N]);reshape(vely',[1,nsamp,N])];
    F.rforce = [reshape(forcex',[1,nsamp,N]);reshape(forcey',[1,nsamp,N])];
    original.pos = [reshape(original.posx',[1,original.lngth,N]);reshape(original.posy',[1,original.lngth,N])];
    original.pathlength = pathlength;
    
    clear  posx posy velx vely forcex forcey ind1 ind2 state trialtime pathlength
    
    %% Rotate kinematics from each start position
    startangle = -D.HomeAngle.*(pi/180);
    startpos   = repmat(reshape(D.StartPosition(:,1:2)', [2,1,N]),[1,nsamp,1]);
    original.startpos = repmat(reshape(D.StartPosition(:,1:2)', [2,1,N]),[1,original.lngth,1]);
    
    rotate = @(x,theta,N) cell2mat(( arrayfun(@(k)(reshape([cos(theta(k)), -sin(theta(k)); sin(theta(k)), cos(theta(k))]' * squeeze(x(:,:,k)),[2,1,size(x,2)])), 1:N, 'UniformOutput', false)));
    
    original.pos = rotate((original.pos-original.startpos),startangle,N);
    pos   = rotate((F.rpos-startpos),startangle,N);
    vel   = rotate(F.rvel,startangle,N);
    force = rotate(F.rforce,startangle,N);
    
    % rotate misstrial eye and hand data (not resampled)
    if group<G.IMAGERYNOFIX
        misstrialdata = D.MissTrialData;
        misstrialstartangle = -misstrialdata.HomeAngle.*(pi/180);
        viapos = misstrialdata.ViaPosition(:,[1,2]);
        viapos = permute(repmat(viapos,[1,1,size(misstrialdata.FrameData.EyeTrackerEyeXY,2)]),[1,3,2]);
        eyepos = rotate(permute((misstrialdata.FrameData.EyeTrackerEyeXY-viapos),[3,2,1]), misstrialstartangle, misstrialdata.Trials);
        handpos = rotate(permute((misstrialdata.FrameData.RobotPosition(:,:,[1:2])-viapos),[3,2,1]), misstrialstartangle, misstrialdata.Trials);
        F.mistrialeyepos = eyepos;
        F.mistrialhandpos = handpos;
    end
    
    F.posx = squeeze(pos(1,:,:));
    F.posy = squeeze(pos(2,:,:));
    F.velx = squeeze(vel(1,:,:));
    F.vely = squeeze(vel(2,:,:));
    F.forcex = squeeze(force(1,:,:));
    F.forcey = squeeze(force(2,:,:));
    F.rstartpos = unique(squeeze(startpos(:,1,:))','rows');
    original.posx = squeeze(original.pos(1,:,:));
    original.posy = squeeze(original.pos(2,:,:));
    
    %% Calculate MPE and adaptation
    
    % calculate mpe
    [test,ind] = max(abs(F.posx'));
    ind = sub2ind(size(F.posx),1:N,ind);
    mpe = F.posx(ind)'.*(fdir==1) - F.posx(ind)'.*(fdir==-1);     % accounts for counterbalance subjects and switching field directions
    
    % calculate maximum signed lateral deviation
    F.lateraldev = F.posx(ind)';
    
    % Calculate adaptation (% complete)
    ph1 = find(field(:,S.FIELD_VISCOUS)==1);
    fieldconstant = D.FieldConstants(ph1(1),1);                   % same for all subjects
    perffx = fieldconstant.*repmat(fdir,[1,1000]).*F.vely;        % ideal force. NB: we can only measure forces in the x-direction
    
    % find +-200ms window either side of time of max speed (NB: often takes entire length of trial, but limits very long trials)
    speed = sqrt((F.velx.^2)+(F.vely.^2));
    [~,vmaxi] = max(speed,[],2);                                  % index of max speed
    tsmpl = 0.2;
    lvmaxi = sub2ind(size(F.trialtime),[1:N]',vmaxi);
    tvmax  = F.trialtime(lvmaxi);
    t_low   = tvmax-tsmpl;                                        % limits of time interval
    t_up    = tvmax+tsmpl;
    for i=1:N
        smpl_L(i,1)  = length(find(F.trialtime(i,:) < tvmax(i) & F.trialtime(i,:) > t_low(i)));
        smpl_U(i,1)  = length(find(F.trialtime(i,:) > tvmax(i) & F.trialtime(i,:) < t_up(i)));
    end
    p1 = vmaxi-smpl_L;
    p2 = vmaxi+smpl_U;
    
    % regress actual force onto ideal force (zero intercept) for adaptation measure. (lambda/anonymous function here takes <half computing time of 'for' loop)
    adaptregress = @(x,y,p1,p2,N) cell2mat( arrayfun(@(k)(100.*regress(y(k,p1(k):p2(k))',x(k,p1(k):p2(k))') ), 1:N, 'UniformOutput', false));
    adaptation = adaptregress(perffx,F.forcex,p1,p2,N);
    
    %------------------------
    %% Calculate mistrials
    % The mistrial data is saved in different forms for different groups -
    % some groups save trials with breaks in fixation in a separate
    % misstrialdata structure, so search through both this structure and
    % regular trials for mistrials.
    
    % Total number of mistrials of all types, on trials without fixation breaks (not a major metric)
    for i=1:nsubj
        ind = [ntrial.*(i-1)+1,ntrial.*i];
        misstrialrate(i) = max(D.MissTrials(ind(1):ind(2))) ./ (max(D.MissTrials(ind(1):ind(2)))+ntrial);
    end
    
    tolerance = 3;   % cm. Same metric for hand overshoot and breaks in fixation.
   
    % percentage of overshot trials per subject on regular exposure trials (not fixation break trials)
    %(NB: overshot = if subject travels 3cm+ in the rotated +y direction after the viapoint)
    targetdist = 12;
    overdistance = original.posy - repmat(targetdist,size(original.posy));
    
    ch = (original.state >= S.STATE_MOVING1);                     % after hitting the via-point
    overshot = sign(sum(ch & (overdistance >= tolerance),2));     % binary overshot or not per trial
    overshot(field(:,S.FIELD_VISCOUS)~=1)=0;                      % ignore count on channel trials (subjects supposed to follow through) and null trials
    nincludedtrials = sum(field(1:ntrial,S.FIELD_VISCOUS)==1);    % same number of planning-only trials per subject (planning-only group only)
    for i=1:nsubj
        ind = [ntrial.*(i-1)+1,ntrial.*i];
        overshootrate(i) = sum(overshot(ind(1):ind(2))) ./ nincludedtrials;
        novershot(i) =  sum(overshot(ind(1):ind(2)));
    end
    timing.overshootrate = overshootrate;
    
    % Fixation break mistrials (for fixation groups only). (consider only exposure phase)
     nimagine = sum((field(:,S.FIELD_VISCOUS)==1) & D.subj==1);  % same for each subject, so use subj==1 as example
     
    if group<G.IMAGERYNOFIX 
        ind = misstrialdata.FieldType==1;                       % imagining trials in exposure phase (by induction since FieldType==1 means curl field)
        
        % measure percent overshoot relative to all attempts at imagining trials
        for i=1:nsubj
            nimaginesub(i) = nimagine + sum(ind & misstrialdata.subj==i);    
        end
        
        eyebrokefixation = findfirst((sqrt(eyepos(2,:,:).^2 + eyepos(1,:,:).^2) > tolerance),3);  % broke fixation by radial 3cm
        handovershot = findfirst(handpos(2,:,:) > tolerance,3);             % overshot of via target in rotated y direction
       
        % ignore trials that were blinks, where eyetracker loses eye, or if fixation broken before 'go' cue 
        blink = squeeze(findfirst(eyepos(1,:,:)>100 | eyepos(2,:,:)>100,3));
        eyebrokefixation(blink>0) = 0;
        eyebrokefixation(eyebrokefixation<300) = 0;                        % don't count instances where fixation is broken before 'go' cue (pre-300ms)
        
        % per subject misstrial errors on imagining trials
        for i=1:nsubj
            trials = find((misstrialdata.subj==i) & ind);
            nerrsaccadeovershoot(i) = sum(eyebrokefixation(trials)~=0);    % number of times eye broke fixation on imagining trials
            nerrhandovershoot(i) = sum(handovershot(trials)~=0);           % number of times hand overshot central target on imagining trials
        end
             
        % percentage of imagining trials with breaks in fixation 
        % (excl. blinks, fixation broken before 'go' cue)
        timing.percerreyebreakfixation = nerrsaccadeovershoot./nimaginesub;
        timing.nerrsaccadeovershoot = nerrsaccadeovershoot;
        F.misstrialdata = misstrialdata;
    else
        % a few subjects in the imagery no fixation group have misstrialdata saved separately, so make
        % sure this is included for hand overshoot data
        if isfield(D,'MissTrialData')
            misstrialdata = D.MissTrialData;
            ind = misstrialdata.FieldType==1;                              % imagining trials in exposure phase
            
            % measure percent overshoot relative to all attempts at imagining exposure trials
            for i=1:nsubj
                nimaginesub(i) = nimagine + sum(ind & misstrialdata.subj==i);
            end
            
            % find trials in which hand overshot via point by 3cm
            viapos = misstrialdata.ViaPosition(:,[1,2]);
            viapos = permute(repmat(viapos,[1,1,size(misstrialdata.FrameData.RobotPosition,2)]),[1,3,2]);
            handpos = rotate(permute((misstrialdata.FrameData.RobotPosition(:,:,[1:2])-viapos),[3,2,1]),startangle,misstrialdata.Trials);
            handovershot = findfirst(handpos > tolerance,2);   % hand overshot central target by 3cm+ in y direction
            
            % per subject misstrial errors on imagining/aborting trials
            for i=1:nsubj
                trials = find((misstrialdata.subj==i) & ind);
                nerrhandovershoot(i) = sum(handovershot(trials)~=0);       % number of times hand overshot on imagining trials
            end
        else
            nerrhandovershoot = zeros(nsubj,1);
            nimaginesub = nimagine;
        end
    end
    
    % percentage of imagining trials that eye-moves or hand-moves to
    % second target triggered errors on
    timing.percerrhandovershoot = (nerrhandovershoot+novershot)./nimaginesub;
    timing.nerrhandovershot = nerrhandovershoot + novershot;
    timing.misstrialrate = misstrialrate;
    
    %------------------------
    %% Calculate timing measures e.g. planning time, duration...
    
    % pre-movement time (from when target first appears)
    ph  = find((field(:,S.FIELD_CHANNEL)==1));
    appeardist = [0,0,10,6,0];
    appeardist = appeardist(group);
    
    % find duration and timing of target appearances
    iRX = sub2ind(size(original.state),[1:N],findfirst(original.state' >= S.STATE_MOVING0));
    iVia = sub2ind(size(original.state),[1:N],findfirst(original.state' >= S.STATE_MOVING1));
    iGo = sub2ind(size(original.state),[1:N],findfirst(original.state' >= S.STATE_GO));
    iFin = sub2ind(size(original.state),[1:N],findfirst((original.state >= S.STATE_GO),2,1,'last')');
    rxtime = original.trialtime(iRX);                                           % reaction time (incl. 300ms delay period)
    
    ind = sub2ind(size(original.posy),[1:N],findfirst(original.posy' >= appeardist));
    appeartime = original.trialtime(ind) - rxtime;                          % time of target appear after movement onset
    appeartimeaftergo = original.trialtime(ind) - original.trialtime(iGo);  % time of target appear after go cue
    
    % nan appearance times for channel trials (should be 0 for all except planning group)
    appeartime(ph) = nan;
    appeartimeaftergo(ph) = nan;
    
    duration = original.trialtime(iFin) - original.trialtime(iRX);          % duration of whole trial
    FTduration = original.trialtime(iFin) - original.trialtime(iVia);          % duration of follow-through
    
    if group<G.IMAGERYNOFIX
        ind = find(D.recordfbtime);
        duration(ind) = duration(ind) - 0.50;             % if an experiment where we record all the feedback period, remove this period from the 'moving' duration
        FTduration(ind) = FTduration(ind) - 0.50;
    end
    % find dwell time at via-point
    viapos = repmat([0;targetdist],[1,size(original.posx)]);
    viaradius = 1.25;
    posfromvia = (original.pos - viapos);
    distvia = sqrt(squeeze(posfromvia(1,:,:)).^2 + squeeze(posfromvia(2,:,:)).^2);
    ind1 = sub2ind(size(original.trialtime),[1:N]',findfirst((distvia<viaradius),2));           % subject enters via point
    ind2 = sub2ind(size(original.trialtime),[1:N]',findfirst((distvia<viaradius),2,1,'last'));  % subject leaves via point
    dwell = original.trialtime(ind2)-original.trialtime(ind1);
    peakspeed  = speed(lvmaxi);
    peakspeedch  = peakspeed(ph);  % channel trials
    peakspeedex  = peakspeed(ph1); % field trials
    
    % movement time to via and in follow-through
    movetimetovia = original.trialtime(iVia) - original.trialtime(iRX);
    movetimefollow = original.trialtime(iFin) - original.trialtime(iVia);
    
    % perpendicular error at halfway to via-point (6cm)
    halfwayy = 6;
    ind = sub2ind(size(F.posx),[1:N]',findfirst(F.posy>=halfwayy,2));
    halfwayperr = F.posx(ind).*fdir;
    
    timing.rxtime = rxtime;
    timing.appeartime = appeartime;
    timing.appeartimeaftergo = appeartimeaftergo;
    timing.peakspeed = peakspeed;
    timing.peakspeedch = peakspeedch;
    timing.peakspeedex = peakspeedex;
    timing.dwell = dwell;
    timing.duration = duration;
    timing.FTduration = FTduration;
    timing.halfwayperr = halfwayperr;
    timing.overshootrate = overshootrate;
    timing.movetimetovia = movetimetovia;
    timing.movetimefollow = movetimefollow;
    
    %------------------------
    %% Calculate motor imagery measures
    imagery = {};
    % read in table of MI measures
    if group<G.PLANNING
        imagery.MIQscore = D.MIscore;
        imagery.MIQvisual = D.MIQvisual;
        imagery.MIQkinesthetic = D.MIQkinesthetic;
        
        imagery.MImaintain = permute(reshape(D.maintenanceMIscore',[2,8,nsubj]),[3,2,1]);
        
        % button press timing for completion of imagined follow-throughs
        ind = find(field(:,S.FIELD_CHANNEL)~=1);
        imagineduration = FTduration(ind);
        timing.imagineduration = reshape(imagineduration, length(imagineduration)./nsubj,nsubj);
        
        % timing for completion of executed follow-throughs
        trialendwait = 0.1;
        ind = find(field(:,S.FIELD_CHANNEL)==1);
        executionduration = FTduration(ind)-trialendwait;                         % execution trials required 100ms wait at target before trial ended, so remove for a 'target hit' timing comparison.
        timing.executionduration = reshape(executionduration, length(executionduration)./nsubj,nsubj);
    end
    %------------------------
    F.original = original;
    %% save the data we care about according to group name
    MPE{group} = mpe;
    Adaptation{group} = adaptation;
    FrameData{group} = F;
    Timing{group} = timing;
    Experiment{group} = usefulthings;
    MI{group} = imagery;
    
    clearvars -except MPE FLAGS Adaptation FrameData Timing Experiment ngroups S MI AdaptationDetail G
end

%-----------------------------------------------------------
%% Plot our data across subjects and groups
% (and do a little data smoothing)
fh = [];
fh2 = [];
cfall = [];
speed_series = [];
speed_error = [];
imduration = [];
execduration = [];   
%colours = ColourSelect('ColourSelect.jpg',6); % selects colours from jpg
 colours(1,:) = [0.5453    0.6654    0.0451];  % manually fix colours
 colours(2,:) = [0.6673    0.2870    0.1166];
 colours(3,:) = [0.8843    0.5526    0.0197];
 colours(4,:) = [0.5072    0.4968    0.4941];
 colours(5,:) = [0.0876    0.4175    0.9734];
 colours(6,:) = [0.9462    0.1401    0.2563];

% choose a way to plot the individual subjects on the barplots
colours2 = colours.*0.6;   % make dots dark

% create plotting patch for exposure phase
P.Vertices = [ 4.5 -100; 155.5 -100; 155.5 100; 4.5 100];
P.Faces = [1 2 3 4];
P.FaceColor = [.3, .3, .3];
P.FaceAlpha = 0.08;
P.EdgeColor = 'white';
P.LineWidth = 0.1;

for group = 1:ngroups
    mpe = MPE{group};
    adaptation = Adaptation{group};
    F = FrameData{group};
    timing = Timing{group};
    imagery = MI{group};
    usefulthings = Experiment{group};
    S = usefulthings.S;
    N = usefulthings.N;
    D = usefulthings.D;
    field  = usefulthings.field;
    ntrial = usefulthings.ntrial;
    nsubj  = usefulthings.nsubj;
    fdir   = usefulthings.fdir;
    indB   = usefulthings.indB;
    indE   = usefulthings.indE;
    indP   = usefulthings.indP;
    %-----
    
    %% Plot MPE for each subject, and group average
    statsblocks = 6;

    pre=5;  exp=150; post=3;
    basephase = zeros(N,1); basephase(indB) = 1;
    exposurephase = zeros(N,1); exposurephase(indE) = 1;
    postphase = zeros(N,1); postphase(indP) = 1;
    indbase = find((field(:,S.FIELD_CHANNEL)~=1) & basephase);
    indexp = find((field(:,S.FIELD_CHANNEL)~=1) & exposurephase);
    indpost = find((field(:,S.FIELD_CHANNEL)~=1) & postphase);
    
    % smooth MPE by block (c=8)
    smooth_factor = 3;                                            % number of blocks to smooth across when plotting exposure phase
    ind = find(field(:,S.FIELD_CHANNEL)~=1);                      % null and exposure trials give mpe measures
    c = 8;                                                        % non-channel trials per block
    mpeblock = mean(reshape(mpe(ind), c, length(mpe(ind))/c),1);
    mpeblock = reshape(mpeblock,length(mpeblock)/nsubj,nsubj);
    if group~=G.FOLLOWTHROUGH                                     % if it's not the followthrough group, get rid of some extra data on the ends
        mpeblock = mpeblock(2:end-1,:);
    end
    
    mpebase = mean(reshape(mpe(indbase), c, length(mpe(indbase))/c),1);
    mpebase = reshape(mpebase,length(mpebase)/nsubj,nsubj);
  
    stats.mpebaseline{group} = mpebase;
    stats.firstmpeblock{group} = mpeblock(pre+1:pre+statsblocks,:);
    stats.finalmpeblock{group} = mpeblock(pre+exp+1-statsblocks:pre+exp,:);
    
    % smooth MPE by 3 blocks in exposure phase (c=24)
    c = c*smooth_factor;
    mpeblockexp = mean(reshape(mpe(indexp), c, length(mpe(indexp))/c),1);
    mpeblockexp = reshape(mpeblockexp,length(mpeblockexp)/nsubj,nsubj);
   
    % plot mpe per subject (smoothed per block only)
    if FLAGS.plotextra
        figure();
        for i=1:nsubj
            subplot(1,nsubj,i);
            patch(P); hold all;
            plot(linspace(1,ntrial,size(mpeblock,1)),mpeblock(:,i),'Color',colours(group,:)); hold on;
            plot([0 ntrial],[0,0],'k');
            axis([0 ntrial -1.5 4]);
            if i==1
                xlabel('Trials');
                ylabel('MPE (cm)');
            end
        end
        strng = sprintf('MPE per subject, group %d', group);
        suptitle(strng);
    end
    
    % average across subjects
    mpeblock_mean = mean(mpeblock,2);
    mpeblock_se = std(mpeblock,0,2) ./ sqrt(size(mpeblock,2));
    mpeblockexp_mean = mean(mpeblockexp,2);
    mpeblockexp_se = std(mpeblockexp,0,2) ./ sqrt(size(mpeblockexp,2));
    aftereffects{group} = mean(mpeblock(end-post+1:end-post+2,:),1);
    
    if group==G.PLANNING   % add in the transfer group data 
        figure(1000);
        x=4;
        subplot(1,2,1);
        load transferlearning
        shadeplot(1:pre, transfergrouplearning.mpeblock_mean(1:pre), transfergrouplearning.mpeblock_se(1:pre),'-',colours(6,:),0.3); hold all; % pre-exposure
        shadeplot(linspace(pre+1,exp+pre,length(transfergrouplearning.mpeblockexp_mean)), transfergrouplearning.mpeblockexp_mean, transfergrouplearning.mpeblockexp_se,'-',colours(6,:),0.3); hold all;  % exposure
        shadeplot(exp+pre+1:exp+pre+post, mpeblock_mean(end-post+1:end), mpeblock_se(end-post+1:end),'-',colours(6,:),0.3); hold all;  % pre-exposure
        plot([0 pre+exp+post],[0 0],'k');
        ylabel('MPE (cm)');
        xlabel('Block');
        axis([0 pre+exp+post+30 -1.5 4]);
        
        % plot the average after effects in a subpanel
        errorbar(pre+exp+post+3+x*4, transfergrouplearning.aftereffects_mean,transfergrouplearning.aftereffects_se,'k'); hold on;
        plot(pre+exp+post+3+x*4, transfergrouplearning.aftereffects_mean, 'o','MarkerSize',7,'MarkerEdgeColor',colours(6,:), 'MarkerFaceColor',colours(6,:)); hold on;
        
         % plot another figure of learning that separates out groups and only
        % includes the exposure phase
        figure(1100)
        subplot(2,ngroups+1,ngroups+1);
        shadeplot(linspace(1,exp,length(transfergrouplearning.mpeblockexp_mean)), transfergrouplearning.mpeblockexp_mean, transfergrouplearning.mpeblockexp_se,'-',[0.4 0.4 0.4],0.3); hold all;  % exposure
        shadeplot(linspace(exp+1,exp+transfergrouplearning.trans,length(transfergrouplearning.mpeblocktrans_mean)), transfergrouplearning.mpeblocktrans_mean, transfergrouplearning.mpeblocktrans_se,'-',[0.4 0.4 0.4],0.3); hold all;  % exposure
        plot([0 exp],[0 0],'k');
        ylabel('MPE (cm)');
        xlabel('Exposure Block');
        axis([0 exp -1 3.5])
    end

    % fix the order to plot the groups in, to the right of the learning plot
    switch group
        case G.NOIMAGERY
             x = 1;
        case G.IMAGERY
            x = 2;
        case G.IMAGERYNOFIX
            x = 3;
        case G.PLANNING
            x = 5;
        case G.FOLLOWTHROUGH
            x = 6;
    end
    
    % plot average across all subjects in group
    figure(1000);
    subplot(1,2,1);
    shadeplot(1:pre, mpeblock_mean(1:pre), mpeblock_se(1:pre),'-',colours(group,:),0.3); hold all; % pre-exposure
    shadeplot(linspace(pre+1,exp+pre,length(mpeblockexp_mean)), mpeblockexp_mean, mpeblockexp_se,'-',colours(group,:),0.3); hold all;  % exposure
    shadeplot(exp+pre+1:exp+pre+post, mpeblock_mean(end-post+1:end), mpeblock_se(end-post+1:end),'-',colours(group,:),0.3); hold all;  % pre-exposure
    plot([0 pre+exp+post],[0 0],'k');
    ylabel('MPE (cm)');
    xlabel('Block');
    axis([0 pre+exp+post+30 -1.5 4]);
    
    % plot the average after effects in a subpanel
    aftereffects_mean = mean(mean(mpeblock(end-post+1:end-post+2,:),1));
    aftereffects_se = std(mean(mpeblock(end-post+1:end-post+2,:),1)) ./ sqrt(nsubj);
    errorbar(pre+exp+post+3+x*4, aftereffects_mean,aftereffects_se,'k'); hold on;
    plot(pre+exp+post+3+x*4, aftereffects_mean, 'o','MarkerSize',7,'MarkerEdgeColor',colours(group,:), 'MarkerFaceColor',colours(group,:)); hold on;
    
    
    % plot another figure that separates out the learning by group
    figure(1100)
    subplot(2,ngroups+1,group);
    shadeplot(linspace(1,exp,length(mpeblockexp_mean)), mpeblockexp_mean, mpeblockexp_se,'-',colours(group,:),0.3); hold all;  % exposure
    plot([0 exp],[0 0],'k');
    ylabel('MPE (cm)');
    xlabel('Exposure Block');
    axis([0 exp -1 3.5])
    
    % plot the planning only group mean MPE response on all subplots as a
    % reference (exposure phase only)
    if group==G.PLANNING
       for i=1:ngroups+1
           subplot(2,ngroups+1,i)
           plot(linspace(1,exp,length(mpeblockexp_mean)), mpeblockexp_mean,'Color',colours(group,:)); hold all;
       end
    end
   
   % fix the desired order for each group in the barplot
    switch group
        case G.NOIMAGERY  
             x = 3;
        case G.IMAGERY
            x = 4;
        case G.IMAGERYNOFIX
            x = 5;
        case G.PLANNING
            x = 2;
        case G.FOLLOWTHROUGH
            x = 1;
    end
    
    % plot the after-effects as a barplot
    figure(1112)
    subplot(1,5,4:5);
    aftereffects_mean = mean(mean(mpeblock(end-post+1:end-post+2,:),1));
    aftereffects_se = std(mean(mpeblock(end-post+1:end-post+2,:),1)) ./ sqrt(nsubj);
    errorbar(x, -aftereffects_mean,aftereffects_se,'Color', colours(group,:),'CapSize',20, 'LineWidth',1); hold on;
    bar(x, -aftereffects_mean,'FaceColor',colours(group,:),'EdgeColor',colours(group,:)); hold on;
    for i=1:nsubj
       plot(x, -mean(mpeblock(end-post+1:end-post+2,i),1), 'o','MarkerSize',4,'MarkerEdgeColor',colours2(group,:),'MarkerFaceColor',colours2(group,:)); hold on; 
    end
    ylabel('After-effects (cm)');
    axis([0 7 -.3 1.8]);
    clear xticks; xticks([1:5])
    xticklabels({'Follow through','Planning only','No MI','MI','MI no fixation'});
    xtickangle(45)
    box off;
    
 
    %% Plot adaptation for each subject on the same figure (1 per group)
    pre=5;  exp=150; post=3;
    
    % smooth adaptation by block
    if group~=G.FOLLOWTHROUGH                                              % take channel trials from 0 and 180 degrees
        inclchannel = (round(D.HomeAngle)==-180) | (round(D.HomeAngle==0));
        ind = find((field(:,S.FIELD_CHANNEL)==1) & inclchannel);
        indexp = find((field(:,S.FIELD_CHANNEL)==1) & (exposurephase & inclchannel));
        basechannels =  find(((field(:,S.FIELD_CHANNEL)==1) & inclchannel) & basephase);
    else
        ind = find(field(:,S.FIELD_CHANNEL)==1);                            
        indexp = find((field(:,S.FIELD_CHANNEL)==1) & exposurephase);
        basechannels =  find(((field(:,S.FIELD_CHANNEL)==1)) & basephase);
    end
    
    incl0channel = (round(D.HomeAngle==0));
    basechannels0 =  find(((field(:,S.FIELD_CHANNEL)==1) & incl0channel) & basephase);
    ind0 = find((field(:,S.FIELD_CHANNEL)==1) & incl0channel);
    ind0exp = find((field(:,S.FIELD_CHANNEL)==1) & (exposurephase & incl0channel));
    
    c = 2;                                                                   % channel trials per block (for block = 12 trials, trials balanced across 24 trials)
    
     % find the baseline amount of adaptation
    adaptationbase = mean(reshape(adaptation(basechannels), c, length(adaptation(basechannels))/c),1);  % smooth by block
    adaptationbase = reshape(adaptationbase,length(adaptationbase)/nsubj,nsubj);
    adaptationbase0 = mean(reshape(adaptation(basechannels0), c, length(adaptation(basechannels0))/c),1);  % smooth by block
    adaptationbase0 = reshape(adaptationbase0,length(adaptationbase0)/nsubj,nsubj);
   
    % adaptation throughout training, baseline and post
    adaptationblock = mean(reshape(adaptation(ind), c, length(adaptation(ind))/c),1);  % smooth by block
    adaptationblock = reshape(adaptationblock,length(adaptationblock)/nsubj,nsubj);
    
    if group~=G.FOLLOWTHROUGH           % in the planning only and imagine conditions, get rid of some extra data on the ends so group plotting is comparable
        adaptationblock = adaptationblock(2:end-1,:);
        adaptation0block = reshape(adaptation(ind0),length(adaptation(ind0))/nsubj,nsubj);  % only 1 channel of this type per block
        adaptation0block = adaptation0block(2:end-1,:);
    else
        adaptation0block = mean(reshape(adaptation(ind0), c, length(adaptation(ind0))/c),1);  % smooth by block
        adaptation0block = reshape(adaptation0block,length(adaptation0block)/nsubj,nsubj);  % 2 channels of this type per block, and doesnt need end blocks clipped
    end
    
    % save data for statistical comparisons
    stats.adaptbaseline{group} = adaptationbase;
    stats.firstadaptblock{group} = adaptationblock(pre+1:pre+statsblocks,:);
    stats.finaladaptblock{group} = adaptationblock(pre+exp+1-statsblocks:pre+exp,:);
    
    % save data for statistical comparisons: check that the statistical conclusions are robust to inclusion of
    % only 0degree channel trials
    stats.adaptbaseline0{group} = adaptationbase0;
    stats.firstadapt0block{group} = adaptation0block(pre+1:pre+statsblocks,:);
    stats.finaladapt0block{group} = adaptation0block(pre+exp+1-statsblocks:pre+exp,:);
    stats.final4adaptblock{group} = adaptationblock(pre+exp+1-4:pre+exp,:);
    
    % smooth adaptation by 3 blocks in exposure phase 
    c = c*smooth_factor;
    adaptationblockexp = mean(reshape(adaptation(indexp), c, length(adaptation(indexp))/c),1);
    adaptationblockexp = reshape(adaptationblockexp,length(adaptationblockexp)/nsubj,nsubj);
    
    % plot by subject in group
    if FLAGS.plotextra
    figure();
    for i=1:nsubj
        subplot(1,nsubj,i);
        patch(P); hold all;
        plot(linspace(1,ntrial,size(adaptationblock,1)),adaptationblock(:,i),'Color',colours(group,:)); hold on;
        plot([0 ntrial],[0,0],'k');
        axis([0 ntrial -20 100]);
        if i==1
            xlabel('Trials');
            ylabel('Adapatation (%)');
        end
    end
    strng = sprintf('Adaptation per subject, group %d', group);
    suptitle(strng);
    end
    
    % save final adaptation measure for statistical comparisons 
    clear adaptationsubjfinalblock
    for i=1:nsubj
        adaptationsubjfinalblock(i) = mean(adaptationblock(pre+exp+1-statsblocks:pre+exp,i));
        adaptationsubjfinalblock_se(i) = std(adaptationblock(pre+exp+1-statsblocks:pre+exp,i))./sqrt(length(adaptationblock(pre+exp+1-statsblocks:pre+exp,i)));
    end
    AdaptationSubs{group} = adaptationsubjfinalblock;
    
    % Plot average ? se. adaptation across all subjects in group
    % average across subjects
    adaptblock_mean = mean(adaptationblock,2);
    adaptblock_se = std(adaptationblock,0,2) ./ sqrt(size(adaptationblock,2));
    adaptblockexp_mean = mean(adaptationblockexp,2);
    adaptblockexp_se = std(adaptationblockexp,0,2) ./ sqrt(size(adaptationblockexp,2));
    
    % plot the motor imagery transfer group adaptation response
    if group==G.PLANNING
        load transferlearning;
        x = 4;
        figure(1000);
        subplot(1,2,2);
        shadeplot(1:pre, transfergrouplearning.adaptblock_mean(1:pre), transfergrouplearning.adaptblock_se(1:pre),'-',colours(6,:),0.3); hold all; % pre-exposure
        shadeplot(linspace(pre+1,exp+pre,length(transfergrouplearning.adaptblockexp_mean)), transfergrouplearning.adaptblockexp_mean, transfergrouplearning.adaptblockexp_se,'-',colours(6,:),0.3); hold all;  % exposure
        h = shadeplot(exp+pre+1:exp+pre+post, transfergrouplearning.adaptblock_mean(end-post+1:end), transfergrouplearning.adaptblock_se(end-post+1:end),'-',colours(6,:),0.3); hold all;  % post-exposure
        plot([0 pre+exp+post],[0 0],'k');
        ylabel('Adaptation (%)');
        xlabel('Block');
        fh = [fh,h];
        
        % plot the average final adaptation level in a subpanel
        % plot transfer trials (FT) final adaptation
        imfinal_mean = mean(transfergrouplearning.adaptationsubjfinalblock);
        imfinal_se = std(transfergrouplearning.adaptationsubjfinalblock) ./ sqrt(nsubj);
        errorbar(pre+exp+post+3+x*4, imfinal_mean,imfinal_se,'k'); hold on;
        plot(pre+exp+post+3+x*4, imfinal_mean, 'o','MarkerSize',7,'MarkerEdgeColor',colours(6,:), 'MarkerFaceColor',colours(6,:)); hold on;
        
        % plot imagining trials final adaptation
        imfinal_mean = mean(transfergrouplearning.adaptationsubjfinalblockim);
        imfinal_se = std(transfergrouplearning.adaptationsubjfinalblockim) ./ sqrt(nsubj);
        errorbar(pre+exp+post+3+x*4, imfinal_mean,imfinal_se,'k'); hold on;
        plot(pre+exp+post+3+x*4, imfinal_mean, 'o','MarkerSize',7,'MarkerEdgeColor',colours(6,:), 'MarkerFaceColor','w'); hold on;
        axis([0 pre+exp+post+30 -20 60]);
        
        
         % plot another figure of learning that separates out groups and only
        % includes the exposure phase
        figure(1100)
        subplot(2,ngroups+1,(ngroups+1)*2);
        shadeplot(linspace(1,exp,length(transfergrouplearning.adaptblockexp_mean)), transfergrouplearning.adaptblockexp_mean, transfergrouplearning.adaptblockexp_se,'-',[0.4 0.4 0.4],0.3); hold all;  % exposure
        shadeplot(linspace(exp+1,exp+transfergrouplearning.trans,length(transfergrouplearning.adaptblocktrans_mean)), transfergrouplearning.adaptblocktrans_mean, transfergrouplearning.adaptblocktrans_se,'-',[0.4 0.4 0.4],0.3); hold all;  % exposure
        shadeplot(linspace(exp+1,exp+transfergrouplearning.trans,length(transfergrouplearning.adaptblocktransim_mean)), transfergrouplearning.adaptblocktransim_mean, transfergrouplearning.adaptblocktransim_se,'-',[0.1 0.9 0.9],0.3); hold all;  % exposure
        plot([0 exp],[0 0],'k');
        ylabel('Adaptation (%)');
        xlabel('Exposure Block');
        axis([0 exp -10 60])
        
        % plot imagining force expression
        figure(1112)
        subplot(1,5,3)
        bar(2, mean(transfergrouplearning.adaptationsubjfinalblockim),'EdgeColor',colours(6,:), 'FaceColor',colours(6,:)); hold on
        errorbar(2, mean(transfergrouplearning.adaptationsubjfinalblockim),std(transfergrouplearning.adaptationsubjfinalblockim)./sqrt(8), 'Color',colours(6,:),'CapSize',20, 'LineWidth',1); hold on
        
        % plot follow through force expression
        subplot(1,5,1:2)
        finaladapt_mean = mean(transfergrouplearning.adaptationsubjfinalblock);
        finaladapt_se = std(transfergrouplearning.adaptationsubjfinalblock) ./ sqrt(nsubj);
        errorbar(5, finaladapt_mean,finaladapt_se,'Color',colours(6,:),'CapSize',20, 'LineWidth',1); hold on;
        bar(5, finaladapt_mean,'FaceColor',colours(6,:),'EdgeColor',colours(6,:)); hold on;
        
        % plot each subject and lines connecting imagined and FT movement
        % force expression
        for i=1:nsubj
            subplot(1,5,3)
            plot(2, transfergrouplearning.adaptationsubjfinalblockim(i),'o','MarkerSize',4,'MarkerEdgeColor',colours2(6,:), 'MarkerFaceColor',colours2(6,:)); hold on
            subplot(1,5,1:2)
            plot(5, transfergrouplearning.adaptationsubjfinalblock(i),'o','MarkerSize',4,'MarkerEdgeColor',colours2(6,:), 'MarkerFaceColor',colours2(6,:)); hold on
        end
        box off;
    end
    
    % fix the group order in which the final adaptation measures are
    % plotted to the right of the learning figure
    switch group
        case G.NOIMAGERY 
            x = 1;
        case G.IMAGERY
            x = 2;
        case G.IMAGERYNOFIX
            x = 3;
        case G.PLANNING
            x = 5;
        case G.FOLLOWTHROUGH
            x = 6;
    end
    
    % plot the learning figure with all groups overlaid
    figure(1000);
    subplot(1,2,2);
    shadeplot(1:pre, adaptblock_mean(1:pre), adaptblock_se(1:pre),'-',colours(group,:),0.3); hold all; % pre-exposure
    shadeplot(linspace(pre+1,exp+pre,length(adaptblockexp_mean)), adaptblockexp_mean, adaptblockexp_se,'-',colours(group,:),0.3); hold all;  % exposure
    h = shadeplot(exp+pre+1:exp+pre+post, adaptblock_mean(end-post+1:end), adaptblock_se(end-post+1:end),'-',colours(group,:),0.3); hold all;  % post-exposure
    plot([0 pre+exp+post],[0 0],'k');
    ylabel('Adaptation (%)');
    xlabel('Block');
    fh = [fh,h];
    
    % plot the average final adaptation level in a subpanel
    finaladapt_mean = mean(adaptationsubjfinalblock);
    finaladapt_se = std(adaptationsubjfinalblock) ./ sqrt(nsubj);
    errorbar(pre+exp+post+3+x*4, finaladapt_mean,finaladapt_se,'k'); hold on;
    plot(pre+exp+post+3+x*4, finaladapt_mean, 'o','MarkerSize',7,'MarkerEdgeColor',colours(group,:), 'MarkerFaceColor',colours(group,:)); hold on;
    axis([0 pre+exp+post+30 -20 60]);
    
    % plot another figure of learning that separates out groups and only
    % includes the exposure phase
    figure(1100)
    subplot(2,ngroups+1,ngroups+1+group);
    shadeplot(linspace(1,exp,length(adaptblockexp_mean)), adaptblockexp_mean, adaptblockexp_se,'-',colours(group,:),0.3); hold all;  % exposure
    plot([0 exp],[0 0],'k');
    ylabel('Adaptation (%)');
    xlabel('Exposure Block');
    axis([0 exp -10 60])
    
    % repeat the planning only group mean adaptation response on all
    % subplots as a reference (exposure phase only)
    if group==G.PLANNING
        for i = 1:ngroups+1
            subplot(2,ngroups+1,ngroups+1+i);
            plot(linspace(1,exp,length(adaptblockexp_mean)), adaptblockexp_mean,'Color',colours(group,:)); hold all;
       end
    end
   
    
    % plot the average final adaptation level in a barplot
    switch group
        case G.NOIMAGERY 
             x = 3;
        case G.IMAGERY
            x = 4;
        case G.IMAGERYNOFIX
            x = 6;
        case G.PLANNING
            x = 2;
        case G.FOLLOWTHROUGH
            x = 1;
            poolFT = [AdaptationSubs{4}, AdaptationSubs{5}];  % pool data from FT and planning only groups
            
            % plot the pooled FT and planning only groups final adaptation
            % and compare to MI transfer channels
            figure(1112)
            subplot(1,5,3)  % this should come out dark blue, which is appropriate since coming a grey and a blue group
            bar(1, mean(poolFT),'EdgeColor',colours(5,:), 'FaceColor',colours(5,:)); hold on
            errorbar(1, mean(poolFT),std(poolFT)./sqrt(length(poolFT)), 'Color',colours(5,:),'CapSize',20, 'LineWidth',1); hold on
            
            for i=1:length(poolFT)
                plot(1, poolFT(i),'o','MarkerSize',4,'MarkerEdgeColor','k', 'MarkerFaceColor','k'); hold on
            end
            axis([0.2 2.7 -10 80]);
            clear xticks; xticks([1:2])
            xticklabels({'FT channels (pooled FT)','MI channels (MI transfer)'});
            xtickangle(45)
            ylabel('Adaptation (%)');
            box off;
    end
    figure(1112)
    subplot(1,5,1:2)
    finaladapt_mean = mean(adaptationsubjfinalblock);
    finaladapt_se = std(adaptationsubjfinalblock) ./ sqrt(nsubj);
    errorbar(x, finaladapt_mean,finaladapt_se,'Color',colours(group,:),'CapSize',20, 'LineWidth',1); hold on;
    bar(x, finaladapt_mean, 'EdgeColor',colours(group,:), 'FaceColor',colours(group,:)); hold on;
    for i=1:nsubj
           plot(x, adaptationsubjfinalblock(i),'o','MarkerSize',4,'MarkerEdgeColor',colours2(group,:), 'MarkerFaceColor',colours2(group,:)); hold on 
    end
    axis([0 7 -10 80])
    ylabel('Adaptation (%)'); 
    clear xticks; xticks([1:6])
    xticklabels({'Follow through','Planning only','No MI','MI','MI transfer','MI no fixation'});
    xtickangle(45)
    box off;
    
      
    %% Plot the trajectories for each group, at different experiment phases
    % Note: this takes a long time to execute
    if FLAGS.plotextra
        if group~=G.FOLLOWTHROUGH
            pre=[5,6];
            earlyexp=[1,2];
            lateexp=[149,150];
            post=[7,8];
        else
            pre=[4,5];
            earlyexp=[1,2];
            lateexp=[149,150];
            post=[6,7];
        end
        
        % plot the targets
        figure();
        colour_grey = [0.5, 0.5, 0.5];
        colour_lightgrey = [0.8, 0.8, 0.8];
        r = 1.25;
        F.rstartpos(5,:) = [0,-4]; % via position
        for k=1:4                  % phases
            for j=1:size(F.rstartpos,1)
                subplot(1,4,k);
                circle(F.rstartpos(j,1)+1000,F.rstartpos(j,2)+1000,r,'Color',colour_grey,'FillColor',colour_lightgrey); hold on;
                axis equal
            end
        end
        sz = size(F.rpos,2);
        % for each starting position, extract left and right trajectories for 2-block epochs in each phase
        homeangles = unique(round(D.HomeAngle));
        labels = {'pre','early','late','post'};
        for angle=1:length(homeangles)
            i = homeangles(angle);
            
            n1 = find(((field(:,S.FIELD_NULL)==1) & (fdir==1))  & (round(D.HomeAngle)==i)); %leftward target    0 degrees
            n2 = find(((field(:,S.FIELD_NULL)==1) & (fdir==-1))  & (round(D.HomeAngle)==i)); %rightward target
            n1 = reshape(n1,length(n1)/nsubj,nsubj);
            n2 = reshape(n2,length(n2)/nsubj,nsubj);
            trials.pre(:,1) = n1(pre,:);
            trials.post(:,1) = n1(post,:);
            trials.pre(:,2) = n2(pre,:);
            trials.post(:,2) = n2(post,:);
            
            e1 = find(((field(:,S.FIELD_VISCOUS)==1) & (fdir==1)) & (round(D.HomeAngle)==i)); %leftward target    0 degrees
            e2 = find(((field(:,S.FIELD_VISCOUS)==1) & (fdir==-1)) & (round(D.HomeAngle)==i)); %rightward target
            e1 = reshape(e1,length(e1)/nsubj,nsubj);
            e2 = reshape(e2,length(e2)/nsubj,nsubj);
            trials.early(:,1) = e1(earlyexp,:);
            trials.late(:,1) = e1(lateexp,:);
            trials.early(:,2) = e2(earlyexp,:);
            trials.late(:,2) = e2(lateexp,:);
            
            % average the trajectories across blocks and subjects
            for j=1:length(labels)
                for k=1:2    % number of different force fields
                    tmp = trials.(labels{j})(:,k);
                    x = mean(reshape(squeeze(F.rpos(1,:,tmp(:)))',[2,nsubj,sz]),1);
                    y = mean(reshape(squeeze(F.rpos(2,:,tmp(:)))',[2,nsubj,sz]),1);
                    trajectory_x.(labels{j})(k,:,:) = [squeeze(mean(x)),squeeze(std(x)./sqrt(nsubj))];
                    trajectory_y.(labels{j})(k,:,:) = [squeeze(mean(y)),squeeze(std(y)./sqrt(nsubj))];
                    
                    subplot(1,4,j);
                    if k==1
                        shadedTrajectory(squeeze(trajectory_x.(labels{j})(k,:,:))', squeeze(trajectory_y.(labels{j})(k,:,:))', [0.8 0.2 0.2], 0.1); hold on;
                    else
                        shadedTrajectory(squeeze(trajectory_x.(labels{j})(k,:,:))', squeeze(trajectory_y.(labels{j})(k,:,:))', [0.2 0.2 0.8], 0.1); hold on;
                    end
                    xlabel('Position (cm)');
                    ylabel('Position (cm)');
                    axis equal;
                    
                    if j==1
                        title('Pre-exposure');
                    elseif j==2
                        title('Early-exposure');
                    elseif j==3
                        title('Late-exposure');
                    else
                        title('Post-exposure');
                    end
                end
            end
            clear trajectory_x trajectory_y
        end
        strng = sprintf('Average hand paths, group %d ',group);
        suptitle(strng);
    end
    
    %% Plot timing of imagined and executed movements
    executecolour = [0.5,0.5,0.5];
    imaginecolour = [0,1,1];
    
    if group==G.PLANNING
        % lets just overwrite the planning only group's duration data here
        % with the transfer group because we dont care about planning only
        load transfergroup;
        v2struct(transferimagerydata);
        imagery.MIQscore = transferimagerydata.MIQscore;
        imagery.MImaintain = transferimagerydata.MImaintain;
        imagery.MIQkinesthetic = transferimagerydata.MIQkinesthetic;
        imagery.MIQvisual = transferimagerydata.MIQvisual;
    end
    
    % since we are temporarily using group=4 to mean the transfer group for now, compare all other groups for movement timing
    if group~=G.FOLLOWTHROUGH  
        % plot by subject in group
        if FLAGS.plotextra
            figure();
            for i=1:nsubj
                subplot(2,nsubj,i);
                h1 = plot(linspace(1,ntrial,size(timing.imagineduration,1)),timing.imagineduration(:,i),'Color',imaginecolour); hold on;
                h2 = plot(linspace(1,ntrial,size(timing.executionduration,1)),timing.executionduration(:,i),'Color',executecolour); hold on;
                plot([0 ntrial],[0,0],'k');
                axis([0 ntrial 0 2]);
                if i==1
                    xlabel('Trials');
                    ylabel('Trial Duration (s)');
                    legend([h1',h2'],{'Imagined','Executed'});
                end
                
                subplot(2,nsubj,nsubj+i);
                bar(1,mean(timing.imagineduration(:,i)),'FaceColor',imaginecolour); hold on;
                errorbar(1,mean(timing.imagineduration(:,i)),std(timing.imagineduration(:,i)),'Color',imaginecolour); hold on;
                bar(2,mean(timing.executionduration(:,i)),'FaceColor',executecolour); hold on;
                errorbar(2,mean(timing.executionduration(:,i)),std(timing.executionduration(:,i)),'Color',executecolour); hold on;
                ylabel('Movement duration (s)');
                ax=gca;
                ax.XTick = [1,2];
                ax.XTickLabel = {'Imagined','Executed'};
                ax.XTickLabelRotation = 50;
            end
            suptitle('Movement duration per subject');
        end
        
        % average for each subject
        averageimagineduration = mean(timing.imagineduration,1);
        averageexecutionduration = mean(timing.executionduration,1);
        fieldnm = {'','fixation','nofixation','transferimagery'};
        if ((group==G.IMAGERY) ||(group==G.IMAGERYNOFIX)) || (group==G.PLANNING)
            imduration = [imduration,averageimagineduration];
            execduration = [execduration,averageexecutionduration];
            
            % plot average mental chronometry across all subjects in group
            figure(2000);
            cf = plot(averageexecutionduration,averageimagineduration,'.', 'MarkerSize',15, 'Color', colours(group,:)); hold on;
            cfall = [cfall,cf];
            xlabel('Execution duration (s)');
            ylabel('Imagining duration (s)');
            axis([0.1 1.7 .1 1.7])
            
            imaginedurationgroups(group-1,:,:) = averageimagineduration;
            executiondurationgroups(group-1,:,:) = averageexecutionduration;
            if group==G.PLANNING  % on final loop of the imagery groups, plot collapsed results
                
                % plot correlation between imagined and executed follow
                % through durations
                y = reshape(imduration,1,size(imduration,1)*size(imduration,2))';   % imagining durations
                exdata = reshape(execduration,1,size(execduration,1)*size(execduration,2))';  % execution durations
                X = [exdata];                       % slope, no intercept
                [b,bint,~,~,imagestats]=regress(y,X);                % b = [slope]
                [r,pim] = corr(y,X);  % significant correlation, p=0.0000...
                
                x = linspace(min(execduration), max(execduration),10);
                cf = plot(x, b.*x,'b');
                legend([cfall(1)',cfall(2)',cfall(3)',cf'],{'Motor imagery with fixation','Motor imagery','Motor imagery transfer','Regression fit'});
                strng = sprintf('R = %.2f, p = %.7f', r, pim);
                text(0.3,1,strng);
                axis([0 1.5 0 1.5]);
                axis equal
            end
        end
    end
    
    %% Plot correlations between final adaptation level and MI scores
    if ((group==G.IMAGERY) ||(group==G.IMAGERYNOFIX)) || (group==G.PLANNING)  % again using planning group as placeholder for imagery transfer data
        figure(3000);
        subplot(1,3,1);  % plot final adaptation level vs. MIQ-score
        h1=plot(imagery.MIQscore,adaptationsubjfinalblock,'xk'); hold on;
        h2=plot(imagery.MIQkinesthetic, adaptationsubjfinalblock, '.b'); hold on;
        h3=plot(imagery.MIQvisual, adaptationsubjfinalblock, '.c'); hold on;
        
        XMIQ(:,group-1) = imagery.MIQscore;
        xlabel('Mean MIQ-score');
        ylabel('Adaptation level (%)');
        legend([h1',h2',h3'], {'Both MIQ','MIQ kinesthetic','MIQ visual'});
        
        subplot(1,3,2);  % plot final adaptation level vs. MIQ-maintenance frequency.
        plot(squeeze(mean(imagery.MImaintain(:,:,2),2)),adaptationsubjfinalblock,'.','MarkerSize',15, 'Color', colours(group,:)); hold on;
        for i=1:nsubj
            errorbar(squeeze(mean(imagery.MImaintain(i,:,2),2)),adaptationsubjfinalblock(i), adaptationsubjfinalblock_se(i),'Color', colours(group,:)); hold on;
        end
        
        tmp = mean(imagery.MImaintain(:,:,2),2);
        Xmaintain(:,group-1) = tmp;
        xlabel('Mean MI-maintenance frequency');
        ylabel('Adaptation level (%)');
        
        subplot(1,3,3);  % plot final adaptation level vs. MIQ-maintenance ease.
        plot(squeeze(mean(imagery.MImaintain(:,:,1),2)),adaptationsubjfinalblock,'.','MarkerSize',15, 'Color', colours(group,:)); hold on;
        for i=1:nsubj
            errorbar(squeeze(mean(imagery.MImaintain(i,:,1),2)),adaptationsubjfinalblock(i), adaptationsubjfinalblock_se(i),'Color', colours(group,:)); hold on;
        end
        
        tmp = mean(imagery.MImaintain(:,:,1),2);
        Xmaintainease(:,group-1) = tmp;
        xlabel('Mean MI-maintenance ease');
        ylabel('Adaptation level (%)');
        
        finaladapt(:,group-1) = adaptationsubjfinalblock';
        
        % regress and plot correlation on final loop through imagery data
        if group==G.PLANNING
            XMIQtotal = reshape(XMIQ, numel(XMIQ),1);
            Xmaintaintotal = reshape(Xmaintain, numel(Xmaintain),1);
            Xmaintaineasetotal = reshape(Xmaintainease, numel(Xmaintainease),1);
            y = reshape(finaladapt, numel(finaladapt),1);
            
            subplot(1,3,1);
            [b,bint,~,~,imagestats]=regress(y,[XMIQtotal, ones(size(XMIQtotal))]);               
            [r,pim] = corrcoef(y,XMIQtotal);
            plot(0:10,b(1).*[0:10]+b(2), '--b'); hold on;
            strng = sprintf('R = %.2f, p = %.4f', r(1,2),pim(1,2));
            text(4,20,strng);
            
            subplot(1,3,2);
            [b,bint,~,~,imagestats]=regress(y,[Xmaintaintotal, ones(size(Xmaintaintotal))]);                
            [r,pim] = corrcoef(y,Xmaintaintotal);
            plot(0:3,b(1).*[0:3]+b(2), '--b'); hold on;
            strng = sprintf('R = %.2f, p = %.4f', r(1,2), pim(1,2));
            text(1,20,strng);
            
            subplot(1,3,3);
            [b,bint,~,~,imagestats]=regress(y,[Xmaintaineasetotal, ones(size(Xmaintaineasetotal))]);                
            [r,pim] = corrcoef(y,Xmaintaineasetotal);
            plot(0:10,b(1).*[0:10]+b(2), '--b'); hold on;
            strng = sprintf('R = %.2f, p = %.4f', r(1,2), pim(1,2));
            text(1,20,strng);
                            
            % restore the planning only group's timing data
            timing = Timing{group};
        end
    end
   
    clearvars -except finaladapt  AdaptationSubs transferimagerydata colours2 G  Xmaintainease Xmaintain XMIQ Adaptation AdaptationDetail P MPE stats ngroups nsubj ntrial N S FLAGS Timing Experiment colours FrameData aftereffects  fh fh2 A M MI speed_series speed_error imduration execduration cfall chronometry imaginedurationgroups executiondurationgroups
end
%%
figure(1000);                       % plot some background shaded areas
subplot(1,2,1)
pre=5;  exp=150; post=3;
Pg = P;
Pg.Vertices = [ pre+.5 -100; exp+pre+.5 -100; exp+pre+.5 100; pre+.5 100];
patch(Pg); hold on;
subplot(1,2,2)
patch(P); hold on;
legend([fh(1)',fh(2)',fh(3)',fh(4)',fh(5)',fh(6)'],{'No imagery','Imagery','Imagery no fixation','Imagery transfer','Planning only','Follow-through'});

figure(1200);
subplot(1,2,1);
legend([fh2(1)',fh2(2)'],{'Imagery (fixation)','Imagery (free-eye movement)'});
patch(Pg); hold on;
subplot(1,2,2);
patch(P); hold on;

clear pre exp post fh fh2 P Pg ntrial N

%------------------------------
%% For the imagining groups, plot the correlation between imagined and executed movements
y = reshape(imduration,1,size(imduration,1)*size(imduration,2))';   % imagining durations
exdata = reshape(execduration,1,size(execduration,1)*size(execduration,2))';  % execution durations
X = [exdata];                       % slope, no intercept
[b,bint,~,~,imagestats]=regress(y,X);                % b = [slope]
[r,pim] = corr(y,X);  % significant correlation, p=0.0000...

figure(2000);
subplot(1,2,1);
x = linspace(min(execduration), max(execduration),10);
cf = plot(x, b.*x,'b');
legend([cfall(1)',cfall(2)',cfall(3)',cf'],{'Motor imagery with fixation','Motor imagery','Motor imagery transfer','Regression fit'});
strng = sprintf('R = %.2f', r);
text(1.7,1.7,strng);
axis([0 2 0 2]);
axis equal

%------------------------------
%% Statistical analysis - learning within each group
% Used epochs of 6 blocks (= 48 exposure trials, 12 channel trials), and
% compare pre-exposure to late-exposure performance.

% load learning stats data from transfer group
v2struct(transferimagerydata);
load transferlearning;

% Assess within group learning using adaptation, MPE and after-effects
groupname = {'noimagery','imageryfix','imagerynofix','planning','fullfollow'};
for group=1:ngroups        
    [~, statstable.mpe{group}] = anova_rm([squeeze(mean(stats.firstmpeblock{group}))', squeeze(mean(stats.finalmpeblock{group}))'],'off');
   
    % Compare MPE and adaptation pre- and late-exposure phase with paired ttests
    [~, p.adaptation{group},~,statstable.adaptationtable{group}] = ttest(squeeze(mean(stats.adaptbaseline{group}))', squeeze(mean(stats.finaladaptblock{group}))');
    [~, p.adaptation0{group},~,statstable.adaptation0onlytable{group}] = ttest(squeeze(mean(stats.adaptbaseline0{group}))', squeeze(mean(stats.finaladapt0block{group}))');
    adaptstats.(groupname{group}) = [squeeze(mean(stats.adaptbaseline{group}))', squeeze(mean(stats.finaladaptblock{group}))']; 
    adaptstats0only.(groupname{group}) = [squeeze(mean(stats.adaptbaseline0{group}))', squeeze(mean(stats.finaladapt0block{group}))']; 
    mpestats.(groupname{group}) = [squeeze(mean(stats.firstmpeblock{group}))', squeeze(mean(stats.finalmpeblock{group}))'];
    
    % Assess after-effect magnitude (paired ttest compared to pre-exposure)
    [~,tmp_p,~,tmp_s] = ttest(aftereffects{group}',mean(stats.mpebaseline{group})');
    p.aftereffects.(groupname{group}) = tmp_p;
    statstable.aftereffects.(groupname{group}) = tmp_s;
    
    % Print mean difference pre- to late-exposure for each learning measure
    mn = mean(mean(stats.finaladaptblock{group})- mean(stats.adaptbaseline{group}));
    se = std(mean(stats.finaladaptblock{group})- mean(stats.adaptbaseline{group})) ./ sqrt(8);
    %sprintf('Diff adaptation %s: %.2f  +/- %.2f', groupname{group}, mn, se)
    
    mn = mean(aftereffects{group}'- mean(stats.mpebaseline{group})');
    se = std(aftereffects{group}'- mean(stats.mpebaseline{group})') ./ sqrt(8);
    %sprintf('Diff adaptation %s: %.2f  +/- %.4f', groupname{group}, mn, se)    
end

% compare adaptation within motor imagery transfer group: adaptation on FT vs imagery channels
[H,P,CI,STATS] = ttest2(transfergrouplearning.adaptationsubjfinalblock',transfergrouplearning.adaptationsubjfinalblockim');

%% Statistical analysis - adaptation between groups
% Perform planned comparisons between groups using unpaired t-tests.

% pooled group follow through and planning only
poolft = [mean(stats.finaladaptblock{4}), mean(stats.finaladaptblock{5})];

% pooled fixation imagery groups
poolimfix = [mean(stats.finaladaptblock{2}), mean(transfergrouplearning.finaladaptblock{1})];

% compare adaptation: pooled FT vs pooled imagery fixation groups
[H,P,CI,STATS] = ttest2(poolft, poolimfix)

% compare adaptation: pooled imagery fixation vs no imagery group
[H,P,CI,STATS] = ttest2(poolimfix, mean(stats.finaladaptblock{1}))

% compare adaptation: pooled imagery fixation vs imagery group
[H,P,CI,STATS] = ttest2(poolimfix, mean(stats.finaladaptblock{3}))

% compare adaptation: pooled FT vs imagery channel trials (transfer group)
[H,P,CI,STATS] = ttest2(poolft, transfergrouplearning.adaptationsubjfinalblockim)


%% Statistical analysis - compare other behavioural factors

% Compare mistrial percentages for hand overshoots of the central target
% compare overshoots: pooled imagery fixation group vs no MI group
poolimfix_handovershoot = [transferimagerydata.percerrhandovershoot, Timing{2}.percerrhandovershoot];
[H,P,CI,STATS] = ttest2(poolimfix_handovershoot, Timing{1}.percerrhandovershoot)

% Compare mistrial percentages for breaks in fixation (excl. blinks or breaks before 'go' cue to move)
% compare fixation breaks: pooled imagery fixation group vs no MI group
poolimfix_eyefixation = [transferimagerydata.percerrsaccadetotarget, Timing{2}.percerreyebreakfixation];
[H,P,CI,STATS] = ttest2(poolimfix_eyefixation, Timing{1}.percerreyebreakfixation)

% Compare wait time at the central target on imagery trials
% time at central target: no imagery vs imagery fixation group
[H,P,CI,STATS] = ttest2(mean(Timing{1}.imagineduration), mean(Timing{2}.imagineduration))


%% Statistical analysis - Baseline kinematics by Left vs Right followthrough
% check for baseline difference for L vs R follow targets for 4
% kinematic variables in each group. Note that this give 10 datapoints per
% Left or Right condition per subject in groups1-4, and 8 datapoints per L or R per sub in group 5.
% (use mean value per subject to contribute to group stats & consider just the 0 degree starting position)
for i=1:ngroups
   FIELD_NULL = Experiment{i}.S.FIELD_NULL; % same for all groups
   ind = (Experiment{i}.field(:,FIELD_NULL) & Experiment{i}.D.HomeAngle==0) & (Experiment{i}.D.HomePosition(:,2)==[-16]); % hack for choosing 0degree start position. Seems redundant?
   indL = ind & Experiment{i}.D.TargetAngle==45;
   indR = ind & Experiment{i}.D.TargetAngle==-45;
   
   % baseline peak speed by FT target
   baseLeft.(groupname{i}) = mean(reshape(Timing{i}.peakspeed(indL), length(Timing{i}.peakspeed(indL))./nsubj, nsubj),1)';
   baseRight.(groupname{i}) = mean(reshape(Timing{i}.peakspeed(indR), length(Timing{i}.peakspeed(indR))./nsubj, nsubj),1)'; 
   [~, statstable.peakspeedtarget{i}] = anova_rm([baseLeft.(groupname{i}),baseRight.(groupname{i})],'off');    
   difference.peakspeedtarget(i) = mean(baseLeft.(groupname{i})) - mean(baseRight.(groupname{i}));
  
   % baseline duration by FT target
   baseLeft.(groupname{i}) = mean(reshape(Timing{i}.duration(indL), length(Timing{i}.duration(indL))./nsubj, nsubj),1)';
   baseRight.(groupname{i}) = mean(reshape(Timing{i}.duration(indR), length(Timing{i}.duration(indR))./nsubj, nsubj),1)'; 
   [~, statstable.durationtarget{i}] = anova_rm([baseLeft.(groupname{i}),baseRight.(groupname{i})],'off');    
   difference.durationtarget(i) = mean(baseLeft.(groupname{i})) - mean(baseRight.(groupname{i}));
   
   % baseline path length by FT target
   baseLeft.(groupname{i}) = mean(reshape(FrameData{i}.original.pathlength(indL), length(FrameData{i}.original.pathlength(indL))./nsubj, nsubj),1)';
   baseRight.(groupname{i}) = mean(reshape(FrameData{i}.original.pathlength(indR), length(FrameData{i}.original.pathlength(indR))./nsubj, nsubj),1)'; 
   [~, statstable.pathlengthtarget{i}] = anova_rm([baseLeft.(groupname{i}),baseRight.(groupname{i})],'off');    
   difference.pathlengthtarget(i) = mean(baseLeft.(groupname{i})) - mean(baseRight.(groupname{i}));
   
   % baseline dwell time by FT target (Note that this measure only makes sense for the full follow-through group)
   baseLeft.(groupname{i}) = mean(reshape(Timing{i}.dwell(indL), length(Timing{i}.dwell(indL))./nsubj, nsubj),1)';
   baseRight.(groupname{i}) = mean(reshape(Timing{i}.dwell(indR), length(Timing{i}.dwell(indR))./nsubj, nsubj),1)'; 
   [~, statstable.dwelltarget{i}] = anova_rm([baseLeft.(groupname{i}),baseRight.(groupname{i})],'off');    
   difference.dwelltarget(i) = mean(baseLeft.(groupname{i})) - mean(baseRight.(groupname{i}));

   % baseline maximum lateral deviation by FT target
   baseLeft.(groupname{i}) = mean(reshape(FrameData{i}.lateraldev(indL), length(FrameData{i}.lateraldev(indL))./nsubj, nsubj),1)';
   baseRight.(groupname{i}) = mean(reshape(FrameData{i}.lateraldev(indR), length(FrameData{i}.lateraldev(indR))./nsubj, nsubj),1)'; 
   [~, statstable.lateraldevtarget{i}] = anova_rm([baseLeft.(groupname{i}),baseRight.(groupname{i})],'off');    
   difference.lateraldevtarget(i) = mean(baseLeft.(groupname{i})) - mean(baseRight.(groupname{i}));
end

%%

toc
