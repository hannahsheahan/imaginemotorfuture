%% Motor Planning, Not Execution, Separates Motor Memories
% Sheahan, Ingram, Zalalyte and Wolpert (2018)
% - This is an analysis script for assessing learning in a motor imagery transfer group (n=8).
% To assess how much learning from imagining follow-throughs transfers to
% when planning and executing them.

% Group 1: imagery transfer (fixation)
% Group 2: no imagery transfer (fixation)


% Author: Hannah Sheahan, sheahan.hannah@gmail.com
% Date:   Feb 2018
%--------------------------------------------------------------------------

%% settings
clc; clear all; close all;
tic
FontSize = 14;
set(0, 'DefaultAxesFontSize',FontSize); clear FontSize;
set(0,'DefaultFigureWindowStyle','docked');
FLAGS.normMPE = 0;              % normalise MPE by group peak speed
FLAGS.plotextra = 0;            % plot additional figures (subject-by-subject analysis, hand paths, target appearance time distributions)

%% load subject data
ngroups = 2;
for group = 1:ngroups  
    switch group  % load in pre-analysed data so you don't have to do this again
        case 1
            datname = 'NoImageryFixationTransferExl1_analysed.mat';
            if ~exist(datname,'file')
                load NoImageryFixationTransferExl1_combined;
            end
        case 2
            datname = 'ImagineFixationTransfer_analysed.mat';
            if ~exist(datname,'file')
                load ImagineFixationTransfer_combined;
            end
    end
    
    % load the pre-analysed data if it exists so we dont have to reanalyse it every time
    if exist(datname,'file')
        load(datname);
    else
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
        
        % State Names
        S.STATE_GO       = 5;
        S.STATE_MOVEWAIT = 6;
        S.STATE_MOVING0  = 7;
        S.STATE_MOVING1  = 8;
        S.STATE_FINISH   = 9;
        
        % binary field-type matrix
        fieldlist = unique(D.FieldType);
        field = zeros(N,length(fieldlist)+1);
        for i=1:length(fieldlist)
            FN = fieldlist(i)+1;
            field(:,FN) = D.FieldType==fieldlist(i);
        end
        
        % Experiment phases
        phaselist = cumsum([3 75 5 2]);
        
        Baseline = D.PhaseIndex<=phaselist(1);  indB = find(Baseline==1);
        exposurephase = D.PhaseIndex>phaselist(1) & D.PhaseIndex<=phaselist(2);  indE = find(exposurephase==1);
        transferphase = D.PhaseIndex>phaselist(2) & D.PhaseIndex<=phaselist(3);  indT = find(transferphase==1);
        Post     = D.PhaseIndex>phaselist(3) & D.PhaseIndex<=phaselist(4);  indP = find(Post==1);
        
        clear i FN fieldlist phaselist Baseline Exposure Post
        HomeAngle    = D.HomeAngle;
        HomePosition = D.HomePosition;
        TargetAngle  = D.TargetAngle;
        ContextType  = D.ContextType;
        usefulthings  = v2struct(S,N,field,ntrial,nsubj,indB,indE,indT,indP,HomeAngle,HomePosition,ContextType,TargetAngle);    % save some general access variables             % save some general access variables
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
        rvel = [reshape(velx',[1,nsamp,N]);reshape(vely',[1,nsamp,N])];
        rforce = [reshape(forcex',[1,nsamp,N]);reshape(forcey',[1,nsamp,N])];
        original.pos = [reshape(original.posx',[1,original.lngth,N]);reshape(original.posy',[1,original.lngth,N])];
        original.pathlength = pathlength;
        
        clear  posx posy velx vely forcex forcey ind1 ind2 state trialtime pathlength
        
        %% Rotate kinematics from each start position
        startangle = -HomeAngle.*(pi/180);
        startpos   = repmat(reshape(D.StartPosition(:,1:2)', [2,1,N]),[1,nsamp,1]);
        original.startpos = repmat(reshape(D.StartPosition(:,1:2)', [2,1,N]),[1,original.lngth,1]);
        
        rotate = @(x,theta,N) cell2mat(( arrayfun(@(k)(reshape([cos(theta(k)), -sin(theta(k)); sin(theta(k)), cos(theta(k))]' * squeeze(x(:,:,k)),[2,1,size(x,2)])), 1:N, 'UniformOutput', false)));
        
        original.pos = rotate((original.pos-original.startpos),startangle,N);
        pos   = rotate((F.rpos-startpos),startangle,N);
        vel   = rotate(rvel,startangle,N);
        force = rotate(rforce,startangle,N);
        
        % rotate misstrial eye and hand data (not resampled)
        misstrialdata = D.MissTrialData;
        misstrialstartangle = -misstrialdata.HomeAngle.*(pi/180);
        viapos = misstrialdata.ViaPosition(:,[1,2]);
        viapos = permute(repmat(viapos,[1,1,size(misstrialdata.FrameData.EyeTrackerEyeXY,2)]),[1,3,2]);
        eyepos = rotate(permute((misstrialdata.FrameData.EyeTrackerEyeXY-viapos),[3,2,1]), misstrialstartangle, misstrialdata.Trials);
        handpos = rotate(permute((misstrialdata.FrameData.RobotPosition(:,:,[1:2])-viapos),[3,2,1]), misstrialstartangle, misstrialdata.Trials);
        F.mistrialeyepos = eyepos;
        F.mistrialhandpos = handpos;
        
        posx = squeeze(pos(1,:,:));
        posy = squeeze(pos(2,:,:));
        velx = squeeze(vel(1,:,:));
        vely = squeeze(vel(2,:,:));
        forcex = squeeze(force(1,:,:));
        forcey = squeeze(force(2,:,:));
        F.rstartpos = unique(squeeze(startpos(:,1,:))','rows');
        original.posx = squeeze(original.pos(1,:,:));
        original.posy = squeeze(original.pos(2,:,:));
        
        %% Calculate MPE and adaptation
        
        % calculate mpe
        [test,ind] = max(abs(posx'));
        ind = sub2ind(size(posx),1:N,ind);
        mpe = posx(ind)'.*(fdir==1) - posx(ind)'.*(fdir==-1);     % accounts for counterbalance subjects and switching field directions
        
        % calculate maximum signed lateral deviation
        F.lateraldev = posx(ind)';
        
        % Calculate adaptation (% complete)
        ph1 = find(field(:,S.FIELD_VISCOUS)==1);
        fieldconstant = D.FieldConstants(ph1(1),1);                   % same for all subjects
        perffx = fieldconstant.*repmat(fdir,[1,1000]).*vely;        % ideal force. NB: we can only measure forces in the x-direction
        
        % find +-200ms window either side of time of max speed (NB: often takes entire length of trial, but limits very long trials)
        speed = sqrt((velx.^2)+(vely.^2));
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
        adaptation = adaptregress(perffx,forcex,p1,p2,N);
        
        %------------------------
        %% Calculate misstrials
        
        % Arm movement misstrials (for planning only group mainly):
        % percentage of misstrials per subject
        for i=1:nsubj
            ind = [ntrial.*(i-1)+1,ntrial.*i];
            misstrialrate(i) = max(D.MissTrials(ind(1):ind(2))) ./ (max(D.MissTrials(ind(1):ind(2)))+ntrial);
        end
        
        % percentage of overshot trials per subject (for planning only group)
        %(NB: overshot = if subject travels 3cm+ (missthresh) in the rotated +y direction after the viapoint)
        missthresh = 3;
        targetdist = 12;
        overdistance = original.posy - repmat(targetdist,size(original.posy));
        
        ch = (original.state >= S.STATE_MOVING1);                               % after hitting the via-point
        overshot = sign(sum(ch & (overdistance >= missthresh),2));
        overshot(field(:,S.FIELD_VISCOUS)~=1)=0;                                % ignore count on channel trials (subjects supposed to follow through) and null trials
        nincludedtrials = sum(field(1:ntrial,S.FIELD_VISCOUS)==1 & exposurephase(1:ntrial));  % consider only main exposure phase
        for i=1:nsubj
            ind = [ntrial.*(i-1)+1,ntrial.*i];
            overshootrate(i) = sum(overshot(ind(1):ind(2))) ./ nincludedtrials;
            novershot(i) =  sum(overshot(ind(1):ind(2)));
        end
        
        % Fixation misstrials (for eyetracked groups only)
        tolerance = 3;   % cm
        nimagine = sum(((field(:,S.FIELD_VISCOUS)==1) & D.subj==1) & exposurephase);  % same for each subject, so use subj==1 as example
        
        if group<3
            ind = (misstrialdata.FieldType==1) & (misstrialdata.PhaseIndex<78);       % imagining trials in exposure phase (do not include transfer phase)
            
            % measure %overshoot relative to all attempts at imagining trials
            for i=1:nsubj
                nimaginesub(i) = nimagine + sum(ind & misstrialdata.subj==i);
            end
            
            eyehittarget = findfirst((sqrt(eyepos(2,:,:).^2 + eyepos(1,:,:).^2) > tolerance),3);  % overshoot of via target in rotated y direction
            handhittarget = findfirst(handpos(2,:,:) > tolerance,3);
            
            % remove trials that were blinks or where eyetracker loses eye
            blink = squeeze(findfirst(eyepos(1,:,:)>100 | eyepos(2,:,:)>100,3));
            eyehittarget(blink>0) = 0;
            eyehittarget(eyehittarget<300) = 0;  % don't count instances where fixation is broken before 'go' cue.
            
            % per subject misstrial errors on imagining trials
            for i=1:nsubj
                trials = find((misstrialdata.subj==i) & ind);
                nerrsaccadeovershoot(i) = sum(eyehittarget(trials)~=0);                     % number of times eye moved to target on imagining trials
                nerrhandovershoot(i) = sum(handhittarget(trials)~=0);                       % number of times hand reached to target on imagining trials
            end
            
            % percentage of imagining trials that eye-moves or hand-moves
            % overshot central target (consider both regular overshoots and
            % overshoots on saccade error trials).
            timing.percerrsaccadetotarget = nerrsaccadeovershoot./nimaginesub;
            timing.nerrsaccadeovershoot = nerrsaccadeovershoot;
            F.misstrialdata = misstrialdata;
            
        else
            % a few subjects here have misstrialdata saved separately, so make
            % sure this is included for hand overshoot data
            if isfield(D,'MissTrialData')
                
                misstrialdata = D.MissTrialData;
                ind = misstrialdata.FieldType==1;   % imagining trials in exposure phase
                
                % measure %overshoot relative to all attempts at imagining trials
                for i=1:nsubj
                    nimaginesub(i) = nimagine + sum(ind & misstrialdata.subj==i);
                end
                
                viapos = misstrialdata.ViaPosition(:,[1,2]);
                viapos = permute(repmat(viapos,[1,1,size(misstrialdata.FrameData.RobotPosition,2)]),[1,3,2]);
                handpos = rotate(permute((misstrialdata.FrameData.RobotPosition(:,:,[1:2])-viapos),[3,2,1]),startangle,misstrialdata.Trials);
                
                % hand overshot central target by 3cm+ in y direction
                handhittarget = findfirst(handpos > tolerance,2);
                
                % per subject misstrial errors on imagining/aborting trials
                for i=1:nsubj
                    trials = find((misstrialdata.subj==i) & ind);
                    nerrhandovershoot(i) = sum(handhittarget(trials)~=0);                       % number of times hand reached to target on imagining trials
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
        
        if group<3
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
        ind = sub2ind(size(posx),[1:N]',findfirst(posy>=halfwayy,2));
        halfwayperr = posx(ind).*fdir;
        
        timing.peakspeed = peakspeed;
        timing.dwell = dwell;
        timing.duration = duration;
        timing.FTduration = FTduration;
        timing.halfwayperr = halfwayperr;
        timing.overshootrate = overshootrate;
        %timing.peakspeedch = peakspeedch;
        %timing.peakspeedex = peakspeedex;
        %timing.rxtime = rxtime;
        %timing.appeartime = appeartime;
        %timing.appeartimeaftergo = appeartimeaftergo;
        %timing.movetimetovia = movetimetovia;
        %timing.movetimefollow = movetimefollow;
        %------------------------
        %% Calculate motor imagery measures
        imagery = {};
        % read in table of MI measures
        if group<4
            imagery.MIQscore = D.MIscore;
            imagery.MIQvisual = D.MIQvisual;
            imagery.MIQkinesthetic = D.MIQkinesthetic;
            
            imagery.MImaintain = permute(reshape(D.maintenanceMIscore',[2,9,nsubj]),[3,2,1]);
            
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
        %F.original = original;
        F.original.pathlength = original.pathlength;
        
        %% save the data to a pre-analysed file so we don't have to do this again
        save(datname, 'mpe', 'adaptation', 'F', 'timing', 'usefulthings', 'imagery','-v7.3');
    end
    %% save the data we care about according to group name
    MPE{group} = mpe;
    Adaptation{group} = adaptation;
    FrameData{group} = F;
    Timing{group} = timing;
    Experiment{group} = usefulthings;
    MI{group} = imagery;
    
    clearvars -except MPE FLAGS Adaptation FrameData Timing Experiment ngroups S MI
end

%-----------------------------------------------------------
%% Plot our data across subjects and groups
% (and do a little data smoothing)
fh = [];
fh2 = [];
speed_series = [];
speed_error = [];
imduration = [];
execduration = [];
colours = ColourSelect('ColourSelect.jpg',ngroups+2);
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
    %    adaptationdetails = AdaptationDetail{group};
    F = FrameData{group};
    timing = Timing{group};
    imagery = MI{group};
    usefulthings = Experiment{group};
    S = usefulthings.S;
    N = usefulthings.N;
    HomeAngle  = usefulthings.HomeAngle;
    HomePosition  = usefulthings.HomePosition;
    ContextType  = usefulthings.ContextType;
    TargetAngle  = usefulthings.TargetAngle;
    field  = usefulthings.field;
    ntrial = usefulthings.ntrial;
    nsubj  = usefulthings.nsubj;
    fdir   = usefulthings.fdir;
    indB   = usefulthings.indB;
    indE   = usefulthings.indE;
    indT   = usefulthings.indT;
    indP   = usefulthings.indP;
    %-----
    
    %% Plot MPE for each subject, and group average
    statsblocks = 6;
    
    pre=5;  exp=150; trans=20; post=3;
    basephase = zeros(N,1); basephase(indB) = 1;
    exposurephase = zeros(N,1); exposurephase(indE) = 1;
    transferphase = zeros(N,1); transferphase(indT) = 1;
    postphase = zeros(N,1); postphase(indP) = 1;
    indbase = find((field(:,S.FIELD_CHANNEL)~=1) & basephase);
    indexp = find((field(:,S.FIELD_CHANNEL)~=1) & exposurephase);
    indtrans = find((field(:,S.FIELD_CHANNEL)~=1) & transferphase);
    indpost = find((field(:,S.FIELD_CHANNEL)~=1) & postphase);
    
    % smooth MPE by block (c=8)
    smooth_factor = 3;
    ind = find(field(:,S.FIELD_CHANNEL)~=1);                      % null and exposure trials give mpe measures
    c = 8;                                                        % non-channel trials per block
    mpeblock = mean(reshape(mpe(ind), c, length(mpe(ind))/c),1);
    mpeblock = reshape(mpeblock,length(mpeblock)/nsubj,nsubj);
    if group<5                                                   % if it's the planning only or imagine condition, get rid of some extra data on the ends
        mpeblock = mpeblock(2:end-1,:);
    end
    
    mpebase = mean(reshape(mpe(indbase), c, length(mpe(indbase))/c),1);
    mpebase = reshape(mpebase,length(mpebase)/nsubj,nsubj);
    stats.mpebaseline{group} = mpebase;
    stats.firstmpeblock{group} = mpeblock(pre+1:pre+statsblocks,:);
    stats.finalmpeblock{group} = mpeblock(pre+exp+1-statsblocks:pre+exp,:);
    
    % smooth MPE by 2 blocks in exposure phase (c=16)
    c = c*smooth_factor;
    mpeblockexp = mean(reshape(mpe(indexp), c, length(mpe(indexp))/c),1);
    mpeblockexp = reshape(mpeblockexp,length(mpeblockexp)/nsubj,nsubj);
    
    c = 8*4;  % can't smooth every 3 in transfer phase
    mpeblocktrans = mean(reshape(mpe(indtrans), c, length(mpe(indtrans))/c),1);
    mpeblocktrans = reshape(mpeblocktrans,length(mpeblocktrans)/nsubj,nsubj);
    
    % plot mpe per subject (smoothed per block only)
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
    
    % average across subjects
    mpeblock_mean = mean(mpeblock,2);
    mpeblock_se = std(mpeblock,0,2) ./ sqrt(size(mpeblock,2));
    mpeblockexp_mean = mean(mpeblockexp,2);
    mpeblockexp_se = std(mpeblockexp,0,2) ./ sqrt(size(mpeblockexp,2));
    mpeblocktrans_mean = mean(mpeblocktrans,2);
    mpeblocktrans_se = std(mpeblocktrans,0,2) ./ sqrt(size(mpeblocktrans,2));
    aftereffects{group} = mean(mpeblock(end-post+1:end-post+2,:),1);
    
    % plot average across all subjects in group
    figure(1000);
    subplot(1,2,1);
    shadeplot(1:pre, mpeblock_mean(1:pre), mpeblock_se(1:pre),'-',colours(group,:),0.3); hold all; % pre-exposure
    shadeplot(linspace(pre+1,exp+pre,length(mpeblockexp_mean)), mpeblockexp_mean, mpeblockexp_se,'-',colours(group,:),0.3); hold all;  % exposure
    shadeplot(linspace(pre+exp+1,exp+pre+trans,length(mpeblocktrans_mean)), mpeblocktrans_mean, mpeblocktrans_se,'-',colours(group,:),0.3); hold all;  % exposure
    shadeplot(exp+pre+trans+1:exp+pre+post+trans, mpeblock_mean(end-post+1:end), mpeblock_se(end-post+1:end),'-',colours(group,:),0.3); hold all;  % pre-exposure
    plot([0 pre+exp+trans+post],[0 0],'k');
    ylabel('MPE (cm)');
    xlabel('Block');
    axis([0 pre+exp+post+40 -1.5 3.5]);
    
    % plot the average after effects in a subpanel
    aftereffects_mean = mean(mean(mpeblock(end-post+1:end-post+2,:),1));
    aftereffects_se = std(mean(mpeblock(end-post+1:end-post+2,:),1)) ./ sqrt(nsubj);
    errorbar(pre+exp+trans+post+3+group*3, aftereffects_mean,aftereffects_se,'k'); hold on;
    plot(pre+exp+trans+post+3+group*3, aftereffects_mean, 'o','MarkerSize',7,'MarkerEdgeColor',colours(group,:), 'MarkerFaceColor',colours(group,:)); hold on;
    
    %% Plot adaptation for each subject on the same figure (1 per group)
    
    % smooth adaptation by block
    inclchannel = ((round(HomeAngle)==-180) | (round(HomeAngle==0)));
    inclchannelim = ((round(HomeAngle)==-180) | (round(HomeAngle==0)));
    basechannels =  find(((field(:,S.FIELD_CHANNEL)==1) & inclchannel) & (ContextType==0 & basephase));
    ind = find(((field(:,S.FIELD_CHANNEL)==1) & inclchannel) & (ContextType==0));
    indexcltrans = find(((field(:,S.FIELD_CHANNEL)==1) & inclchannel) & ((ContextType==0) & ~transferphase));
    indexp = find((field(:,S.FIELD_CHANNEL)==1) & (exposurephase & inclchannel));
    indtrans = find(((field(:,S.FIELD_CHANNEL)==1) & (ContextType==0)) & (transferphase & inclchannel));
    indtransim = find(((field(:,S.FIELD_CHANNEL)==1) & (ContextType==2)) & (transferphase & inclchannelim));
    
    incl0channel = (round(HomeAngle==0) & (ContextType==0));
    basechannels0 =  find(((field(:,S.FIELD_CHANNEL)==1) & incl0channel) & (ContextType==0 & basephase));
    ind0 = find((field(:,S.FIELD_CHANNEL)==1) & incl0channel);
    ind0excltrans = find(((field(:,S.FIELD_CHANNEL)==1) & incl0channel) & ~transferphase);
    ind0exp = find((field(:,S.FIELD_CHANNEL)==1) & (exposurephase & incl0channel));
    
    c = 2;      % channel trials per block (for block = 12 trials, trials balanced across 24 trials). Note that this is not including the transfer phase
    adaptationblock = mean(reshape(adaptation(indexcltrans), c, length(adaptation(indexcltrans))/c),1);  % smooth by block
    adaptationblock = reshape(adaptationblock,length(adaptationblock)/nsubj,nsubj);
    
    % find the baseline amount of adaptation
    adaptationbase = mean(reshape(adaptation(basechannels), c, length(adaptation(basechannels))/c),1);  % smooth by block
    adaptationbase = reshape(adaptationbase,length(adaptationbase)/nsubj,nsubj);
    adaptationbase0 = mean(reshape(adaptation(basechannels0), c, length(adaptation(basechannels0))/c),1);  % smooth by block
    adaptationbase0 = reshape(adaptationbase0,length(adaptationbase0)/nsubj,nsubj);
    
    % in the planning only and imagine conditions, get rid of some extra data on the ends so group plotting is comparable
    adaptationblock = adaptationblock(2:end-1,:);
    adaptation0block = reshape(adaptation(ind0excltrans),length(adaptation(ind0excltrans))/nsubj,nsubj);  % only 1 channel of this type per block
    adaptation0block = adaptation0block(2:end-1,:);
    
    stats.adaptbaseline{group} = adaptationbase;
    stats.adaptbaseline0{group} = adaptationbase0;
    stats.firstadaptblock{group} = adaptationblock(pre+1:pre+statsblocks,:);
    stats.finaladaptblock{group} = adaptationblock(pre+exp+1-statsblocks:pre+exp,:);
    stats.firstadapt0block{group} = adaptation0block(pre+1:pre+statsblocks,:);
    stats.finaladapt0block{group} = adaptation0block(pre+exp+1-statsblocks:pre+exp,:);
    
    % smooth adaptation by 2 blocks in exposure phase
    c = c*smooth_factor;
    adaptationblockexp = mean(reshape(adaptation(indexp), c, length(adaptation(indexp))/c),1);
    adaptationblockexp = reshape(adaptationblockexp,length(adaptationblockexp)/nsubj,nsubj);
    
    c = 4;   % each block in the transfer phase (20 of them) is 2 included channel trials in length
    adaptationblocktrans = mean(reshape(adaptation(indtrans), c, length(adaptation(indtrans))/c),1);
    adaptationblocktrans = reshape(adaptationblocktrans,length(adaptationblocktrans)/nsubj,nsubj);
    adaptationblocktransim = mean(reshape(adaptation(indtransim), c, length(adaptation(indtransim))/c),1);
    adaptationblocktransim = reshape(adaptationblocktransim,length(adaptationblocktransim)/nsubj,nsubj);
    
    adtrans = reshape(adaptation(indtrans),length(adaptation(indtrans))/nsubj,nsubj);
    adtransim = reshape(adaptation(indtransim),length(adaptation(indtransim))/nsubj,nsubj);
    
    % plot by subject in group
    figure();
    clear adaptationsubjfinalblock
    for i=1:nsubj
        subplot(1,nsubj,i);
        patch(P); hold all;
        plot(linspace(1,pre+exp,size(adaptationblock(1:pre+exp,:),1)),adaptationblock(1:pre+exp,i),'Color',colours(group,:)); hold on;
        plot(linspace(pre+exp+1,pre+exp+trans,size(adaptationblocktrans,1)),adaptationblocktrans(:,i),'Color',colours(group,:)); hold on;
        plot(linspace(pre+exp+1,pre+exp+trans,size(adaptationblocktransim,1)),adaptationblocktransim(:,i),'Color',colours(group+1,:)); hold on;
        plot(linspace(pre+exp+trans+1,pre+exp+trans+post,size(adaptationblock(pre+exp+1:end,:),1)),adaptationblock(pre+exp+1:end,i),'Color',colours(group,:)); hold on;
        plot([0 pre+exp+trans+post],[0,0],'k');
        axis([0 pre+exp+trans+post -20 100]);
        if i==1
            xlabel('Trials');
            ylabel('Adapatation (%)');
        end
        adaptationsubjfinalblock(i) = mean(adaptationblock(pre+exp+1-statsblocks:pre+exp,i));
        adaptationsubjfinalblock_sd(i) = std(adaptationblock(pre+exp+1-statsblocks:pre+exp,i));
        adaptationsubjfinalblocktrans(i) = mean(adtrans(1:4,i));  % first 4 blocks of transfer phase
        adaptationsubjfinalblocktrans_sd(i) = std(adtrans(1:4,i));
        adaptationsubjfinalblockim(i) = mean(adtransim(1:4,i));  % first 4 blocks of transfer phase
        adaptationsubjfinalblockim_sd(i) = std(adtransim(1:4,i));
    end
    strng = sprintf('Adaptation per subject, group %d', group);
    suptitle(strng);
    
    % average across subjects
    adaptblock_mean = mean(adaptationblock,2);
    adaptblock_se = std(adaptationblock,0,2) ./ sqrt(size(adaptationblock,2));
    adaptblockexp_mean = mean(adaptationblockexp,2);
    adaptblockexp_se = std(adaptationblockexp,0,2) ./ sqrt(size(adaptationblockexp,2));
    adaptblocktrans_mean = mean(adaptationblocktrans,2);
    adaptblocktrans_se = std(adaptationblocktrans,0,2) ./ sqrt(size(adaptationblocktrans,2));
    adaptblocktransim_mean = mean(adaptationblocktransim,2);
    adaptblocktransim_se = std(adaptationblocktransim,0,2) ./ sqrt(size(adaptationblocktransim,2));
    
    % plot average across all subjects in group
    figure(1000);
    subplot(1,2,2);
    shadeplot(1:pre, adaptblock_mean(1:pre), adaptblock_se(1:pre),'-',colours(group,:),0.3); hold all; % pre-exposure
    shadeplot(linspace(pre+1,exp+pre,length(adaptblockexp_mean)), adaptblockexp_mean, adaptblockexp_se,'-',colours(group,:),0.3); hold all;  % exposure
    h1=shadeplot(linspace(pre+exp+1,exp+trans+pre,length(adaptblocktrans_mean)), adaptblocktrans_mean, adaptblocktrans_se,'-',colours(group,:),0.3); hold all;  % exposure
    h2=shadeplot(linspace(pre+exp+1,exp+trans+pre,length(adaptblocktransim_mean)), adaptblocktransim_mean, adaptblocktransim_se,'-',colours(group+2,:),0.3); hold all;  % exposure
    shadeplot(exp+pre+trans+1:exp+pre+trans+post, adaptblock_mean(end-post+1:end), adaptblock_se(end-post+1:end),'-',colours(group,:),0.3); hold all;  % post-exposure
    plot([0 pre+exp+post],[0 0],'k');
    ylabel('Adaptation (%)');
    xlabel('Block');
    
    % plot the average final adaptation level in a subpanel
    
    % execution trials - end of exposure phase
    finaladapt_mean = mean(adaptationsubjfinalblock);
    finaladapt_se = std(adaptationsubjfinalblock) ./ sqrt(nsubj);
    errorbar(pre+exp+trans+post+3+group*3, finaladapt_mean,finaladapt_se,'k'); hold on;
    plot(pre+exp+trans+post+3+group*3, finaladapt_mean, 'o','MarkerSize',7,'MarkerEdgeColor',colours(group,:), 'MarkerFaceColor',colours(group,:)); hold on;
    
   
    % imagining trials - first 4 blocks of transfer phase
    finaladapt_mean = mean(adaptationsubjfinalblockim);
    finaladapt_se = std(adaptationsubjfinalblockim) ./ sqrt(nsubj);
    errorbar(pre+exp+trans+post+6+group*3, finaladapt_mean,finaladapt_se,'k'); hold on;
    plot(pre+exp+post+trans+6+group*3, finaladapt_mean, 'o','MarkerSize',7,'MarkerEdgeColor',colours(group,:), 'MarkerFaceColor','w'); hold on;
    
    axis([0 pre+exp+trans+post+20 -20 60]);
    legend([h1',h2'],{'Execution trials','Imagery trials'})
    
    %% Plot the trajectories for each group, at different experiment phases
    % Note: this takes a long time to execute
    
    if FLAGS.plotextra
        prephase=[5,6];
        earlyexp=[1,2];
        lateexp=[149,150];
        postphase=[7,8];
        
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
        homeangles = unique(round(HomeAngle));
        labels = {'pre','early','late','post'};
        for angle=1:length(homeangles)
            i = homeangles(angle);
            
            n1 = find(((field(:,S.FIELD_NULL)==1) & (fdir==1))  & (round(HomeAngle)==i)); %leftward target    0 degrees
            n2 = find(((field(:,S.FIELD_NULL)==1) & (fdir==-1))  & (round(HomeAngle)==i)); %rightward target
            n1 = reshape(n1,length(n1)/nsubj,nsubj);
            n2 = reshape(n2,length(n2)/nsubj,nsubj);
            trials.pre(:,:,1) = n1(prephase,:);
            trials.post(:,:,1) = n1(postphase,:);
            trials.pre(:,:,2) = n2(prephase,:);
            trials.post(:,:,2) = n2(postphase,:);
            
            e1 = find(((field(:,S.FIELD_VISCOUS)==1) & (fdir==1)) & (round(HomeAngle)==i)); %leftward target    0 degrees
            e2 = find(((field(:,S.FIELD_VISCOUS)==1) & (fdir==-1)) & (round(HomeAngle)==i)); %rightward target
            e1 = reshape(e1,length(e1)/nsubj,nsubj);
            e2 = reshape(e2,length(e2)/nsubj,nsubj);
            trials.early(:,:,1) = e1(earlyexp,:);
            trials.late(:,:,1) = e1(lateexp,:);
            trials.early(:,:,2) = e2(earlyexp,:);
            trials.late(:,:,2) = e2(lateexp,:);
            
            % average the trajectories across blocks and subjects
            for j=1:length(labels)
                for k=1:2    % number of different force fields
                    tmp = trials.(labels{j})(:,:,k);
                    x = mean(reshape(squeeze(F.rpos(1,:,tmp(:)))',[2,nsubj,sz]),1);
                    y = mean(reshape(squeeze(F.rpos(2,:,tmp(:)))',[2,nsubj,sz]),1);
                    if nsubj~=1
                        trajectory_x.(labels{j})(k,:,:) = [squeeze(mean(x)),squeeze(std(x)./sqrt(nsubj))];
                        trajectory_y.(labels{j})(k,:,:) = [squeeze(mean(y)),squeeze(std(y)./sqrt(nsubj))];
                        subplot(1,4,j);
                        if k==1
                            shadedTrajectory(squeeze(trajectory_x.(labels{j})(k,:,:))', squeeze(trajectory_y.(labels{j})(k,:,:))', [0.8 0.2 0.2], 0.1); hold on;
                        else
                            shadedTrajectory(squeeze(trajectory_x.(labels{j})(k,:,:))', squeeze(trajectory_y.(labels{j})(k,:,:))', [0.2 0.2 0.8], 0.1); hold on;
                        end
                    else % plot just the excluded subject's mean trajectory
                        trajectory_x.(labels{j})(k,:) = squeeze(x)';
                        trajectory_y.(labels{j})(k,:) = squeeze(y)';
                        subplot(1,4,j);
                        if k==1
                            plot((1000+trajectory_x.(labels{j})(k,:))', (1000+trajectory_y.(labels{j})(k,:))', 'Color',[0.8 0.2 0.2]); hold on;
                        else
                            plot((1000+trajectory_x.(labels{j})(k,:))', (1000+trajectory_y.(labels{j})(k,:))', 'Color',[0.2 0.2 0.8]); hold on;
                        end
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
    fieldnm = {'fixationtransfer'};
    
    % plot average mental chronometry across all subjects in group
    figure(2000);
    subplot(1,2,1);
    cf1 = plot([0 2.5],[0 2.5],'-','Color',[0.5 0.5 0.5]); hold on;
    cf2 = plot(averageexecutionduration,averageimagineduration,'.', 'MarkerSize',15, 'Color', colours(group,:)); hold on;
    xlabel('Execution duration (s)');
    ylabel('Imagining duration (s)');
    axis([0 1.5 0 1.5])
    
    imaginedurationgroups(group,:,:) = averageimagineduration;
    executiondurationgroups(group,:,:) = averageexecutionduration;
    
    im = mean(imaginedurationgroups,1);
    ex = mean(executiondurationgroups,1);
    subplot(1,2,2);
    bar(1,mean(im),'FaceColor',imaginecolour); hold on;
    errorbar(1,mean(im),std(im)./sqrt(nsubj),'Color',imaginecolour); hold on;
    bar(2,mean(ex),'FaceColor',executecolour); hold on;
    errorbar(2,mean(ex),std(ex),'Color',executecolour); hold on;
    ylabel('Movement duration (s)');
    ax=gca;
    ax.XTick = [1,2];
    ax.XTickLabel = {'Imagined','Executed'};
    ax.XTickLabelRotation = 50;
    
    %% Plot correlations between final adaptation level and MI scores
    figure(3000);
    subplot(1,3,1);  % plot final adaptation level vs. MIQ-score
    h1=plot(imagery.MIQscore,adaptationsubjfinalblock,'xk'); hold on;
    h2=plot(imagery.MIQkinesthetic, adaptationsubjfinalblock, '.b'); hold on;
    h3=plot(imagery.MIQvisual, adaptationsubjfinalblock, '.c'); hold on;
    
    XMIQ(:,group) = imagery.MIQscore;
    xlabel('Mean MIQ-score');
    ylabel('Adaptation level (%)');
    legend([h1',h2',h3'], {'Both MIQ','MIQ kinesthetic','MIQ visual'});
    
    subplot(1,3,2);  % plot final adaptation level vs. MIQ-maintenance frequency. This looks good.
    plot(squeeze(mean(imagery.MImaintain(:,:,2),2)),adaptationsubjfinalblock,'.','MarkerSize',15, 'Color', colours(group,:)); hold on;
    for i=1:nsubj
        errorbar(squeeze(mean(imagery.MImaintain(i,:,2),2)),adaptationsubjfinalblock(i), adaptationsubjfinalblock_sd(i),'Color', colours(group,:)); hold on;
    end
    
    tmp = mean(imagery.MImaintain(:,:,2),2);
    Xmaintain(:,group) = tmp;
    xlabel('Mean MI-maintenance frequency');
    ylabel('Adaptation level (%)');
    
    subplot(1,3,3);  % plot final adaptation level vs. MIQ-maintenance ease.
    plot(squeeze(mean(imagery.MImaintain(:,:,1),2)),adaptationsubjfinalblock,'.','MarkerSize',15, 'Color', colours(group,:)); hold on;
    for i=1:nsubj
        errorbar(squeeze(mean(imagery.MImaintain(i,:,1),2)),adaptationsubjfinalblock(i), adaptationsubjfinalblock_sd(i),'Color', colours(group,:)); hold on;
    end
    
    tmp = mean(imagery.MImaintain(:,:,1),2);
    Xmaintainease(:,group) = tmp;
    xlabel('Mean MI-maintenance ease');
    ylabel('Adaptation level (%)');
    finaladapt(:,group) = adaptationsubjfinalblock';
    
    % regress and plot correlation
    XMIQtotal = reshape(XMIQ, numel(XMIQ),1);
    Xmaintaintotal = reshape(Xmaintain, numel(Xmaintain),1);
    Xmaintaineasetotal = reshape(Xmaintainease, numel(Xmaintainease),1);
    y = reshape(finaladapt, numel(finaladapt),1);
    
    subplot(1,3,1);
    [b,bint,~,~,imagestats]=regress(y,[XMIQtotal, ones(size(XMIQtotal))]);                % b = [slope, intercept]
    [r,pim] = corrcoef(y,XMIQtotal);
    plot(0:10,b(1).*[0:10]+b(2), '--b'); hold on;
    strng = sprintf('R = %.2f, p = %.4f', r(1,2),pim(1,2));
    text(4,20,strng);
    
    subplot(1,3,2);
    [b,bint,~,~,imagestats]=regress(y,[Xmaintaintotal, ones(size(Xmaintaintotal))]);                % b = [slope, intercept]
    [r,pim] = corrcoef(y,Xmaintaintotal);
    plot(0:3,b(1).*[0:3]+b(2), '--b'); hold on;
    strng = sprintf('R = %.2f, p = %.4f', r(1,2), pim(1,2));
    text(1,20,strng);
    
    subplot(1,3,3);
    [b,bint,~,~,imagestats]=regress(y,[Xmaintaineasetotal, ones(size(Xmaintaineasetotal))]);                % b = [slope, intercept]
    [r,pim] = corrcoef(y,Xmaintaineasetotal);
    plot(0:10,b(1).*[0:10]+b(2), '--b'); hold on;
    strng = sprintf('R = %.2f, p = %.4f', r(1,2), pim(1,2));
    text(1,20,strng);
    
    %% Save the data for comparison to the other groups
    if group==1 % we've already saved the rest
        v2struct(stats);
        noimtransfergrouplearning = v2struct(pre,post,exp,trans,adaptblock_mean, ...
            adaptblock_se, adaptblocktrans_mean, adaptblocktrans_se,adaptblocktransim_mean, ...
            adaptblocktransim_se, finaladapt_mean, finaladapt_se, adaptationsubjfinalblock, ...
            adaptationsubjfinalblockim, adaptationsubjfinalblocktrans, adaptblockexp_mean, ...
            adaptblockexp_se, mpeblock_mean, mpeblock_se, mpeblocktrans_mean, mpeblocktrans_se, ...
            mpeblockexp_mean, mpeblockexp_se, mpeblock, mpebaseline, adaptationbase0, adaptationbase, ...
            adaptbaseline, adaptbaseline0, aftereffects, aftereffects_mean, aftereffects_se,...
            firstadaptblock, finaladaptblock, firstadapt0block, finaladapt0block);
        save('noimtransferlearning','noimtransfergrouplearning');
    end
    
    clearvars -except finaladapt  Xtiming Xmaintainease Xmaintain XMIQ Adaptation AdaptationDetail P MPE stats ngroups nsubj ntrial N S FLAGS Timing Experiment colours FrameData nmpe_scalefactor aftereffects fh fh2 A M MI speed_series speed_error imduration execduration cf1 cf2 cf3 chronometry imaginedurationgroups executiondurationgroups adaptationsubjfinalblockim adaptationsubjfinalblocktrans
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

clear pre exp post fh fh2 P Pg ntrial N
%------------------------------
%% Statistical Analysis - learning
% Used pre-exposure (first 6) blocks and last 6 blocks of exposure data for epochs
groupname = {'noimagerytransfer','imagerytransfer'};
for group=1:ngroups                    % within-group learning
    
    % compare pre and late exposure phase
    [~, p.adaptation{group},~,statstable.adaptationtable{group}] = ttest(squeeze(mean(stats.adaptbaseline{group}))', squeeze(mean(stats.finaladaptblock{group}))');
    [~, p.adaptation0{group},~,statstable.adaptation0onlytable{group}] = ttest(squeeze(mean(stats.adaptbaseline0{group}))', squeeze(mean(stats.finaladapt0block{group}))');
    adaptstats.(groupname{group}) = [squeeze(mean(stats.adaptbaseline{group}))', squeeze(mean(stats.finaladaptblock{group}))'];
    adaptstats0only.(groupname{group}) = [squeeze(mean(stats.adaptbaseline0{group}))', squeeze(mean(stats.finaladapt0block{group}))'];
    
    mpestats.(groupname{group}) = [squeeze(mean(stats.firstmpeblock{group}))', squeeze(mean(stats.finalmpeblock{group}))'];
    
    % Assess after-effect magnitude (paired ttest compared to pre-exposure)
    [~,tmp_p,~,tmp_s] = ttest(aftereffects{group}',mean(stats.mpebaseline{group})');
    p.aftereffects.(groupname{group}) = tmp_p;
    statstable.aftereffects.(groupname{group}) = tmp_s;
    
    
    mn = mean(mean(stats.finaladaptblock{group})- mean(stats.adaptbaseline{group}));
    se = std(mean(stats.finaladaptblock{group})- mean(stats.adaptbaseline{group})) ./ sqrt(8);
    sprintf('Diff adaptation %s: %.2f  +/- %.2f', groupname{group}, mn, se)
    
    mn = mean(aftereffects{group}'- mean(stats.mpebaseline{group})');
    se = std(aftereffects{group}'- mean(stats.mpebaseline{group})') ./ sqrt(8);
    sprintf('Diff aftereffects %s: %.2f  +/- %.4f', groupname{group}, mn, se)
end

%% Statistical analysis - Baseline kinematics by Left vs Right followthrough
% check for baseline difference for L vs R follow targets for for 5
% kinematic variables in each group. Note that this give 10 datapoints per
% Left or Right condition per subject in groups1-4, and 8 datapoints per L or R per sub in group 5.
% (use mean value per subject to contribute to group stats & consider just the 0 degree starting position)
for i=1:ngroups
    FIELD_NULL = Experiment{i}.S.FIELD_NULL; % same for all groups
    ind = (Experiment{i}.field(:,FIELD_NULL) & Experiment{i}.HomeAngle==0) & (Experiment{i}.HomePosition(:,2)==[-16]); % hack for choosing 0degree start position. Seems redundant?
    indL = ind & Experiment{i}.TargetAngle==45;
    indR = ind & Experiment{i}.TargetAngle==-45;
    
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

%% Compare initial after-effect magnitudes between groups to assess transfer
% we don't have channel trials on imagining trials, so compare learning on
% imagining trials to planning trials by using after-effects


toc
