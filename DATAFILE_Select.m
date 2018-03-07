
function D=DATAFILE_Select(P,samp,trials)
%D= DATAFILE_Select(P,samp,trials)
% Extracts the trial numbers specified in the vector samp from the data structure P
% trials are renumbered sequentially  - specify total  number of trials if there
% is not P.Trials as part of structure
% samp can be a logical inex or set of trial numbers

if nargin==2
    trials=P.Trials; 
end

%get fieldnames
fu=fieldnames(P);
%loop over fieldnames
for k=1:length(fu)
  
    %only perform extraction on non-structures
    if ~isstruct(P.(fu{k}))
        
        %pull out matrices
        c1= P.(fu{k});
        s1=size(c1);
        
        %find the trial dimension
        i=find(s1==trials,1,'first');
        
        %now extract depending on whether the matrix is 2 or 3D
        if length(s1)==2
          if isempty(i)
            D.(fu{k})=c1;
          elseif i==1
            D.(fu{k})=c1(samp,:);
          elseif i==2
            D.(fu{k})=c1(:,samp);
          end
        elseif length(s1)==3
          if i==1
            D.(fu{k})=c1(samp,:,:);
          elseif i==2
            D.(fu{k})=c1(:,samp,:);
          elseif i==3
            D.(fu{k})=c1(:,:,samp);
          end
        end
    else %if a structure, then use this function recursively but give the number of trial
    D.(fu{k})=DATAFILE_Select(P.(fu{k}),samp,trials);
    
    end
end

%finally set the number of trials and trial numbers
if nargin==2
    if isequal(unique(samp),[0 1]')
        D.Trials=sum(samp);
        D.ResampTrialNumber=(1: D.Trials)';
        
    else
        D.Trials=length(samp);
        D.ResampTrialNumber=(1: D.Trials)';
    end
end
