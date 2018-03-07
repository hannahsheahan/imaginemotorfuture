function D=DATAFILE_Append(P,Q,pn)
%D=DATAFILE_Append(P,Q)
%Appends trial structures P and Q into output D throwing away any
%components that only appear in one structure or that have different
%dimensions in the two structures. The append acts so that D appears as if
%experiments P and Q were run sequentially within the same session. If one
%of the structures is empty the other is simply returned - useful for looping  
%e.g.
%
% FileList = { 'file1' 'file2' 'file3'}
% FileCount = size(FileList,2);
%
%    D={};
%    for k=1:FileCount
%        tmp=DATAFILE_Load(sprintf('%s',FileList{k}));
%        D = DATAFILE_append(D,tmp);
%        clear tmp;
%    end
%
% To avoid appending any substructures e.g. FrameData (structures can get very large)
% use D=DATAFILE_append(P,Q,-1) in the loop 
  
%if we are inside a structure we need to tell it the trial length pn of one
%file so it can find correct dimension to search for in concatenation
if isempty(P),  D=Q;, return, end
if isempty(Q),  D=P;, return, end

%determines whether the frame data is concatenated
if nargin==3 &  pn<0
    frame=0;
else
    frame=1;
end

if nargin==2 | frame==0
    pn=P.Trials;
    pq=Q.Trials;
    D.Trials=pn+pq;  %total number of trials
end

%get fieldnames 
fp=fieldnames(P);
fq=fieldnames(Q);

%find unique set of fieldnames
fu= unique([fq; fp]);

%turn structures to cells
Pc=struct2cell(P);
Qc=struct2cell(P);

%loop over fieldnames
for k=1:length(fu)
    %find index of field names 
    
    ip=find(ismember(fp, fu{k})==1);
    iq=find(ismember(fq, fu{k})==1);
    
    %only concatenate if exist in both
    if isempty(ip) | isempty (iq)
        fprintf('WARNING: %s is  not in one of the structures ... skipping\n',fu{k})
        continue
    end
    
    %only perform concatenation on non-structures
    if eval(sprintf('~isstruct(%s.%s)','P',fu{k}))
        
        %pull out matrices
        c1= eval(sprintf('%s.%s','P',fu{k}));
        c2= eval(sprintf('%s.%s','Q',fu{k}));
        s1=size(c1);
        s2=size(c2);
        
        if numel(s1)~=numel(s2)
            fprintf('WARNING: %s different dimensions ... skipping\n',fu{k})
        elseif numel(c1)>1 %only concatenate matrices
            %find the trial dimension
            i=find(s1==pn);
            
            %set to zero in case they are different as we will pad for
            %trials of different length
            s1(i)=0;
            s2(i)=0;
                    
            %find dimensions that disagree
            f=find((s1==s2)==0);
            
            %pad with zero to same size
            if ~isempty(f)
                pad=max(s1-s2,s1*0);
                c2=padarray(c2,pad,NaN,'post');
                
                pad=max(s2-s1,s1*0);
                c1=padarray(c1,pad,NaN,'post');
            end
            
            %concatenate on trial dimension
            
            %assign output
            eval(sprintf('D.%s=cat(i,c1,c2);',fu{k}));
        end
    elseif frame %if a structure, then use this function recursively but give the number of trial for first strcture
        if strcmp(fu{k},'Files')
            eval(sprintf('D.%s=[P.%s Q.%s];',fu{k},fu{k},fu{k}));
        else
            eval(sprintf('D.%s=DATAFILE_Append(P.%s,Q.%s,%i);',fu{k},fu{k},fu{k},pn));
        end
    end
end



