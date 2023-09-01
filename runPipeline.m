% We will be working in this directory, which contains a subdirectory for
% each subject
ind = '/Volumes/Quattro/ct_alex/acute2/rerun';
addpath(genpath('/Volumes/Quattro/ct_alex/acute2'))
% First, we will run the prebet pipeline...
d = dir(ind);
for i = i:length(d)
    disp(['Working on prebetting of subject ' num2str(i) ' of ' num2str(length(d))])
    % get files within directory
    fd = [d(i).folder '/' d(i).name];
    d2 = dir([fd '/*.nii']);
    tmp = {d2.name};
    
    % reorganize files so that head is first
    id = find(contains(tmp,'Head'));
    if isempty(id)
        id = find(contains(tmp,'CTP'));
        if isempty(id)
            id = find(contains(tmp,'CT'));
        end
    end
        
    id2 = setdiff([1:length(tmp)],id);
    fnm = [tmp(id) tmp(id2)];
    
    % run prebet
    ct_rgb(strcat([fd '/'],fnm),'prebet');    
end

% Now let's do some manual betting for participants that failed...
% From some testing I did earlier, it looks like 0.005 works well for some
% folks that had most of their brain missing. Just below that can be our
% lower limit. Upper limit can be just under 0.1 (which is what the
% original data used)...skip below this matlab code for something faster.
% Also note, I ended up running some extra -f values not shown below so
% ymmv.
bets = [0.001 0.003 0.005 0.007 0.009 0.01 0.03 0.05 0.06 0.07 0.08 0.09];
for i = i:length(d)
    disp(['Working on betting of subject ' num2str(i) ' of ' num2str(length(d))])
    fd = [d(i).folder '/' d(i).name];
    d2 = dir([fd '/*.nii']);
    tmp = {d2.name};
    idz = find(startsWith(tmp,'z'));
    fnm = strcat([fd '/'],tmp{idz});
    tmpo = num2str(bets(j));
    fnmo = strcat([fd '/'],['b' tmp{idz}(1:end-4) '_f_' tmpo(3:end) '.nii']);
    for j = 1:length(bets)
        cmd = ['/usr/local/fsl/bin/bet ' fnm ' ' fnmo ' -f ' num2str(bets(j)) ' -g 0'];
        system(cmd);
        fprintf(cmd);
    end
end
   
% Geeze, this is taking forever, it's much faster in bash. Here's some bash
% code...the bottom portion includes radius option example for tough subs.
% bd=/Volumes/Quattro/ct_alex/acute2/rerun
% cd $bd
% for subs in *; do
% echo ${subs}
%     fnm=`ls $bd/${subs}/z*.nii`
%     on=`basename $fnm .nii`
%     for j in 0.002 0.003 0.005 0.007 0.009 0.01 0.03 0.05 0.06 0.07 0.08 0.09; do
%         echo ${j}
%         fnmo=$bd/${subs}/b"$on"_f_${j:2}_R.nii
%         bet $fnm $fnmo -f ${j} -g 0 -R
%     done
% done
%
% Here is code for a subject that required some additional
% intevention (see bet info for the options here)
% bd=/Volumes/Quattro/ct_alex/acute2/rerun
% cd $bd
% for subs in RGE7333; do 
%      echo subject is ${subs};
%      fnm=`ls $bd/${subs}/z*.nii`;     
%      on=`basename $fnm .nii`;     
%      for j in 0.000001; do 
%          for k in 175; do             
% 	     echo threshold ${j} using radius ${k};
% 	     fnmo=$bd/${subs}/b"$on"_f_${j:2}_R_${k}_r_TESTING.nii;             
% 	     bet $fnm $fnmo -f ${j} -g 0 -r ${k} -c 234 327 133; 
% 	 done;     
%      done; 
% done;

% The next step is to id the best bet result for each subject. In some
% cases, I used neck cleanup (-B) but this typically took a while to run.
% Using robust bet (-R) worked just as well so I stuck with that.

% And now that we have sorted through our betted brains and we have the
% best ones for each subject, we just need to run the post-bet portion
% of our pipeline

% This is some code just to cleanup our directories and remove files we don't need...
% bd=/Volumes/Quattro/ct_alex/acute2/rerun
% cd $bd
% for subs in *; do
%     echo ${subs}
%     rm $bd/${subs}/q*.nii
%     rm $bd/${subs}/*.tab
%     rm $bd/${subs}/scalar*.nii
%     rm $bd/${subs}/atlas*
% done
    
% And now getting to the majority of the pipeline...
for i = 1:length(d)
    disp(['Working on subject ' num2str(i) ' of ' num2str(length(d))])
    % get files within directory
    fd = [d(i).folder '/' d(i).name];
    d2 = dir([fd '/*.nii']);
    tmp = {d2.name};
    
    % reorganize files so that head is first
    tf = startsWith(tmp,'z');
    %z = tmp(tf);
    tmp(tf) = [];

    tf = startsWith(tmp,'bz');
    bz = tmp{tf};
    tmp(tf) = [];
    
    id = find(contains(tmp,'Head'));
    if isempty(id)
        id = find(contains(tmp,'CTP'));
        if isempty(id)
            id = find(contains(tmp,'CT'));
        end
    end
        
    id2 = setdiff([1:length(tmp)],id);
    fnm = [tmp(id) tmp(id2)];
    
    ct_rgb(strcat([fd '/'],fnm),'betted',strcat([fd '/'],bz));
end


for i = 1:length(d)
    disp(['Working on subject ' num2str(i) ' of ' num2str(length(d))])
    % get files within directory
    fd = [d(i).folder '/' d(i).name];
    d2 = dir([fd '/*.nii']);
    tmp = {d2.name};

    tf = startsWith(tmp,'bz');
    bz = tmp{tf};
    tmp(tf) = [];

    ct_rgb({strcat([fd '/'],bz)},'normalize',4,8);
end

d = dir('/Volumes/Quattro/ct_alex/acute2/normStragglers');
for i = 1:length(d) 
    disp(['Working on subject ' num2str(i) ' of ' num2str(length(d))])
    % get files within directory
    fd = [d(i).folder '/' d(i).name];
    d2 = dir([fd '/*.nii']);
    tmp = {d2.name};

    tf = startsWith(tmp,'bz');
    bz = tmp{tf};
    tmp(tf) = [];
    ct_rgb({strcat([fd '/'],bz)},'normalize',4,8);
end

d = dir('/Volumes/Quattro/ct_alex/acute2/normStragglers');
for i = 1:length(d) % 27 failed...
    disp(['Working on subject ' num2str(i) ' of ' num2str(length(d))])
    % get files within directory
    fd = [d(i).folder '/' d(i).name];
    d2 = dir([fd '/*.nii']);
    tmp = {d2.name};

    tf = startsWith(tmp,'bz');
    bz = tmp{tf};
    tmp(tf) = [];
    % run prebet
    ct_rgb({strcat([fd '/'],bz)},'normalize',4,8);
end

d = dir('/Volumes/Quattro/ct_alex/acute2/snrStragglers');
for i = 1:length(d)
    disp(['Working on subject ' num2str(i) ' of ' num2str(length(d))])
    % get files within directory
    fd = [d(i).folder '/' d(i).name];
    d2 = dir([fd '/*.nii']);
    tmp = {d2.name};
    
    % reorganize files so that head is first
    tf = startsWith(tmp,'z');
    %z = tmp(tf);
    tmp(tf) = [];

    tf = startsWith(tmp,'bz');
    bz = tmp{tf};
    tmp(tf) = [];
    
    id = find(contains(tmp,'Head'));
    if isempty(id)
        id = find(contains(tmp,'CTP'));
        if isempty(id)
            id = find(contains(tmp,'CT'));
        end
    end
        
    id2 = setdiff([1:length(tmp)],id);
    fnm = [tmp(id) tmp(id2)];
    fnm{1} = bz;
    
    % run prebet
    ct_rgb(strcat([fd '/'],fnm),'betted');
end

