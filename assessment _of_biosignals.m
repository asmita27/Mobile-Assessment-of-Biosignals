% Author: Asmita Narang, Martina Köhler, Nehchal Kaur
% Please start running the EEGlab by changing the workplace to wherever you
% have installed eeglab, so you can see the eeglab.m file. click on the
% command window below and type in eeglab so that the application can
% start. Once it is running (it has stopped running if the run button above in
% the editor tab is green again). Now the code below is executable. 

%% The three E's: EEG Neural processing of emojis, human faces, emoticons and words 
%  the following steps were used for pre processing the EEG data 
%  Step 1 : Read Files from the directory (you have to always specify the directory which has all the subject files)
%  Step 2 : Read EEG Files recorded in brain vision recorder (dynamic)
%  Step 3 : Filter data - Band pass filter (0.5Hz - 30Hz)
%  Step 4 : Artefact Removal (done automatically)
%           Output which channels are removed and keep record of the eye channels  
%  Step 5 : Epochs
%  Step 6 : Baseline Correction
%  Step 7 : Ocular Correction using Gratton algorithm, see step for
%           blinkcriterion and blinklength (hyperparameters that can be
%           changed
%  Step 8 : Average for every subject for its differnt variables 
%  Step 9 : Re-referencing using average of all the channels excluding eye
%           electrode 
%  Step 10: Grand average for all the subjects 

%% Step 1: Read Filenames
% Change directory  to wherever your *.eeg and *.vhdr files are: It's
% important that all files are in the same folder, otherwise they are not
% processed. Double check in the workplace files variable, whether all data
% files were picked up! If you use windows, you can simply type in the
% absolute file and put it into a string like: 'U:\11.Semester\EEG_Project'
% If the path has (a) space(s) somewhere, use double backslashes like:
% 'C:\\Users\\Martina Koehler\\Documents\\MATLAB' or try using "" to
% enclose the string. 

files = dir(fullfile('/Users/asmitanarang/Desktop/untitled folder 2/ee/Hands-On/','*vhdr'));
num_files = length(files) + 1; % number of files + 1 because of structure of while loops later

%% Step 2: Read files dynamically 
% Reading all the files and making nameArray which will have EEG1, EEG2...
% files according to the number of files found in directory 

count = 1;
while(count ~= num_files)% loops through the data of all participants
    N = num2str(count);
    nameArray(count,1) =  {strcat('eeg',N)};
    Data.(nameArray{count,1})= pop_loadbv(files(count,1).folder, files(count,1).name,...
                               [] ,[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 ...
                               22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 ...
                               42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 ...
                               62 63 64 65 66]);
    
    count = count+1;
end


%% Step 3: Filter data ~ Bandpass filtering

count = 1;
while(count ~= num_files)
    % creates a field within Data called eeg1, eeg2,... for each
    % participant: Data.eeg1, Data.eeg2,....
    Data.(nameArray{count,1})= eeg_checkset(Data.(nameArray{count,1}));
    % creates a new field in each Data.eeg# that contains the bandbass
    % filtering: saved in Data.eeg1.filtered, Data.eeg2.filtered,...
    Data.(nameArray{count,1}).filtered = pop_eegfiltnew(Data.(nameArray{count,1}),...
        'locutoff',0.5,'hicutoff',30,'plotfreqz',1);
    count = count + 1;
end


%% Step 4: Artefact Removal 

count = 1;
while(count ~= num_files)
    s = strcat('Clearning up file...' ,nameArray{count,1});
    disp(s)
    % creates fields: Data.eeg1.artefact_rem, Data.eeg2.artefact_rem,....
    Data.(nameArray{count,1}).artefact_rem = clean_artifacts(Data.(nameArray{count,1}).filtered); 
    count = count + 1; 
end

%% outputs which channels were removed:
count = 1;
while(count ~= num_files)
    filtered_chan = string(char(Data.(nameArray{count,1}).filtered.chanlocs(:).labels));
    artefact_rem_channels = string(char(Data.(nameArray{count,1}).artefact_rem.chanlocs(:).labels));

    count2 = 1;
    s = strcat('removed channels from file...', nameArray{count,1});
    disp(s);
    removed_chan = [];
    while (count2 <= length(filtered_chan))
        comp = filtered_chan{count2};
        %disp(comp)
        if(~ismember(artefact_rem_channels, comp))
            disp(filtered_chan{count2});
            removed_chan = [string(removed_chan), string(filtered_chan{count2})];
        end
        count2 = count2 + 1;
    end
    
    % find index of either HEOG or VEOG:
    % takes VEOG-channel as a default
    if (any(strcmp({Data.(nameArray{count,1}).artefact_rem.chanlocs.labels}, 'VEOG')==1))
        eog_index.(nameArray{count,1}) = ...
            find(strcmp({Data.(nameArray{count,1}).artefact_rem.chanlocs.labels}, 'VEOG')==1);
    elseif (any(strcmp({Data.(nameArray{count,1}).artefact_rem.chanlocs.labels}, 'HEOG')==1))
        eog_index.(nameArray{count,1}) = ...
            find(strcmp({Data.(nameArray{count,1}).artefact_rem.chanlocs.labels}, 'HEOG')==1);
    else
        disp('No EOG channel left after automatic artefact removal!')
        return
    end
    count = count + 1;
end


%% Step 5 + 6 + 7: Epochs + Baseline Correction + Occular Correction (RIGHT NOW: WRONG ORDER!)
% OUTPUT = It makes separate structure inside the structure of the cleaned data file 
% Segments epochs on the basis of markers and then performs occular correction and 
% then baseline correction 

% everything above 10volt can be a possible blink
blinkcritvolt = 10;
% everything above a duration of 400ms is detected as blink
blinkcritwin = 400;
count = 1;

while(count ~= num_files)
   s = strcat('Creating epochs, baseline correction + occular correction for file...' ,nameArray{count,1});
   disp(s)
   
   % fear:
   Data.(nameArray{count,1}).artefact_rem.fear = ...
       pop_epoch( Data.(nameArray{count,1}).artefact_rem, {  'S  1'  }, [-0.2         0.4], 'epochinfo', 'yes');
   % remove baseline from eeg data:
   Data.(nameArray{count,1}).artefact_rem.fear = ... 
       pop_rmbase(Data.(nameArray{count,1}).artefact_rem.fear, [-200 0] ,[]);
   % get eye-movement electrode:
   HEOG_fear = Data.(nameArray{count,1}).artefact_rem.fear.data(eog_index.(nameArray{count,1}),:,:); 
   HEOG_fear = squeeze(HEOG_fear);
   % perform gratton algorithm to remove blinking artefacts:
   Data.(nameArray{count,1}).artefact_rem.fear = ... 
       gratton(Data.(nameArray{count,1}).artefact_rem.fear.data, HEOG_fear, ...
       blinkcritvolt, blinkcritwin ); 

   % repeat for all emotions: 
   
   % anger:
   Data.(nameArray{count,1}).artefact_rem.anger =... 
       pop_epoch( Data.(nameArray{count,1}).artefact_rem, {  'S  2'  }, [-0.2         0.4], 'epochinfo', 'yes');
   Data.(nameArray{count,1}).artefact_rem.anger = ...
       pop_rmbase(Data.(nameArray{count,1}).artefact_rem.anger, [-200 0] ,[]);
   HEOG_anger = Data.(nameArray{count,1}).artefact_rem.anger.data(eog_index.(nameArray{count,1}),:,:);
   HEOG_anger = squeeze(HEOG_anger);
   Data.(nameArray{count,1}).artefact_rem.anger = ...
       gratton(Data.(nameArray{count,1}).artefact_rem.anger.data, HEOG_anger, ...
       blinkcritvolt, blinkcritwin );
   
   % sadness: 
   Data.(nameArray{count,1}).artefact_rem.sad = ...
       pop_epoch(Data.(nameArray{count,1}).artefact_rem, {  'S  3'  }, [-0.2         0.4], 'epochinfo', 'yes');
   Data.(nameArray{count,1}).artefact_rem.sad = ...
       pop_rmbase(Data.(nameArray{count,1}).artefact_rem.sad, [-200 0] ,[]);
   HEOG_sad = Data.(nameArray{count,1}).artefact_rem.sad.data(eog_index.(nameArray{count,1}),:,:);
   HEOG_sad = squeeze(HEOG_sad);
   Data.(nameArray{count,1}).artefact_rem.sad = ...
       gratton(Data.(nameArray{count,1}).artefact_rem.sad.data, HEOG_sad, ...
       blinkcritvolt, blinkcritwin );

   
   % happiness:
   Data.(nameArray{count,1}).artefact_rem.happy = ...
       pop_epoch(Data.(nameArray{count,1}).artefact_rem, {  'S  4'  }, [-0.2         0.4], 'epochinfo', 'yes');
   Data.(nameArray{count,1}).artefact_rem.happy = ...
       pop_rmbase(Data.(nameArray{count,1}).artefact_rem.happy, [-200 0] ,[]);
   HEOG_happy = Data.(nameArray{count,1}).artefact_rem.happy.data(eog_index.(nameArray{count,1}),:,:);
   HEOG_happy = squeeze(HEOG_happy);
   Data.(nameArray{count,1}).artefact_rem.happy = ...
       gratton(Data.(nameArray{count,1}).artefact_rem.happy.data, HEOG_happy, ...
       blinkcritvolt, blinkcritwin );
   
   % surprise
   Data.(nameArray{count,1}).artefact_rem.surprise = ...
       pop_epoch(Data.(nameArray{count,1}).artefact_rem, {  'S  5'  }, [-0.2         0.4], 'epochinfo', 'yes');
   Data.(nameArray{count,1}).artefact_rem.surprise = ...
       pop_rmbase(Data.(nameArray{count,1}).artefact_rem.surprise, [-200 0] ,[]);
   HEOG_surprise = Data.(nameArray{count,1}).artefact_rem.surprise.data(eog_index.(nameArray{count,1}),:,:);
   HEOG_surprise = squeeze(HEOG_surprise);
   Data.(nameArray{count,1}).artefact_rem.surprise = ...
       gratton(Data.(nameArray{count,1}).artefact_rem.surprise.data, HEOG_surprise, ...
       blinkcritvolt, blinkcritwin );
   
   % neutral 
   Data.(nameArray{count,1}).artefact_rem.neutral = ...
       pop_epoch(Data.(nameArray{count,1}).artefact_rem, {  'S  6'  }, [-0.2         0.4], 'epochinfo', 'yes');
   Data.(nameArray{count,1}).artefact_rem.neutral = ...
       pop_rmbase(Data.(nameArray{count,1}).artefact_rem.neutral, [-200 0] ,[]);
   HEOG_neutral = Data.(nameArray{count,1}).artefact_rem.neutral.data(eog_index.(nameArray{count,1}),:,:);
   HEOG_neutral = squeeze(HEOG_neutral);
   Data.(nameArray{count,1}).artefact_rem.neutral = ...
       gratton(Data.(nameArray{count,1}).artefact_rem.neutral.data, HEOG_neutral, ...
       blinkcritvolt, blinkcritwin );
   
   % fear ambiguous
   Data.(nameArray{count,1}).artefact_rem.fear_amb = ...
       pop_epoch( Data.(nameArray{count,1}).artefact_rem, {  'S  7'  }, [-0.2         0.4], 'epochinfo', 'yes');
   Data.(nameArray{count,1}).artefact_rem.fear_amb = ...
       pop_rmbase(Data.(nameArray{count,1}).artefact_rem.fear_amb, [-200 0] ,[]);
   HEOG_fear_amb = Data.(nameArray{count,1}).artefact_rem.fear_amb.data(eog_index.(nameArray{count,1}),:,:);
   HEOG_fear_amb = squeeze(HEOG_fear_amb);
   Data.(nameArray{count,1}).artefact_rem.fear_amb = ...
       gratton(Data.(nameArray{count,1}).artefact_rem.fear_amb.data, HEOG_fear_amb, ...
       blinkcritvolt, blinkcritwin );
   
   % anger ambiguous
   Data.(nameArray{count,1}).artefact_rem.anger_amb = ...
       pop_epoch(Data.(nameArray{count,1}).artefact_rem, {  'S  8'  }, [-0.2         0.4], 'epochinfo', 'yes');
   Data.(nameArray{count,1}).artefact_rem.anger_amb = ...
       pop_rmbase(Data.(nameArray{count,1}).artefact_rem.anger_amb, [-200 0] ,[]);
   HEOG_anger_amb = Data.(nameArray{count,1}).artefact_rem.anger_amb.data(eog_index.(nameArray{count,1}),:,:);
   HEOG_anger_amb = squeeze(HEOG_anger_amb);
   Data.(nameArray{count,1}).artefact_rem.anger_amb = ...
       gratton(Data.(nameArray{count,1}).artefact_rem.anger_amb.data, HEOG_anger_amb, ...
       blinkcritvolt, blinkcritwin );
   
   % sadness ambiguous
   Data.(nameArray{count,1}).artefact_rem.sad_amb = ...
       pop_epoch( Data.(nameArray{count,1}).artefact_rem, {  'S  9'  }, [-0.2         0.4], 'epochinfo', 'yes');
   Data.(nameArray{count,1}).artefact_rem.sad_amb = ...
       pop_rmbase(Data.(nameArray{count,1}).artefact_rem.sad_amb, [-200 0] ,[]);
   HEOG_sad_amb = Data.(nameArray{count,1}).artefact_rem.sad_amb.data(eog_index.(nameArray{count,1}),:,:);
   HEOG_sad_amb = squeeze(HEOG_sad_amb);
   Data.(nameArray{count,1}).artefact_rem.sad_amb = ...
       gratton(Data.(nameArray{count,1}).artefact_rem.sad_amb.data, HEOG_sad_amb, ...
       blinkcritvolt, blinkcritwin );
   
   % happiness amibguous
   Data.(nameArray{count,1}).artefact_rem.happy_amb = ...
       pop_epoch( Data.(nameArray{count,1}).artefact_rem, {  'S 10'  }, [-0.2         0.4], 'epochinfo', 'yes');
   Data.(nameArray{count,1}).artefact_rem.happy_amb = ...
       pop_rmbase(Data.(nameArray{count,1}).artefact_rem.happy_amb, [-200 0] ,[]);
   HEOG_happy_amb = Data.(nameArray{count,1}).artefact_rem.happy_amb.data(eog_index.(nameArray{count,1}),:,:);
   HEOG_happy_amb = squeeze(HEOG_happy_amb);
   Data.(nameArray{count,1}).artefact_rem.happy_amb = ...
       gratton(Data.(nameArray{count,1}).artefact_rem.happy_amb.data, HEOG_happy_amb, ...
       blinkcritvolt, blinkcritwin );
   
   % surprise ambiguous
   Data.(nameArray{count,1}).artefact_rem.surprise_amb = ...
       pop_epoch( Data.(nameArray{count,1}).artefact_rem, {  'S 11'  }, [-0.2         0.4], 'epochinfo', 'yes');
   Data.(nameArray{count,1}).artefact_rem.surprise_amb = ...
       pop_rmbase(Data.(nameArray{count,1}).artefact_rem.surprise_amb, [-200 0] ,[]);
   HEOG_surprise_amb = Data.(nameArray{count,1}).artefact_rem.surprise_amb.data(eog_index.(nameArray{count,1}),:,:);
   HEOG_surprise_amb = squeeze(HEOG_surprise_amb);
   Data.(nameArray{count,1}).artefact_rem.surprise_amb = ...
       gratton(Data.(nameArray{count,1}).artefact_rem.surprise_amb.data, HEOG_surprise_amb,...
       blinkcritvolt, blinkcritwin );
   
   % neutral ambiguous
   Data.(nameArray{count,1}).artefact_rem.neutral_amb = ...
       pop_epoch( Data.(nameArray{count,1}).artefact_rem, {  'S 12'  }, [-0.2         0.4], 'epochinfo', 'yes');
   Data.(nameArray{count,1}).artefact_rem.neutral_amb = ...
       pop_rmbase(Data.(nameArray{count,1}).artefact_rem.neutral_amb, [-200 0] ,[]);
   HEOG_neutral_amb = Data.(nameArray{count,1}).artefact_rem.neutral_amb.data(eog_index.(nameArray{count,1}),:,:);
   HEOG_neutral_amb = squeeze(HEOG_neutral_amb);
   Data.(nameArray{count,1}).artefact_rem.neutral_amb = ...
       gratton(Data.(nameArray{count,1}).artefact_rem.neutral_amb.data, HEOG_neutral_amb, ...
       blinkcritvolt, blinkcritwin );
   
   
   count = count +1;
    
end


%% Step 8: Calculate Average
%we are summing over the 3rd dimension of the
%Data.(nameArray{count,1}).artefact_rem.fear field, so: 
dim = 3;

count = 1;
while(count ~= num_files)
    % using eeglab function sum(data) or mean(data) which is better?
    % sum and mean both are matlab fucntion and i think mean is better in this case because we have to calculate average  
    ERP.(nameArray{count,1}).fear_sum = sum(Data.(nameArray{count,1}).artefact_rem.fear, dim) /...
        size(Data.eeg1.artefact_rem.fear, dim); % divided by the number of trials
    ERP.(nameArray{count,1}).anger_sum = sum(Data.(nameArray{count,1}).artefact_rem.anger, dim)/...
        size(Data.eeg1.artefact_rem.anger, dim);
    ERP.(nameArray{count,1}).sad_sum = sum(Data.(nameArray{count,1}).artefact_rem.sad, dim) /...
        size(Data.eeg1.artefact_rem.sad, dim);
    ERP.(nameArray{count,1}).happy_sum = sum(Data.(nameArray{count,1}).artefact_rem.happy, dim) /...
        size(Data.eeg1.artefact_rem.happy, dim);
    ERP.(nameArray{count,1}).surprise_sum = sum(Data.(nameArray{count,1}).artefact_rem.surprise, dim) /...
        size(Data.eeg1.artefact_rem.surprise, dim);
    ERP.(nameArray{count,1}).neutral_sum = sum(Data.(nameArray{count,1}).artefact_rem.neutral, dim) /...
        size(Data.eeg1.artefact_rem.neutral, dim);
    ERP.(nameArray{count,1}).fear_amb_sum = sum(Data.(nameArray{count,1}).artefact_rem.fear_amb, dim) /...
        size(Data.eeg1.artefact_rem.fear_amb, dim);
    ERP.(nameArray{count,1}).anger_amb_sum = sum(Data.(nameArray{count,1}).artefact_rem.anger_amb, dim) /...
        size(Data.eeg1.artefact_rem.anger_amb, dim);
    ERP.(nameArray{count,1}).sad_amb_sum = sum(Data.(nameArray{count,1}).artefact_rem.sad_amb, dim) /...
        size(Data.eeg1.artefact_rem.sad_amb, dim);
    ERP.(nameArray{count,1}).happy_amb_sum = sum(Data.(nameArray{count,1}).artefact_rem.happy_amb, dim) /...
        size(Data.eeg1.artefact_rem.happy_amb, dim);
    ERP.(nameArray{count,1}).surprise_amb_sum = sum(Data.(nameArray{count,1}).artefact_rem.surprise_amb, dim) /...
        size(Data.eeg1.artefact_rem.surprise_amb, dim);
    ERP.(nameArray{count,1}).neutral_amb_sum = sum(Data.(nameArray{count,1}).artefact_rem.neutral_amb, dim) /...
        size(Data.eeg1.artefact_rem.neutral_amb, dim);
    count = count + 1;
end

%% visualizing example ERP:

% to find which channel is which number: look at Data.eeg1.chanloc, the row
% corresponds to the channel that is is the labels field

stimulus = 200; % absolute time when stimulus was presented:

% visualizing fear ERP:

% if you would want to plot ALL channels: 
% i= 1;
% while (i~=61)
%     plot(ERP.eeg1.fear_sum(i,:))
%     hold on
%     i = i + 1;
% end

% this example just plots P3 and P4:
plot(ERP.eeg1.fear_sum(12,:)) % P3 channel
title('ERP of fear for one person')
hold on % causes the lines to be in the same plot
plot(ERP.eeg1.fear_sum(27,:)) % P4 Channel
hold off

legend('P3', 'P4')
ylim 'manual'
ylim([-10 10])
line([stimulus stimulus], get(gca, 'ylim')) % draws a vertical line when stimulus is presented
title('ERP of fear for one person')
xlabel('time in ms')
ylabel('voltage') % positive is above x-axis, negative is below x-axis
hold off

% happiness: 

% again, if you want to plot ALL channels: 
% i= 1;
% while (i~=61)
%     plot(ERP.eeg1.happy_sum(i,:))
%     hold on
%     i = i + 1;
% end

plot(ERP.eeg1.happy_sum(12,:)) % P3
title('ERP of happy for one person')
hold on
plot(ERP.eeg1.happy_sum(17,:)) % P4
hold off

legend('P3', 'P4')
ylim 'manual'
ylim([-10 10])
line([stimulus stimulus], get(gca, 'ylim'))
title('ERP of happy for one person')
xlabel('time in ms')
ylabel('voltage')
hold off

%% Step 9: Calculating New Reference:
% when reference channels are not mentioned, then it just takes the average
% as reference, excludes the HEOG channels, which is channel 61

count = 1;
while(count ~= num_files)
    ERP.(nameArray{count,1}).fear_ref = reref(ERP.(nameArray{count,1}).fear_sum, [],...
        'exclude', eog_index.(nameArray{count,1}));
    ERP.(nameArray{count,1}).anger_ref = reref(ERP.(nameArray{count,1}).anger_sum, [], ...
        'exclude', eog_index.(nameArray{count,1}));
    ERP.(nameArray{count,1}).sad_ref = reref(ERP.(nameArray{count,1}).sad_sum, [],...
        'exclude', eog_index.(nameArray{count,1}));
    ERP.(nameArray{count,1}).happy_ref = reref(ERP.(nameArray{count,1}).happy_sum, [],...
        'exclude', eog_index.(nameArray{count,1}));
    ERP.(nameArray{count,1}).surprise_ref = reref(ERP.(nameArray{count,1}).surprise_sum, [],...
        'exclude', eog_index.(nameArray{count,1}));
    ERP.(nameArray{count,1}).neutral_ref = reref(ERP.(nameArray{count,1}).neutral_sum, [],...
        'exclude', eog_index.(nameArray{count,1}));
    ERP.(nameArray{count,1}).fear_amb_ref = reref(ERP.(nameArray{count,1}).fear_amb_sum, [],...
        'exclude', eog_index.(nameArray{count,1}));
    ERP.(nameArray{count,1}).anger_amb_ref = reref(ERP.(nameArray{count,1}).anger_amb_sum, [],...
        'exclude', eog_index.(nameArray{count,1}));
    ERP.(nameArray{count,1}).sad_amb_ref = reref(ERP.(nameArray{count,1}).sad_amb_sum, [],...
        'exclude', eog_index.(nameArray{count,1}));
    ERP.(nameArray{count,1}).happy_amb_ref = reref(ERP.(nameArray{count,1}).happy_amb_sum, [],...
        'exclude', eog_index.(nameArray{count,1}));
    ERP.(nameArray{count,1}).surprise_amb_ref = reref(ERP.(nameArray{count,1}).surprise_amb_sum, [], ...
        'exclude', eog_index.(nameArray{count,1}));
    ERP.(nameArray{count,1}).neutral_amb_ref = reref(ERP.(nameArray{count,1}).neutral_amb_sum, [],...
        'exclude', eog_index.(nameArray{count,1}));
    count = count + 1;
end

%% comparing the signals: before referencing again and after referencing
% this is why it's important to always name new fields to store the data
% because now you can see how the transformation impacted the signal
plot(ERP.eeg1.fear_ref(1,:))
ylim 'manual'
ylim([-10 10])
line([stimulus stimulus], get(gca, 'ylim'))
hold on
plot(ERP.eeg1.fear_sum(1,:))
legend('referenced', 'unreferenced')
title('Comparing unreferenced and referenced averaged fear ERP')
hold off

%% Step 10: Grand Average: 
count = 2; % start at 2, because we initiate the variables with the first subject!
num_subj = num_files - 1;
% initiating variables for sum:
ERP.fear_avg = ERP.eeg1.fear_ref;
ERP.anger_avg = ERP.eeg1.anger_ref;
ERP.sad_avg = ERP.eeg1.sad_ref;
ERP.happy_avg = ERP.eeg1.happy_ref;
ERP.surprise_avg = ERP.eeg1.surprise_ref;
ERP.neutral_avg = ERP.eeg1.neutral_ref;
ERP.fear_amb_avg = ERP.eeg1.fear_amb_ref;
ERP.anger_amb_avg = ERP.eeg1.anger_amb_ref;
ERP.sad_amb_avg = ERP.eeg1.sad_amb_ref;
ERP.happy_amb_avg = ERP.eeg1.happy_amb_ref;
ERP.surprise_amb_avg = ERP.eeg1.surprise_amb_ref;
ERP.neutral_amb_avg = ERP.eeg1.neutral_amb_ref;

while (count ~= num_files)

    ERP.fear_avg = ERP.fear_avg + ERP.(nameArray{count,1}).fear_ref;
    ERP.fear_avg = ERP.fear_avg / num_subj;
    ERP.anger_avg = ERP.anger_avg + ERP.(nameArray{count,1}).anger_ref;
    ERP.anger_avg = ERP.anger_avg / num_subj; 
    ERP.sad_avg = ERP.sad_avg + ERP.(nameArray{count,1}).sad_ref;
    ERP.sad_avg = ERP.sad_avg / num_subj; 
    ERP.happy_avg = ERP.happy_avg + ERP.(nameArray{count,1}).happy_ref;
    ERP.happy_avg = ERP.happy_avg / num_subj;
    ERP.surprise_avg = ERP.surprise_avg + ERP.(nameArray{count,1}).surprise_ref;
    ERP.surprise_avg = ERP.surprise_avg / num_subj;
    ERP.neutral_avg = ERP.neutral_avg + ERP.(nameArray{count,1}).neutral_ref;
    ERP.neutral_avg = ERP.neutral_avg / num_subj;
    
    ERP.fear_amb_avg = ERP.fear_amb_avg + ERP.(nameArray{count,1}).fear_amb_ref;
    ERP.fear_amb_avg = ERP.fear_amb_avg / num_subj;
    ERP.anger_amb_avg = ERP.anger_amb_avg + ERP.(nameArray{count,1}).anger_amb_ref;
    ERP.anger_amb_avg = ERP.anger_amb_avg / num_subj;
    ERP.sad_amb_avg = ERP.sad_amb_avg + ERP.(nameArray{count,1}).sad_amb_ref;
    ERP.sad_amb_avg = ERP.sad_amb_avg / num_subj;
    ERP.happy_amb_avg = ERP.happy_amb_avg + ERP.(nameArray{count,1}).happy_amb_ref;
    ERP.happy_amb_avg = ERP.happy_amb_avg / num_subj;
    ERP.surprise_amb_avg = ERP.surprise_amb_avg + ERP.(nameArray{count,1}).surprise_amb_ref;
    ERP.surprise_amb_avg = ERP.surprise_amb_avg / num_subj;
    ERP.neutral_amb_avg = ERP.neutral_amb_avg + ERP.(nameArray{count,1}).neutral_amb_ref;
    ERP.neutral_amb_avg = ERP.neutral_amb_avg / num_subj;
    
    count = count + 1;
end

%% Looking at ERPs graphically: FEAR
% comparing different locations:
% P3/4
plot(ERP.fear_avg(12,:))
ylim 'manual'
ylim([-10 10])
hold on
plot(ERP.fear_avg(17,:))
line([stimulus stimulus], get(gca, 'ylim'))
legend('P3', 'P4')
title('Grand Average: fear')
xlabel('time in ms')
ylabel('voltage')
saveas(gcf,'fear_P3_P4.png')
hold off

% P7/8
plot(ERP.fear_avg(15,:))
ylim 'manual'
ylim([-10 10])
hold on
plot(ERP.fear_avg(20,:))
line([stimulus stimulus], get(gca, 'ylim'))
legend('P7', 'P8')
title('Grand Average: fear')
xlabel('time in ms')
ylabel('voltage')
saveas(gcf,'fear_P7_P8.png')
hold off

% O1/2
plot(ERP.fear_avg(16,:))
ylim 'manual'
ylim([-10 10])
hold on
plot(ERP.fear_avg(18,:))
line([stimulus stimulus], get(gca, 'ylim'))
legend('O1', 'O2')
title('Grand Average: fear')
xlabel('time in ms')
ylabel('voltage')
saveas(gcf,'fear_O1_O2.png')
hold off

%% Looking at ERPs graphically: HAPPINESS
% P3/4
plot(ERP.happy_avg(12,:))
ylim 'manual'
ylim([-10 10])
hold on
plot(ERP.happy_avg(17,:))
line([stimulus stimulus], get(gca, 'ylim'))
legend('P3', 'P4')
title('Grand Average: happy')
xlabel('time in ms')
ylabel('voltage')
saveas(gcf,'happy_P3_P4.png')
hold off

% P7/8
plot(ERP.happy_avg(15,:))
ylim 'manual'
ylim([-10 10])
hold on
plot(ERP.happy_avg(20,:))
line([stimulus stimulus], get(gca, 'ylim'))
legend('P7', 'P8')
title('Grand Average: happy')
xlabel('time in ms')
ylabel('voltage')
saveas(gcf,'happy_P7_P8.png')
hold off

% O1/2
plot(ERP.happy_avg(16,:))
ylim 'manual'
ylim([-10 10])
hold on
plot(ERP.happy_avg(18,:))
line([stimulus stimulus], get(gca, 'ylim'))
legend('O1', 'O2')
title('Grand Average: happy')
xlabel('time in ms')
ylabel('voltage')
saveas(gcf,'happy_O1_O2.png')
hold off

%% Looking at ERPs graphically: NEUTRAL
% P3/4
plot(ERP.neutral_avg(12,:))
ylim 'manual'
ylim([-10 10])
hold on
plot(ERP.neutral_avg(17,:))
line([stimulus stimulus], get(gca, 'ylim'))
legend('P3', 'P4')
title('Grand Average: neutral')
xlabel('time in ms')
ylabel('voltage')
saveas(gcf,'neutral_P3_P4.png')
hold off

% P7/8
plot(ERP.neutral_avg(15,:))
ylim 'manual'
ylim([-10 10])
hold on
plot(ERP.neutral_avg(20,:))
line([stimulus stimulus], get(gca, 'ylim'))
legend('P7', 'P8')
title('Grand Average: neutral')
xlabel('time in ms')
ylabel('voltage')
saveas(gcf,'neutral_P7_P8.png')
hold off

% O1/2
plot(ERP.neutral_avg(16,:))
ylim 'manual'
ylim([-10 10])
hold on
plot(ERP.neutral_avg(18,:))
line([stimulus stimulus], get(gca, 'ylim'))
legend('O1', 'O2')
title('Grand Average: neutral')
xlabel('time in ms')
ylabel('voltage')
saveas(gcf,'neutral_O1_O2.png')
hold off
