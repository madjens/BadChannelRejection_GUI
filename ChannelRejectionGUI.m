function [badchannels_idx,segment_rating,plot_data] = ChannelRejectionGUI(data_eeg,fs,badchannels,segment_rating,chanlocs,figure_name)
% input
% data_eeg       :   detrended eeg data (samples x channels)
% fs             :   sampling frequency of eeg data (scalar)
% badchannels    :   pass if badchannels have already been found, other otherswise empty (vector of channel nummbers)
% segment_rating :   rating how the quality of the recording, if it has already been rated before, otherswise leave empty (scalar)
% chanlocs       :   location file in mat format as used by EEG lab. 
% figure_name    :   name of the channel rejecton GUI figure

% output          
% badchannels_idx : the channel number of bad channels as selected in the GUI
% segment_rating  : the rating as selected in the GUI
% plot_data       : the eeg data as plotted in the GUI (data is resampled to 128 Hz if the sampling frequency is higher, for display purposes)
%
% jenma/jmad 2019

if nargin < 1
    duration = 120; % seconds
    fs = 512;
    samples = duration*fs;
    
    %load example of channel locations
    options.location_file = './BioSemi64.mat'; % only used for testing purposes
    load(options.location_file,'chanlocs')
    channels = length(chanlocs);
    
    %generate random data for test purposes
    variance = 20;
    data_eeg = randn(samples,channels)*variance;
    
    figure_name = '';
    badchannels = [];
    segment_rating = [];
end

%% read channel information
channel_no = (1:size(chanlocs,2))';
channel_names = {chanlocs.labels}';
channel_plotnames = cellfun(@(c,d) ['(' c ') ' num2str(d)],channel_names,num2cell(channel_no),'UniformOutput',false);

fprintf('Performing manual channel rejection...')

%% create the figure
figureHandle = figure('WindowButtonDownFcn',@myFigureCallback,'Visible','on','Name',figure_name,'units', 'normalized', ...
    'OuterPosition',[0 0.04 1 0.95],'Toolbar','none', 'Resize','on','MenuBar','none', 'color', [1 1 1],'DeleteFcn',@myDeleteFcn,'CloseRequestFcn',@closeFcn,'KeyPressFcn',@myKeyPressFcn);

%% handle ratings of the segment
if ~exist('segment_rating','var')
    segment_rating = 1;
elseif isnan(segment_rating)
    segment_rating = 1;
elseif isempty(segment_rating)
    segment_rating = 1;
else
    segment_rating = segment_rating+1;
end

%% load all the data into the figure userdata to see betweeen callbacks

%resample data for display purposes
if fs>128
    data_eeg = resample(data_eeg,128,fs);
    fs = 128;
end

orig_data = data_eeg(1:(floor(size(data_eeg,1)/2)*2),:); clear eeg; %only use an even amount of samples
plot_data = orig_data;
Nsamples = size(orig_data,1);
Nchannels = size(orig_data,2);
Nbadchannels = length(badchannels);

if isnumeric(badchannels)
    badchannels_idx = badchannels;
else
    badchannels_idx = find(arrayfun(@(s) find(ismember(s.labels,channel_names)),badchannels));
end
updateBadChannel;

maxSample = Nsamples-1;
minSample = 0;
maxTime = maxSample/fs;
minTime = minSample/fs;
maxVal = max(max(orig_data(:)));
minVal = min(min(orig_data(:)));

TimeAxis = (minSample:maxSample)./fs;

plotLine_upper = [];
plotLine_lower = [];
CurrentChannel = [];

%this is the axis offset used to plot the EEG
axis_diff = 50;
channel_offset = 0:axis_diff:axis_diff*(Nchannels-1);

caxisData = [-axis_diff axis_diff];

%pre define the histograms
histogram_interval = linspace(minVal,maxVal,50);
histograms = nan(Nchannels,length(histogram_interval));

%pre define the spectra for each channel
spectra_Y = nan(Nchannels,Nsamples);
ylim_spectrum = [-80, 10];

%pre define the power in dB
dbPower = nan(Nchannels,1);

%these are the ratings of the signal quality that can be used by the user
ratings = {'','1','2','3','4','5'};

displayType = 'Time series';

%% pre compute histograms and spectra
for iChannel = 1:Nchannels
    histograms(iChannel,:) = hist(plot_data(:,iChannel),histogram_interval);
    spectra_Y(iChannel,:) = fft(plot_data(:,iChannel),Nsamples)/Nsamples;
    dbPower(iChannel) = db(std(plot_data(:,iChannel)));
end

spectra_f = linspace(0,fs/2,Nsamples/2+1);

%% create some uicontrols
pbh = uicontrol(figureHandle,'Style', 'pushbutton','Tag','pushbutton', 'String', 'Select channel','units', 'normalized','FontSize',15,'Position', [0.27 0.92 0.07 0.07], 'Callback',@mySelectbuttonCallback);
pbh_save = uicontrol(figureHandle,'Style', 'pushbutton','Tag','savebutton', 'String', 'Exit','units', 'normalized','FontSize',15,'Position', [0.35 0.92 0.04 0.07], 'Callback',@myExitbuttonCallback);
pbh_zero = uicontrol(figureHandle,'Style', 'pushbutton','Tag','zerobutton', 'String', 'Zero','units', 'normalized','FontSize',15,'Position', [0.27 0.83 0.03 0.07], 'Callback',@myZerobuttonCallback);
pbh_interp = uicontrol(figureHandle,'Style', 'pushbutton','Tag','interpbutton', 'String', 'Interp','units', 'normalized','FontSize',15,'Position', [0.31 0.83 0.03 0.07], 'Callback',@myInterpbuttonCallback);
pbh_reset = uicontrol(figureHandle,'Style', 'pushbutton','Tag','resetbutton', 'String', 'Reset','units', 'normalized','FontSize',15,'Position', [0.35 0.83 0.04 0.07], 'Callback',@myResetbuttonCallback);

%editbox where user can specify the distance between the eeg timeseries
text_scale = uicontrol(figureHandle,'Style', 'text','Tag','textaxisdiff', 'String', 'Scale','units', 'normalized','FontSize',15,'Position', [0.41 0.95 0.04 0.02],'BackGroundColor',[1 1 1]);
edit_scale = uicontrol(figureHandle,'Style', 'edit','Tag','axisdiff', 'String', num2str(axis_diff),'units', 'normalized','FontSize',15,'Position', [0.41 0.92 0.04 0.03], 'Callback',@myAxisdiffCallback);
pbh_axis_up = uicontrol(figureHandle,'Style', 'pushbutton','Tag','axis_up', 'String', '^','units', 'normalized','FontSize',15,'Position', [0.455 0.935 0.01 0.02], 'Callback',@myAxisCallback);
pbh_axis_down = uicontrol(figureHandle,'Style', 'pushbutton','Tag','axis_up', 'String', 'v','units', 'normalized','FontSize',15,'Position', [0.455 0.915 0.01 0.02], 'Callback',@myAxisCallback);

%this is the dropdown menu that the user can rate the eeg data
text_rating = uicontrol(figureHandle,'Style', 'text','Tag','textscale', 'String', 'Rating','units', 'normalized','FontSize',15,'Position', [0.415 0.86 0.04 0.02],'BackGroundColor',[1 1 1]);
dropDownHandle = uicontrol(figureHandle,'Style', 'popupmenu','Tag','segment_rating','Value',segment_rating, 'String', ratings,'units', 'normalized','HorizontalAlignment','center','FontSize',15,'BackgroundColor', [1 1 1],'Position', [0.41 0.83 0.05 0.03], 'Callback',@mySegmentRatingDropdownCallback);

% bad channel listbox, this shows the list of bad channels that the user has selected
text_list = uicontrol(figureHandle,'Style', 'text','Tag','textscale', 'String', 'Bad channels','units', 'normalized','FontSize',15,'Position', [0.455 0.96 0.08 0.02],'BackGroundColor',[1 1 1]);
listboxHandle = uicontrol(figureHandle,'Style', 'listbox','Tag','listbox', 'units', 'normalized','HorizontalAlignment','center','FontSize',15,'Position', [0.47 0.83 0.05 0.13],'String',badchannels_idx);
text_badchannel = uicontrol(figureHandle,'Style', 'text','Tag','textscale', 'String', sprintf('%d/%d (%1.0f%%)',Nbadchannels,Nchannels,Nbadchannels/Nchannels*100),'units', 'normalized','FontSize',15,'Position', [0.455 0.795 0.08 0.025],'BackGroundColor',[1 1 1]);

%
text_displaytype = uicontrol(figureHandle,'Style', 'text','Tag','textscale', 'String', 'Display type','units', 'normalized','FontSize',15,'Position', [0.53 0.93 0.08 0.02],'BackGroundColor',[1 1 1]);
bg = uibuttongroup('Visible','on','Position',[0.53 0.85 0.09 0.08],'SelectionChangedFcn',@myRadiobuttonCallback,'BackGroundColor',[1 1 1]);
r1 = uicontrol(bg,'Style','radiobutton','units', 'normalized','String','Time series','Position',[0.05 0.6 0.8 0.25], 'HandleVisibility','off','FontSize',15,'BackGroundColor',[1 1 1]);
r2 = uicontrol(bg,'Style','radiobutton','units', 'normalized','String','Image','Position',      [0.05 0.1 0.8 0.25], 'HandleVisibility','off','FontSize',15,'BackGroundColor',[1 1 1]);

%% create some axis
axesHandle_EEG = axes('Position',[0.04 0.035 0.95 0.76],'Tag','EEG data','HitTest','off');drawnow
axesHandle_Hist = axes('Position',[0.04 0.83 0.2 0.15],'Tag','Histogram','HitTest','off', 'Xtick',[],'Ytick',[],'Title','Channel histogram');drawnow
axesHandle_Spectrum = axes('Position',[0.65 0.83 0.3 0.15],'Tag','Spectrum','HitTest','off', 'Xtick',[],'Ytick',[],'Title','Spectrum');drawnow

%% plot eeg on the image
plotEEGdata;

xlabel(axesHandle_Spectrum,'Frequency [Hz]')
ylim(axesHandle_Spectrum,ylim_spectrum)
waitfor(figureHandle) %this pauses the command prompt until the figure has been closed

% handles the checkbox if the segment is 'bad' or 'good'
    function mySegmentRatingDropdownCallback(~,~)
        segment_rating = get(dropDownHandle,'Value');
    end

%% this function handles all keypressed
    function myKeyPressFcn(figureHandle,eventdata)
        switch get(figureHandle,'CurrentKey')
            case 'uparrow'
                if ~isempty(CurrentChannel)
                    CurrentChannel = min(CurrentChannel+1,Nchannels);
                    draw_points;
                end
            case 'downarrow'
                if ~isempty(CurrentChannel)
                    CurrentChannel = max(CurrentChannel-1,1);
                    draw_points;
                end
            case {'s','space'}
                childrenHandles = get(figureHandle,'Children');
                buttonHandle = childrenHandles(cellfun(@(s) contains(s,'pushbutton'), arrayfun(@(s) s.Tag,childrenHandles,'UniformOutput',false)));
                mySelectbuttonCallback(buttonHandle,eventdata)
            case {'delete','del'}
                childrenHandles = get(figureHandle,'Children');
                badchannels_idx(badchannels_idx==CurrentChannel) = [];
                plot_data(:,CurrentChannel) = orig_data(:,CurrentChannel);
                listboxHandle.String = badchannels_idx;
                Nbadchannels = length(badchannels_idx);
                updateBadChannelCounter;
                updateBadChannel;
                plotEEGdata;
            case {'escape','e'}
                myExitbuttonCallback;
            case {'1','2','3','4','5'}
                set(dropDownHandle,'Value',str2double(eventdata.Key)+1);
                mySegmentRatingDropdownCallback;
            case {'i'}
                myInterpbuttonCallback;
        end
    end

    function myDeleteFcn(figureHandle,~)
        delete(figureHandle)
    end
    function closeFcn(figureHandle,~)
        %         badchannels = inf;
        delete(figureHandle)
    end
    function myAxisdiffCallback(fieldHandle,~)
        axis_diff = str2double(get(fieldHandle,'String'));
        if axis_diff<=0
            axis_diff=1;
        end
        channel_offset = 0:axis_diff:axis_diff*(Nchannels-1);
        caxisData = [-axis_diff axis_diff];
        plotEEGdata;
    end

    function myAxisCallback(buttonHandle,eventdata)
        if strcmpi(eventdata.Source.String,'^')
            axis_diff = axis_diff+10;
        elseif strcmpi(eventdata.Source.String,'v')
            axis_diff = axis_diff-10;
        end
        if axis_diff<=0
            axis_diff=1;
        end
        
        set(edit_scale,'String',num2str(axis_diff))
        channel_offset = 0:axis_diff:axis_diff*(Nchannels-1);
        caxisData = [-axis_diff axis_diff];
        plotEEGdata;
    end
    function myResetbuttonCallback(~,~)
        % first get all the handles for the different axis in the figure
        
        %reset the data and badchannels
        plot_data = orig_data;
        badchannels_idx = [];
        Nbadchannels = 0;
        updateBadChannelCounter;
        updateBadChannel;
        
        % re-compute all histograms and spectra
        maxVal = max(max(plot_data(:)));
        minVal = min(min(plot_data(:)));
        histogram_interval = linspace(minVal,maxVal,50);
        
        for iChannelf = 1:Nchannels
            histograms(iChannelf,:) = hist(plot_data(:,iChannelf),histogram_interval);
            spectra_Y(iChannelf,:) = fft(plot_data(:,iChannelf),Nsamples)/Nsamples;
            dbPower(iChannelf) = db(std(plot_data(:,iChannelf)));
        end
        
        %plot the reset eeg data
        plotEEGdata;
        
        %reset the listbox with badchannels
        listboxHandle.String = badchannels_idx;
    end
%% simply zeros out the selected channels and replots the eeg channels
    function myZerobuttonCallback(~,~)
        
        % zero out data of the selected channels
        plot_data(:,badchannels_idx) = 0;
        
        %replot the data
        plotEEGdata;
        
        % recomputes the histograms and spectra when they are now zerod out
        maxVal = max(max(plot_data(:)));
        minVal = min(min(plot_data(:)));
        histogram_interval = linspace(minVal,maxVal,50);
        
        for iChannelt = 1:Nchannels
            histograms(iChannelt,:) = hist(plot_data(:,iChannelt),histogram_interval);
            spectra_Y(iChannelt,:) = fft(plot_data(:,iChannelt),Nsamples)/Nsamples;
        end
    end
    function myInterpbuttonCallback(~,~)
        % first get all the handles for the different axis in the figure
        set(pbh_interp,'Enable','off')
        drawnow;
        
        % interpolate out data of the selected channels
        plot_data = fillBadChannels(plot_data,badchannels_idx,'interp',chanlocs);
%                     fillBadChannels(data,badchannels,method,chanlocs)
        %replot the data
        plotEEGdata;
        
        % recomputes the histograms and spectra when they are now zeroed out
        maxVal = max(max(plot_data(:)));
        minVal = min(min(plot_data(:)));
        histogram_interval = linspace(minVal,maxVal,50);
        
        for iChannels = 1:Nchannels
            histograms(iChannels,:) = hist(plot_data(:,iChannels),histogram_interval);
            spectra_Y(iChannels,:) = fft(plot_data(:,iChannels),Nsamples)/Nsamples;
        end
        set(pbh_interp,'Enable','on')
    end
%% save the data and close the figure
    function myExitbuttonCallback(~,~)
        segment_rating = get(dropDownHandle,'Value')-1;
        if segment_rating~=0
            close(figureHandle)
        end
    end

    function mySelectbuttonCallback(~,~)
        badchannels_idx = unique([badchannels_idx; CurrentChannel]);
        listboxHandle.String = badchannels_idx;
        Nbadchannels = length(badchannels_idx);
        updateBadChannel;
        updateBadChannelCounter;
        plotEEGdata;
    end
%% callback handling all button presses on the figure
    function myFigureCallback(~,~)
        % first get all the handles for the different axis in the figure
        
        % get the coordinates of the buttonpress
        coordinates = get(axesHandle_EEG,'CurrentPoint');
        selectedPoint = round(coordinates(1,2)); %can only select an integer channel
        selectedTime = coordinates(1,1);
        
        if strcmpi(displayType,'Time series')
            [~,selectedChannel] = min(abs(channel_offset-selectedPoint));
        elseif strcmpi(displayType,'Image')
            selectedChannel = selectedPoint;
        end
        % only react if the buttonpress is on the eeg axes
        if (selectedTime>minSample && selectedTime<maxSample && selectedChannel>=1 && selectedChannel<=Nchannels)
            CurrentChannel = selectedChannel;
            draw_points
        end
    end

    function myRadiobuttonCallback(~,event)
        displayType = event.NewValue.String;
        plotEEGdata
    end

%% function handling drawing on the eeg axes, histogram and spectra
    function draw_points
        %% plot the selected channel
        % delete plot line if its already in the figure
        if ~isempty(plotLine_upper)
            delete(plotLine_upper)
        end
        if ~isempty(plotLine_lower)
            delete(plotLine_lower)
        end
        hold(axesHandle_EEG,'on')
        axesHandle_Tag = axesHandle_EEG.Tag;
        
        if strcmpi(displayType,'Time series')
            plotLine_upper = plot(axesHandle_EEG,[minTime maxTime],[channel_offset(CurrentChannel)+axis_diff/2 channel_offset(CurrentChannel)+axis_diff/2],'Color',[0 0 1],'LineWidth',1);
            plotLine_lower = plot(axesHandle_EEG,[minTime maxTime],[channel_offset(CurrentChannel)-axis_diff/2 channel_offset(CurrentChannel)-axis_diff/2],'Color',[0 0 1],'LineWidth',1);
        elseif strcmpi(displayType,'Image')
            plotLine_upper = plot(axesHandle_EEG,[minTime maxTime],[channel_no(CurrentChannel)+1/2 channel_no(CurrentChannel)+1/2],'Color',[0 0 1],'LineWidth',1);
            plotLine_lower = plot(axesHandle_EEG,[minTime maxTime],[channel_no(CurrentChannel)-1/2 channel_no(CurrentChannel)-1/2],'Color',[0 0 1],'LineWidth',1);
        end
        
        set(axesHandle_EEG,'Tag',axesHandle_Tag) %the plot command resets the tag field on the handle, so we reset it (there might be a better way, its a bit hacky)
        drawnow % flush the graphics buffer
        
        %% plot the histogram of the selected channel
        axesHandle_Tag = axesHandle_Hist.Tag;
        cla(axesHandle_Hist,'reset') %clear the axis to get ready for the new plot
        bar(axesHandle_Hist,histogram_interval,histograms(CurrentChannel,:),'k')
        title(axesHandle_Hist,sprintf('Histogram of channel %d (%s) power=%1.3f dB',CurrentChannel,channel_names{CurrentChannel},dbPower(CurrentChannel)))
        set(axesHandle_Hist,'Ygrid','on')
        set(axesHandle_Hist,'Tag',axesHandle_Tag) %the plot command resets the tag field on the handle, so we reset it (there might be a better way, its a bit hacky)
        drawnow % flush the graphics buffer
        
        %% plot the spectrum of the selected channel
        axesHandle_Tag = axesHandle_Spectrum.Tag;
        cla(axesHandle_Spectrum,'reset') %clear the axis to get ready for the new plot
        plot(axesHandle_Spectrum,spectra_f,10*log10(2*abs(spectra_Y(CurrentChannel,1:Nsamples/2+1)).^2));
        title(axesHandle_Spectrum,sprintf('Single-Sided Power Spectrum of channel %d (%s)',CurrentChannel,channel_names{CurrentChannel}))
        xlabel(axesHandle_Spectrum,'Frequency [Hz]')
        ylabel(axesHandle_Spectrum,'10log_1_0|Y(F)|^2')
        xlim(axesHandle_Spectrum,[min(spectra_f),64])
        ylim(axesHandle_Spectrum,ylim_spectrum)
        set(axesHandle_Spectrum,'Ygrid','on')
        set(axesHandle_Spectrum,'Tag',axesHandle_Tag) %the plot command resets the tag field on the handle, so we reset it (there might be a better way, its a bit hacky)
        drawnow % flush the graphics buffer
    end

    function plotEEGdata
        if strcmpi(displayType,'Time series')
            plotTimeseriesEEG;
        elseif strcmpi(displayType,'Image')
            plotImageEEG;
        end
        
    end
    function plotTimeseriesEEG
        cla(axesHandle_EEG,'reset')
        hold(axesHandle_EEG,'on');
        
        for iChannel_ = 1:Nchannels
            if ismember(iChannel_,unique(badchannels_idx))
                plot(axesHandle_EEG,TimeAxis,plot_data(:,iChannel_)+channel_offset(iChannel_),'r'), hold on;
            else
                plot(axesHandle_EEG,TimeAxis,plot_data(:,iChannel_)+channel_offset(iChannel_),'k');
            end
            hold(axesHandle_EEG,'on')
        end
        xlim(axesHandle_EEG,[min(TimeAxis) max(TimeAxis)])
        ylim(axesHandle_EEG,[min(channel_offset)-axis_diff max(channel_offset)+axis_diff])
        set(axesHandle_EEG,'YTick',channel_offset)
        set(axesHandle_EEG,'YTickLabels',channel_plotnames)
        set(axesHandle_EEG,'YMinorTick','off')
        set(axesHandle_EEG,'Tag','EEG data')
        set(axesHandle_EEG,'TickLength',[0 0])
        hold(axesHandle_EEG,'on');
        ylabel(axesHandle_EEG,'Channels')
        xlabel(axesHandle_EEG,'Time [sec]')
    end

    function plotImageEEG
        %         cla(axesHandle_EEG,'reset')
        hold(axesHandle_EEG,'on');
        imagesc(axesHandle_EEG,TimeAxis,channel_no,plot_data','HitTest','off');
        
        for iChannel_ = 1:Nchannels
            if ismember(iChannel_,unique(badchannels_idx))
                badchanLine_upper = plot(axesHandle_EEG,[minTime maxTime],[channel_no(iChannel_)+1/2 channel_no(iChannel_)+1/2],'Color',[1 0 0],'LineWidth',0.9,'LineStyle',':');
                badchanLine_lower = plot(axesHandle_EEG,[minTime maxTime],[channel_no(iChannel_)-1/2 channel_no(iChannel_)-1/2],'Color',[1 0 0],'LineWidth',0.9,'LineStyle',':');
                badchanne_dot_end = plot(axesHandle_EEG,maxTime-0.05,channel_no(iChannel_),'MarkerSize',8,'MarkerEdgeColor',[1 0 0],'LineWidth',0.7,'Marker','s','MarkerFaceColor',[1 0 0]);
                badchanne_dot_start = plot(axesHandle_EEG,minTime+0.05,channel_no(iChannel_),'MarkerSize',8,'MarkerEdgeColor',[1 0 0],'LineWidth',0.7,'Marker','s','MarkerFaceColor',[1 0 0]);
            end
            hold(axesHandle_EEG,'on')
        end
        drawnow
        set(axesHandle_EEG,'Tag','EEG data');
        caxis(axesHandle_EEG,caxisData)
        set(axesHandle_EEG,'YTick',channel_no)
        set(axesHandle_EEG,'YMinorTick','off')
        xlim(axesHandle_EEG,[min(TimeAxis) max(TimeAxis)])
        ylim(axesHandle_EEG,[min(channel_no)-0.5 max(channel_no)+0.5])
        set(axesHandle_EEG,'YTickLabels',channel_plotnames)
        set(axesHandle_EEG,'YMinorTick','off')
        set(axesHandle_EEG,'TickLength',[0 0])
        hold(axesHandle_EEG,'on');
        ylabel(axesHandle_EEG,'Channels')
        xlabel(axesHandle_EEG,'Time [sec]')
    end

    function updateBadChannelCounter
        set(text_badchannel,'String',sprintf('%d/%d (%1.0f%%)',Nbadchannels,Nchannels,Nbadchannels/Nchannels*100))
    end

    function updateBadChannel
        badchannels = chanlocs(badchannels_idx);
        
        for iC = 1:length(badchannels)
            badchannels(iC).status = 'bad';
        end
    end
end