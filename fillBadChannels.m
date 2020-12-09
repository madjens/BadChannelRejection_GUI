function data = fillBadChannels(data,badchannels,method,chanlocs)
% data = fillBadchannels(data,badchannels,method,locationfile)
%
%
% data          : samples x channels
% badchannels   : channel number x 1
% method        : either 'interp' or 'zeros'
% chanlocs      : location file in EEG lab format, only used for interp method (channel number, theta, rho, channel name)
%
% jmad/jenma 2018
%
if nargin < 1
    data = randn(1000,64);
end
if nargin < 2
    badchannels = randi(64,3,1);
end
if nargin < 3
    method = 'interp';
end
if nargin < 4
    load('..\data\location_file\BioSemi64.mat');
end

if iscell(badchannels) && size(badchannels,2)==3
    method = badchannels{2};
    chanlocs = badchannels{3};
    badchannels = badchannels{1};
end

Nchannels = size(data,2);
Nsamples = size(data,1);

if ~isempty(data) && ~isempty(badchannels)
    if strcmpi(method,'interp')
        fprintf('Starting interpolation...')
        
        %% read the location file
        %calcuate the X and Y coordinates i.e. the 3D coordinated projected in to 2D space
        theta = [chanlocs.theta];
        theta_rad = (theta+90)/360*2*pi; % rotate the coordinate system and convert to radian
        rho  = [chanlocs.radius]; %get the lenght of the vector
        [X,Y] = pol2cart(theta_rad,rho); %convert from polar to cartesian coordinate system
        
        goodchannels = setdiff(1:Nchannels,badchannels)';
        
        fprintf('channel(s) '), fprintf('%d ',badchannels)
        
        %good channels coordinates
        Xgc = X(goodchannels)';
        Ygc = Y(goodchannels)';
        
        %bad channels coordinates
        Xbc = X(badchannels)';
        Ybc = Y(badchannels)';
        
        % start interpolating
        for iSample = 1:Nsamples
            %get the eeg sample values for the good channels
            Vgc = data(iSample,goodchannels)';
            
            %create a scatter interpolation function
            F = scatteredInterpolant(Xgc,Ygc,Vgc);
            
            %use the function to interpolate missing values
            data(iSample,badchannels) = F(Xbc,Ybc);
        end
        
    elseif strcmpi(method,'zeros')
        fprintf('Filling with zeros...')
        fprintf('channel(s) '), fprintf('%d ',badchannels)
        data(:,badchannels) = 0;
    end
    fprintf('done\n')
end
