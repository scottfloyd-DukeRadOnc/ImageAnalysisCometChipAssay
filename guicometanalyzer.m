function varargout = guicometanalyzer(varargin)
% GUICOMETANALYZER M-file for guicometanalyzer.fig
%      GUICOMETANALYZER, by itself, creates a new GUICOMETANALYZER or raises the existing
%      singleton*.
%
%      H = GUICOMETANALYZER returns the handle to a new GUICOMETANALYZER or the handle to
%      the existing singleton*.
%
%      GUICOMETANALYZER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUICOMETANALYZER.M with the given input arguments.
%
%      GUICOMETANALYZER('Property','Value',...) creates a new GUICOMETANALYZER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before guicometanalyzer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property
%      application
%      stop.  All inputs are passed to guicometanalyzer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help guicometanalyzer

% Last Modified by GUIDE v2.5 29-Sep-2009 10:28:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @guicometanalyzer_OpeningFcn, ...
                   'gui_OutputFcn',  @guicometanalyzer_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
end

% --- Executes just before guicometanalyzer is made visible.
function guicometanalyzer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to guicometanalyzer (see VARARGIN)

% Choose default command line output for guicometanalyzer
handles.output = hObject;

% Initialize handle objects
handles.SelectedFiles = 1;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes guicometanalyzer wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = guicometanalyzer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

analysisok = mwcometanalyzer(handles);

guidata(hObject,handles)
end

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
end

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
end

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
end

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double
end

% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
end


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

FileDirectory = uigetdir('C:\Documents and Settings\Dave Wood\My Documents\Experiments\');
set(handles.edit6,'String',FileDirectory);
dirlist = dir(fullfile(FileDirectory,'*.tif'));
[filenames,dummyvar] = sortrows({dirlist.name}');
set(handles.listbox1,'String',filenames);
handles.FileDirectory = FileDirectory;
guidata(hObject,handles)
end


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1

handles.SelectedFiles = get(hObject,'Value');
guidata(hObject,handles);
end


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double
end


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

FileNames = get(handles.listbox1,'String');
handles.SelectedFilenames = FileNames(handles.SelectedFiles);
set(handles.listbox1,'String',handles.SelectedFilenames);
set(handles.listbox1,'Value',1:length(handles.SelectedFiles));
guidata(hObject,handles);
end


function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double
end

% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function analysisok = mwcometanalyzer(handles)
% Function mwcometanalyzer
% This function analyzes images of microwell comets that are in a grid.
%
% David K. Wood and Drew Regitsky
% Modified 2/19/09
%
% Inputs:
%
% FileDirectory []  - tells the program where to look for comet images. If no
%   directory is set, the program will allow the user to browse for the
%   appropriate folder.
%
% convfac [0.64]    - sets the pixel to distance conversion factor for the imaging
%   microscope and is used to calculate comet parameters. Units are
%   microns/pixel. Default is for a 10X objective.
%
% objective [10]    - sets the magnification of imaging for writing to
%   output file.
%
% hdthresh [0.15]    - sets the threshold used to detect the beginning of a
%   comet in the analysis routine.
%
% tailthresh [0.02] - sets the threshold used to detect the end of a comet
%   in the analysis routine
%
% headdia [40]      - approximate size of microwells. Units are pixels.

%% Assign parameters from handles structure

tempvar = get(handles.edit2,'String');
hdthresh = str2double(tempvar);

tempvar = get(handles.edit3,'String');
tailthresh = str2double(tempvar);

tempvar = get(handles.edit4,'String');
convfac = str2double(tempvar);

tempvar = get(handles.edit5,'String');
objective = str2double(tempvar);

tempvar = get(handles.edit7,'String');
headdia = str2double(tempvar);

tempvar = get(handles.edit9,'String');
imagerot = str2double(tempvar);

tempvar = get(handles.edit10,'String');
scale_intensity = str2double(tempvar);

tempvar = get(handles.edit12,'String');
PIXEL_SPACING = round(str2double(tempvar)/convfac);    %Spacing between comets in array

tempvar = get(handles.edit11,'String');
CROPOFFSET = round(str2double(tempvar)/convfac);   %X-offset for cropped comet picture to include tail (positive = to the left)

tempvar = get(handles.edit13,'String');
MINCOMETAREA = round(str2double(tempvar)/convfac^2);   %Minimum number of pixels in object for it to be a comet

tempvar = get(handles.edit14,'String');
NearestNeighbors = str2double(tempvar);


FileDirectory = handles.FileDirectory;
filenames = handles.SelectedFilenames;


%% Select files for analysis

    if isempty(FileDirectory) % if no directory is specified open gui
        FileDirectory = uigetdir;
    end
    
    NumFiles = length(filenames); % total number of selected files
    
    %% Analysis Loop
    
    h = waitbar(0,'Analyzing File 0 of 0 (0%)'); % initialize waitbar
    % Begin looping through files
    for fileID = 1:NumFiles
        
        % Set Parameters for finding comets in images
        CROPWIDTH = round(PIXEL_SPACING*0.9);    %Width of cropped comet picture to output
        CROPHEIGHT = round(CROPWIDTH/2);  %Height of cropped comet picture to output

        %Set default values for comet analysis parameters
        if isempty(convfac)
            convfac = 0.64; % 0.64 microns/pixel for 10X objective
        end
        if isempty(objective)
            objective = 10; % 10X objective
        end
        if isempty(hdthresh)
            hdthresh = 0.15; % 15% of peak intensity
        end
        if isempty(tailthresh)
            tailthresh = 0.02; % 2% of peak intensity
        end
        if isempty(headdia)
            headdia=40; % 40 pixels
        end
        
        % File housekeeping
        filename = fullfile(FileDirectory,filenames{fileID}); % assign placeholder variable
        display(filename); % Show user what file is being analyzed
        fileinfo = imfinfo(filename);
        Nimages = length(fileinfo); % how many frames in tif
        clear fileinfo % clear unused variable

        %Open file for comet data and write initial info to it
        fid0 = fopen([filename(1:length(filename)-4),'.txt'],'wt');
        fprintf(fid0,'%s \n',['Size Calibration: ',num2str(convfac),' micron/Px']);
        fprintf(fid0,'%s \n',['Objective: ',num2str(objective),'X']);
        fprintf(fid0,'\n');
        fprintf(fid0,'%s \t %s \t %s \t %s \t %s\n','%Head DNA','%Tail DNA','OTM (um)','Tail Len. (um)','Comet Len. (um)');

        % Initialize number of comets in file
        cometsInFile = 0;
        
        % Update wait bar to include number of frames in file
        waitbar(0,h,['Analyzing File ' num2str(fileID) ' of ' num2str(NumFiles) ' (0%)']);
        % Begin looping through image frames in tif
        for nimage = 1:Nimages
            % Update waitbar at each frame
            waitbar(nimage/Nimages,h,['Analyzing File ',num2str(fileID),' of ',num2str(NumFiles),...
                ' (',num2str(round(nimage/Nimages*100)),'%)']);
            % Read in image from file
            cometimage = imread(filename,nimage)/scale_intensity;
            cometimage = imrotate(cometimage,imagerot); %rotate image
            % Find comet in images using findcomets function
            [temp_objs,NumComets] = findcomets(cometimage,MINCOMETAREA,PIXEL_SPACING,NearestNeighbors);
            % Analyze identified comets using runcometanlaysis function and
            % save data to fid0
            if NumComets == 0
                NumAnalyzed = 0;
            else
            NumAnalyzed = runcometanalysis(cometimage,temp_objs,NumComets,CROPWIDTH,...
                CROPHEIGHT,CROPOFFSET,convfac,hdthresh,tailthresh,headdia,fid0);
            end
            clear cometimage % clear comet image
            cometsInFile = cometsInFile + NumAnalyzed; % update number of comets
        end
        
        % Housekeeping and cleanup
        fclose(fid0); % close data file
        clear fid0 I % clear file I/O and image
        display(['Total number of ANALYZED comets in images: ',num2str(cometsInFile)]);
    end
    close(h) % close waitbar
    analysisok = 'Analysis OK'; % assign dummy output
end

function [ObjectsInfo,NumComets] = findcomets(I,MINCOMETAREA,PIXEL_SPACING,NearestNeighbors)
            
    bwim = im2bw(I,graythresh(I)); % Create bw image using otsu thresholding method
    labmat = bwlabel(bwim); % Create label matrix from bw image
    labdata = regionprops(labmat); % Get area and location of objects in label matrix

    %% Sort Objects by size
    ObjectsInfo = zeros(length(labdata),2); % initialize matrix
    NumComets = 0; % initialize number of comets
    % Loope through objects, keeping only objects that are larger than the
    % minimum allowed area
    for nobjects = 1:length(labdata)
        if labdata(nobjects).Area >= MINCOMETAREA && labdata(nobjects).Centroid(1) > 20 && labdata(nobjects).Centroid(2) > 20 && labdata(nobjects).Centroid(2) < size(bwim,1)-20 && labdata(nobjects).Centroid(1) < size(bwim,2)-20
            NumComets = NumComets +1;
            ObjectsInfo(NumComets,:) = labdata(nobjects).Centroid;
        else
            ObjectsInfo(end,:) = [];
        end
    end
    
    if NumComets > 1
    %% Sort objects by nearest neighbor distance
    % Calculate distances between objects and create square matrix of
    % distances
    distances = pdist(ObjectsInfo);
    ObjDist = squareform(distances);
    clear distances
    
    % filter objects based on distance relative to PIXEL_SPACING parameter
    uniquekeepers = [];
    for nmult = 1:NearestNeighbors
        [keepers{nmult}(:,1),keepers{nmult}(:,2)] = find(ObjDist < nmult*PIXEL_SPACING+25 & ObjDist > nmult*PIXEL_SPACING-25);
        uniquekeepers = [uniquekeepers;unique(keepers{nmult})];  % only keep unique objects
    end
    uniquekeepers = unique(uniquekeepers);
    NumComets = length(uniquekeepers);
    temp_objects = ObjectsInfo(uniquekeepers,:);
    clear ObjectsInfo
    ObjectsInfo = round(temp_objects); % round to integer values for indexing
    else
        NumComets = 0;
        ObjectsInfo = [];
    end

    
    %% Create a label image with black dots at each center
    %{
    BWlab2 = labmat;
    XCENT = 1;
    YCENT = 2;
    for obj = 1:NumComets
        BWlab2((ObjectsInfo(obj,YCENT)-20):(ObjectsInfo(obj,YCENT)+20),(ObjectsInfo(obj,XCENT)-1):(ObjectsInfo(obj,XCENT)+1)) = 100;
        BWlab2((ObjectsInfo(obj,YCENT)-1):(ObjectsInfo(obj,YCENT)+1),(ObjectsInfo(obj,XCENT)-20):(ObjectsInfo(obj,XCENT)+20)) = 100;
    end
    RGB = label2rgb(BWlab2, @jet, 'k');
    h1 = figure;
    imshow(RGB,'InitialMagnification',33);
    pause(.5)
    close(h1);
    %}
        
end

function NumAnalyzed = runcometanalysis(OrigCometImage,ObjectsInfo,NumComets,CROPWIDTH,...
                CROPHEIGHT,CROPOFFSET,convfac,hdthresh,tailthresh,headdia,fid0)
            
cometdata = -ones(NumComets,5);
XCENT = 1;
YCENT = 2;
NumAnalyzed = 0;
            for numcomet = 1:NumComets
                cropRect = [ObjectsInfo(numcomet,XCENT)-CROPWIDTH/2+CROPOFFSET,...
                    ObjectsInfo(numcomet,YCENT)-CROPHEIGHT/2,CROPWIDTH,CROPHEIGHT];
                CropCometImage = imcrop(OrigCometImage,cropRect);
                if size(CropCometImage,1)*size(CropCometImage,2) == (CROPWIDTH+1)*(CROPHEIGHT+1)
                    [cometdata(numcomet,1),cometdata(numcomet,2),cometdata(numcomet,3),cometdata(numcomet,4),...
                        cometdata(numcomet,5)]=analyzecomet(CropCometImage,convfac,hdthresh,tailthresh,headdia);
                end
                if cometdata(numcomet,4) > 0
                    fprintf(fid0,'%7.4f \t %7.4f \t %7.4f \t %7.4f \t %7.4f \n',cometdata(numcomet,:));
                    NumAnalyzed = NumAnalyzed+1;
                end
                clear CropCometImage
            end
end

function [headdna,taildna,otm,taillength,cometlength] = analyzecomet(cometim,convfactor,hdthresh,tailthresh,headdiameter)

% Funtion analyzecomet is used to calculate various comet parameters from
% an image of a comet.
%
% David K. Wood
% May 2008
%
% inputs:
% cometim - image of comet
% convfactor - pixel to distance conversion for objective/camera
%   combination used for imaging (pixels/micron)
% hdthresh - % of peak intensity where comet begins
% tailthresh - % of peak intensity where comet ends
% headdiameter - diameter of head in pixels
%
% outpus:
% headdna - % dna in comet head
% taildna - % dna in comet tail
% otm - olive tail moment in microns
% taillength - length of comet tail in microns
% cometlength - total length of comet in microns


%% Calculate Background
%cometim = cometim(:,:,2); % for RGB pictures
backgndwidth = round(0.1*length(cometim(:,1))); % make background 10% of total image width

background = sum(cometim(1:backgndwidth,:))/backgndwidth; % background is the total average intensity in each column from the top 10% of the image
cometsum = sum(cometim)/length(cometim(:,1))-background; % cometsum is the average intensity in each column minus the background

% Reverse polarity of comet image
% cometsum = cometsum([end:-1:1]);

%Default comet params to -1: this will indicate that these comets have not
%been analyzed and should not be written to data file
headdna = -1;
taildna = -1;
otm = -1;
taillength = -1;
cometlength = -1;
    
%% Define Comet Head and Tail
indices = find(cometsum > hdthresh*max(cometsum)); % column where comet begins 

if numel(indices) > 0 && indices(1) < 0.5*length(cometsum) % if beginning found, assign location to startx
    
    startx = indices(1); 
    indices = find(cometsum(startx:end) < tailthresh*max(cometsum)); % find column where comet ends
    
    if numel(indices) > 0 % if end found, assign location to endx
        
        endx = indices(1)+startx-1; 
        
        clear indices % clear unused variable
        
        %Find head/tail division
        temp_htx = startx + headdiameter; % approximate with head diameter
        
        % Search for division by looking for minimum of first derivative
        % around approximate head diameter
       % startsearch = round(temp_htx*0.9); %Original values
      %  endsearch = round(temp_htx*1.1); %Original values
        
         startsearch = round(temp_htx*0.7); % *Lizzie* testing head-tail division %use this for 50um wafer with head diameter "21um"
         endsearch = round(temp_htx*1.2);  % *Lizzie* testing head-tail division

        [trash,headtailx] = min(diff(cometsum(startsearch:endsearch)));
        headtailx = startsearch+headtailx; % assign value found to headtailx
        clear trash temp_htx startsearch endsearch
        
        if headtailx <= endx % if head is not longer than tail, calculate comet parameters 
            
            %Round indices to nearest integer
            headtailx=round(headtailx);
            endx=round(endx);
            startx=round(startx);
            
            tailcom = sum(cometsum(headtailx+1:endx).*(headtailx+1:1:endx))/sum(cometsum(headtailx+1:endx))*convfactor; % tail center-of-mass
            headcom = sum(cometsum(startx:headtailx).*(startx:1:headtailx))/sum(cometsum(startx:headtailx))*convfactor; % head center-of-mass
            headdna = sum(cometsum(startx:headtailx))/sum(cometsum(startx:endx))*100; % percentage dna in comet head
            taildna = sum(cometsum(headtailx:endx))/sum(cometsum(startx:endx))*100; % percentage dna in comet tail
            otm = (tailcom-headcom)*taildna/100; % olive tail moment in microns
            taillength = (endx-headtailx)*convfactor; % tail length in microns
            cometlength = (endx-startx)*convfactor; % total comet length in microns
  
%             % This code is for debugging
%             Always comment out this code after debugging
%             Quick way to comment out: 
%             1. Select all the lines in this section (lines 645-652)
%             % 2. Press "Command" + "/"
%             % To uncomment: Repeat step 1 and press "Command" + "t"  
%  
%                figure
%                subplot(2,1,1)
%                imshow(cometim)
%                line([startx startx],[1 length(cometim(:,1))],'Color','w')
%                line([headtailx headtailx],[1 length(cometim(:,1))],'Color','r')
%                line([endx endx],[1 length(cometim(:,1))],'Color','b')
%                subplot(2,1,2)
%                plot((cometsum)/max(cometsum)) 
% % % 
%                
        end
    end
end
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double
end

% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double
end

% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double
end


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double
end


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double
end


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

