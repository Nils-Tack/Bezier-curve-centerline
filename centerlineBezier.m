%% Extracting centerlines using Bezier curve controlled by handles
% Author: Nils Tack
% Compiled with Matlab 2022a
% Script used to manually extract fish centerlines using Bezier curves
% whose shape is controlled interactively using the original image.

clear all
close all

% Options
opts.checkImage = 1;
opts.resetAxis = 1;
%% Set paths and directories
% Path to image stack
filepath = 'D:\Documents\MATLAB scripts\centerlineBezier\data';
D_image = dir(fullfile(filepath,'images','*.jpg')); % directory with images on interest

% Create directory for new kinematics files
D_kinematics = fullfile(filepath,'kinematics');
mkdir(D_kinematics);

% Set spatial and temporal scales
time = 1/2000; % set time interval between frames in seconds. Used in later sections. was 0.001
scale = 44.3491*1000; %0.000781*1000; % set the scale of the video from mm/px to  px/m

if opts.checkImage
    examplFr = 1;
    tempFig = figure('units','centimeters','Position',[1 1 20 20]); hold on
    A = flipud(imread(quickfilepath(D_image(examplFr))));
    imagesc(A,'XData',[0 size(A,2)/(scale/1000)],'YData',[0 size(A,1)/(scale/1000)]); colormap gray
    axis equal
    hold off
end

%% Extract profiles for all the frames
f1 = 1; % starting frame
incrBetweenFrames = 1; % increment between frames

% Loop through the images
close all % make sure that all the other figures are closed otherwise it may interfere with the scritp

fig = figure('units','centimeters','Position',[1 1 20 20]); % create new figure

% if ~isempty(tempFig)
% close(tempFig)
% end 

for i = f1:incrBetweenFrames:length(D_image)
    hold on
    A = flipud(imread(quickfilepath(D_image(i))));
    imagesc(A,'XData',[0 size(A,2)/(scale/1000)],'YData',[0 size(A,1)/(scale/1000)]); colormap gray
    axis equal

    % Plot the previous profile just to make a basic visual comparison when
    % tracing a new profile
    if i > f1
       plot(centerline(:,1),centerline(:,2),'-y','LineWidth',1')

    end
    
    inloop = 1; % Image loop
    inter = 0; % Initiate loop interruption variable. Set to 0 initially. Stopping the process will change it to 1.

    while inloop == 1
        title(sprintf('Frame %s',D_image(i).name))

        keyChar = waitForKeyPress();
        switch keyChar
            case 'a' % Update profile
                if opts.resetAxis
                % select area where fish needs to be measured - it resized
                % the image for convenience. It is not necessary but can
                % help. Change the size of the field of view 'fov' if
                % needed
                [xcView,ycView] = ginputc(1);
                fov = 15; % tight field of view around the appendage in mm
                axis([xcView-fov/2 xcView+fov/2 ycView-fov/2 ycView+fov/2]) % set the axes for the desirec field of view
                end

                % Define centerline
                n = 5; % number of points that drive the bezier curve (two end points and 2 handles). Add more points if needed. Keep in mind that for a bezier curve the only true points placed on the profiles are the beginning and end. All the other points are handles. 
                n1 = n-1;
                [xTempp,yTempp]=ginputc(n,'ShowPoints',true,'Color','r'); % select three points (2 handles and the distal end of the endo-exo)
                
                p = [xTempp,yTempp]; % coordinates of the points. First and last points are the beginning and end of the curve, other points are the handles

                for    ii=0:1:n1
                sigma(ii+1)=factorial(n1)/(factorial(ii)*factorial(n1-ii));  % for calculating (x!/(y!(x-y)!)) values 
                end
                l=[];
                UB=[];
                for u=0:0.002:1
                for d=1:n
                UB(d)=sigma(d)*((1-u)^(n-d))*(u^(d-1));
                end
                l=cat(1,l,UB); % catenation 
                end
                P=l*p;
           
                h_ax = gca; % define original set of axes serving as the base.

                % Move point handles
                % 1) Define variables
                x = p(:,1); % x coordinates of handles
                y = p(:,2); % y coordinates of handles
                
                X = P(:,1); % x coordinates of Bezier curve
                Y = P(:,2); % y coordinates of Bezier curve
                
                % 2)Set new axes for the curve handles
                ax = axes('Color', 'none','position', get(h_ax, 'position'),'Visible','off');
                axis equal
                ax.XLim = [h_ax.XLim];
                ax.YLim = [h_ax.YLim];
                
                % 3) Link the stacked axes to avoid a mess when panning or zooming
                linkaxes([fig.Children(1) fig.Children(2)],'xy') 
                
                % 4) Plot Bezier curve with handles
                hold on
                pl = plot(x,y,'--r',X,Y,'.r','hittest','on','buttondownfcn',@clickmarker); % vertex construction lines/points
            
            case 's' % Stop the loop if needed
                inter = 1; % Interrupt outer loop
                break;

            otherwise % exists the loop to save profile (hit return)
%                 set(fig, 'windowbuttonmotionfcn', []);  %disable motion callback
%                 goodBezier = [fig.Children(1).Children(1).XData',fig.Children(1).Children(1).YData'];
%                 plot(goodBezier(:,1),goodBezier(:,2),'-g')
                break; % Leave inner loop
        end 
    end

if inter == 1
   close all
   disp('Stopping')
   break
end

% Export profile
goodBezier = [fig.Children(1).Children(1).XData',fig.Children(1).Children(1).YData'];
centerline = goodBezier; % concatenate protopod and endo-exo

filename = sprintf('profile_%05g', i); % Use original number
            fprintf('Exporting: %s\n',filename)
            csvwrite(fullfile(D_kinematics,[filename,'.csv']),centerline)

% Clear all the elements in figure 'fig'
clf
hold off
end


%% Functions to drag points

function clickmarker(src,ev)
set(ancestor(src,'figure'),'windowbuttonmotionfcn',{@dragmarker,src})
set(ancestor(src,'figure'),'windowbuttonupfcn',@stopdragging)
end

function dragmarker(fig,ev,src)
%get current axes and coords
h1=fig.Children(1); % select which axes the data needs to be worked on
coords=get(h1,'currentpoint');

%get all x and y data 
x=h1.Children(2).XData;
y=h1.Children(2).YData;

%check which data point has the smallest distance to the dragged point
x_diff=abs(x-coords(1,1,1));
y_diff=abs(y-coords(1,2,1));
[value index]=min(x_diff+y_diff);

%create new x and y data and exchange coords for the dragged point
x_new=x;
y_new=y;
x_new(index)=coords(1,1,1);
y_new(index)=coords(1,2,1);

%update plot
% set(src,'xdata',x_new,'ydata',y_new);
set(h1.Children(2),'xdata',x_new,'ydata',y_new);


% Update bezier curve (at least try to)
h2=fig.Children(1); % select which axes the data needs to be worked on
xCurv = h2.Children(1).XData;
yCurv = h2.Children(1).YData;
    
% Re-calculate curve
    n = 5;
    n1 = n-1;
    for    i=0:1:n1
    sigma(i+1)=factorial(n1)/(factorial(i)*factorial(n1-i));  % for calculating (x!/(y!(x-y)!)) values 
    end
    l=[];
    UB=[];
    for u=0:0.002:1
        for d=1:n
        UB(d)=sigma(d)*((1-u)^(n-d))*(u^(d-1));
        end
        l=cat(1,l,UB); %catenation 
    end
    P=l*[x_new' y_new'];

%      set(src,'xdata',P(:,1),'ydata',P(:,2));
set(h1.Children(1),'xdata',P(:,1),'ydata',P(:,2));

end 

function stopdragging(fig,ev)
set(fig,'windowbuttonmotionfcn','')
set(fig,'windowbuttonupfcn','')
end

