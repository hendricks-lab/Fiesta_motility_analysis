% Malina K. Iwanski (Hendricks Lab, McGill University) 29 May 2019
% with Adam G. Hendricks (transient MSD over sliding window 10 April 2019)

% function to find local alpha-values using a rolling MSD over a window;
% assigns alpha value to middle point in window (except at end points the first and last calculated alpha-values are repeated)


function tmsd_changepts = tMSD_2D(t,x,y,intens,Lwindow,dt,msd_thresh, msd_step,l_min)
%% inputs: 
% t (timestamps of frames of trajectory [s])
% x (x positions of trajectory [um])
% y (y positions of trajectory [um])
% intens (intensity of each spot in trajectory - only so that it is interpolated to the same number of frames)
% Lwindow (number of frames over which MSD is locally calculated)
% dt (approx. time between frames [s])
% msd_thresh (alpha-value above which is processive, below which is diffusive)
% msd_step (minimum threshold for findchangepts function; minimum improvement in residual error; changepoints currently based on mean values, but can use rms, std, mean and slope)
% l_min (smallest no. of frames a section can be; minimum distance between two changepoints)

%% output:
% [frame number, local alpha-value, processive (1) or paused (0), interpolated velocities [um/s], interpolated intensities)

%% MSD analysis:
origt = t;
%interpolate data so that time between frames is consistent
xdata=interp1(t,x,[t(1):dt:t(end)],'linear')';
ydata=interp1(t,y,[t(1):dt:t(end)],'linear')';
intensdata=interp1(t,intens,[t(1):dt:t(end)],'nearest')';
r = [xdata, ydata];          
t = t(1):dt:t(end);

vels = sqrt(diff(xdata).^2+diff(ydata).^2)./diff(t); %calculate frame-to-frame velocity using interpolated data, strictly positive

delays = floor(Lwindow/4):ceil(Lwindow/2); %calculate which delays to use based on size of averaging window
logdelays = log10(delays.*dt);

Nk = numel(r(:,1)) - Lwindow; %number of windows
mmsd = zeros(numel(t),2); %to store local alpha-values

for k=1:Nk
    jk = k:(k+Lwindow);
    MSDk = MSD_2D({r(jk,:)},delays); %calculate MSD of trajectory within window
    pk = polyfit(logdelays,log10(MSDk),1); %calculate alpha-value using linear fit to log-log plot of MSD over time
    mmsd(k+ceil(Lwindow/2),:) = pk; %slope is alpha-value (first column)
end

for i = 1:ceil(Lwindow/2)
    mmsd(i,:) = mmsd(ceil(Lwindow/2)+1,:); %assign alpha-value to middle point in window
end
for i = size(mmsd,1)-Lwindow+ceil(Lwindow/2)+1:size(mmsd,1)
    mmsd(i,:) = mmsd(size(mmsd,1)-Lwindow+ceil(Lwindow/2),:); %for end points, repeat initial or final calculated alpha-value
end

%% plot displacement of trajectory and corresponding alpha-values before parsing
% figure
% disp = sqrt(sum((r-r(1,:)).^2,2)); 
% subplot(2,1,1), plot(t,disp), ylabel('Displacement')
% subplot(2,1,2), plot(t,mmsd(:,1)), xlabel('Time (s)'), ylabel('Alpha')

%% findchangepts 

proc = zeros(size(mmsd,1),1);
changepts = findchangepts((mmsd(:,1)),'MinThreshold',msd_step,'MinDistance',l_min);%,'statistic','std'); %changepoints identified are only used as putative places where transition paused <--> processive

num_changes = length(changepts);
change_ind = unique([1;changepts;size(mmsd,1)]); %add first and last frame
num_changes = length(change_ind);

%% parsing
for i =1:(num_changes)
    if i == num_changes %for last section
 
        meanmsd = mean((mmsd(change_ind(i):1:change_ind(end),1))); %mean alpha value between last changepoint and last frame
        
        msdquantile = quantile((mmsd(change_ind(i):1:change_ind(end))),[0.01 0.95]); %not currently used
        msdspan = msdquantile(2)-msdquantile(1); %not currently used

        if meanmsd >= msd_thresh % if mean alpha between sequential changepoints is above threshold --> processive
            proc(change_ind(i):1:change_ind(end),1) = 1;
        end
        
    else %for all other sections
        meanmsd = mean((mmsd(change_ind(i):1:change_ind(i+1),1))); %mean alpha value between sequential changepoints
        
        msdquantile = quantile((mmsd(change_ind(i):1:change_ind(i+1))),[0.01 0.95]); %not currently used
        msdspan = msdquantile(2)-msdquantile(1); %not currently used
       if meanmsd >= msd_thresh % if mean alpha between sequential changepoints is above threshold --> processive
            proc(change_ind(i):1:change_ind(i+1),1) = 1;
        end
    end
end
proc_frames = find(proc);
pause_frames = find(~proc);
frames = [proc_frames;pause_frames];
frames = sort(frames);

revert_times = interp1(t,t,origt,'nearest');
revert_times(isnan(revert_times)) = [];
revert_frames = round((revert_times-t(1))./dt) +1;
revert_frames(revert_frames < 1) = 1;
exclude_frames = setdiff(frames,revert_frames);
revert_proc = ismember(revert_frames,proc_frames);
revert_pause = ismember(revert_frames,pause_frames);
revert_proc_frames = find(revert_proc);
revert_pause_frames = find(revert_pause);
local_alphas = mmsd(revert_frames,1);
vels = [0;vels(:,1)];
revert_vels = vels(revert_frames,1);
revert_intens = intensdata(revert_frames,1);
revert_frames = sort([revert_proc_frames;revert_pause_frames]);


%%
tmsd_changepts = [revert_frames,local_alphas,revert_proc,revert_vels,revert_intens];

%plot parsed trajectory: displacement over time and local alpha-values over time; paused in red, processive in blue
figure 
disp = sqrt(sum((r-r(1,:)).^2,2)); 
subplot(2,1,1), hold on
% subplot(2,2,1), hold on
plot(proc_frames,disp(proc_frames),'b.') 
 ylabel('Displacement')
% hold off

% subplot(2,2,2), hold on
plot(pause_frames,disp(pause_frames),'r.') 
ylabel('Displacement')
hold off

subplot(2,1,2), hold on,
% subplot(2,2,3), hold on
plot(proc_frames,mmsd(proc_frames,1),'b.')
 xlabel('Frame'), ylabel('Alpha')
% hold off

% subplot(2,2,4), hold on
plot(pause_frames,mmsd(pause_frames,1),'r.')
xlabel('Frame'), ylabel('Alpha')
hold off


end