%Code to separate phagosome trajectories into runs using change-point analysis
%Use this code when the observed motility is primarily uni-directional

%D. Beaudet, Hendricks lab, 06/22/2021
close all, clear all, warning off, set(0,'DefaultFigureWindowStyle','docked')
addpath('/Users/danielbeaudet/Desktop/Methods Paper/') %directory containing MSD analysis codes

bw=0.1;

delays = [1:1:10];

diff_runs_count=0;
proc_runs_count=0;
junk_runs_count=0;

%% Smoothing Filter Parameters: 
span=10;    %  
pwr=2;      %  

%% Data location

flist{1}='/Users/danielbeaudet/Desktop/Methods Paper/';
DT=0.08;
Mot_file='*.mat';   % A funciton ID to pick only .mat files
Mt_file='*.xlsx';

Dirname='/Users/danielbeaudet/Desktop/Methods Paper/';
count=0;

%% Analysis
for kd=1:numel(flist)

cd(flist{1,kd});    % Load the folders
dmot=dir(Mot_file);
dmt=dir(Mt_file);% Pick the trajectory files in .mat format

%% Data analysis

for k=1:length(dmot) % k is a vector of length dmot
    
    g=load(dmot(k).name); %load everything in dmot
    display(dmot(k).name);  % display
    mt_dat=xlsread(dmt(k).name);    % load everything in dmt
    display(dmt(k).name);   % display
    
        %% Molecule coordinates 
        for ktrack=1:numel(g.Molecule)

        f = g.Molecule(ktrack).Results(:,1); %frames
        x = g.Molecule(ktrack).Results(:,3); %x
        y = g.Molecule(ktrack).Results(:,4); %y
        t = g.Molecule(ktrack).Results(:,2); %t
        intens = g.Molecule(ktrack).Results(:,7); %intensity
        
        %% MT coordinates
        mt_coords=[mt_dat(:,1),mt_dat(:,2)];
        
        %% Change-point analysis input parameters
        exp_time = mean(diff(t)); %approx time between frames [s]
        Lwindow = 12; %length of sliding window
        dt = exp_time;
        msd_thresh = 1.0; %threshold for processive (>) or diffusive (<)
        msd_step = 0.3; %minimum difference in MSD between states
        l_min = 3; %minimum number of points in a trajectory
        
        %% Interpolate 
        jkmt=[find(diff(mt_coords(:,1))~=0);numel(mt_coords(:,1))]; 
        mt_coords=mt_coords(jkmt,:);%[mt_dat(jkmt,6),mt_dat(jkmt,7)].*pixel_size;
        xkm=median(x); %Finds the median value of x-position of the molecule
        ykm=median(y); %Finds the median value of y-position of the molecule
        xmt=min(mt_coords(:,1)):100:max(mt_coords(:,1));
        ymt=interp1(mt_coords(:,1),mt_coords(:,2),xkm,'pchip'); %interp1(x of MT coords,y of MT coords,x-median position of phagosome,piecewise cubic

        
        %% Project displacement along polynomial fit to trajectory
        fit_vector=[]; 
        dat_vector=[0,0;diff(x),diff(y)]; % Create a vector of deltaX and deltaY with initial points being (0,0)

        mt_vector=[0,0;diff(mt_coords)];    % Create a vector of delta MT coordinates
        fit_vector(:,1)=interp1(mt_coords(:,1),mt_vector(:,1),x,'pchip');  %Finds corresponding values of x
        fit_vector(:,2)=interp1(mt_coords(:,1),mt_vector(:,2),x,'pchip');  
        fit_norm=repmat(sqrt(sum(fit_vector.^2,2)),1,2); % Normalize the fit
        p_traj=polyfit(mt_coords(:,1),mt_coords(:,2),2); %polynomial equation for MT coordinates
        fit_vector_norm=fit_vector./fit_norm;   %fit vector normalization factor
        delp=zeros(1,numel(x));    % Create an array of zeros containing 1 column and the row length is equal to the length of x
        delp_off=zeros(1,numel(x));    %Same as delp
      
            for kv=1:numel(x)
                delp(kv)=dat_vector(kv,:)*fit_vector_norm(kv,:)';   %Difference in position between consecutive frames. Multiply the fit vector by a factor of fit normalization 
                delp_off(kv)=dat_vector(kv,:)*(fit_vector_norm(kv,:)*[0,1;-1,0])';
            end

        clear kv;
        
        %% Position index
        position=cumsum(delp);
        position_off=cumsum(delp_off);
        position_index=1:numel(position);
        x_index=1:numel(x);
        y_index=1:numel(y);
       
        %% Position smoothing
            if numel(f)<200
                position=position';
                position=padarray(position,200,'post'); 
                position=position';
                x=x;
                x=padarray(x,200,'post');
                y=y;
                y=padarray(y,200,'post');
            elseif numel(f)>=50
                position=position;
                x=x;
                y=y;
            end
            
        %% Savitsky Golay smoothing
        position_smooth=smooth(position,span,'sgolay',pwr);  % 39 windows size is a good option, 6 is a good power option
        position_smooth=position_smooth(position_index);
        position=position(position_index);

        x_smooth=smooth(x,span,'sgolay',pwr); 
        x_smooth=x_smooth(x_index);
        x=x(x_index);

        y_smooth=smooth(y,span,'sgolay',pwr); 
        y_smooth=y_smooth(y_index);
        y=y(y_index);
       

        %% Plot the position over time of each trajectory
        figure(1),
        xlabel('Frames');ylabel('Position (nm)');
        %xlim([0 4375.*DT]); ylim([-15000 6000]);
        hold on,
        plot(f,position_smooth,'LineWidth',2.5);
  
        
        %% Change-point Analysis
        
        tmsd_changepts = tMSD_2D(t,x_smooth,y_smooth,intens,Lwindow,dt,msd_thresh, msd_step,l_min);
        alphas= tmsd_changepts(:,2);
        run_type= tmsd_changepts(:,3);
        
        %% Separate runs based on change points
        r1{k}= [x_smooth, y_smooth, tmsd_changepts(:,2), tmsd_changepts(:,3)]; %(x,y,alphas,proc(1) or diff(0), intensity)
        chpts{k} = [1, (find(diff(r1{k}(:,4)))+1)', length(t)];% finds frame of each change point for each trajectory
         kchi=0;
            for kch=2:(length(chpts{k})) % finds all runs between change points within each trajectory
                 kchi=kchi+1;

                 jseg{k,kchi}= (find(t>=t(chpts{k}(kchi)) & t<t(chpts{k}(kch)))); %consec time points
                 r{k,kchi}= [x_smooth(jseg{k,kchi}), y_smooth(jseg{k,kchi})]; % x and y points of each run
                 
                 pos_of_run= position(jseg{k,kchi});
                 JSeG{k,kchi}= [x_smooth(jseg{k,kchi}), y_smooth(jseg{k,kchi}), alphas(jseg{k,kchi}), run_type(jseg{k,kchi})]; % x,y,alpha, proc(1) or diff(0) for each run
                 count=count+1;
                 
               
                        if (mean(JSeG{k,kchi}(:,3)) >=-2 & mean(JSeG{k,kchi}(:,4))<1) % finds diffusive runs
                            diff_runs_count= diff_runs_count+1;
                            pos_of_diff_run{diff_runs_count}= pos_of_run'; % finds position of diffusive runs
                            dir_of_run_d{diff_runs_count}=sum(diff(pos_of_diff_run{diff_runs_count})); % total displacement from start of each run, (-) and (+) directions
                            diff_run_time{diff_runs_count}= length(pos_of_diff_run{diff_runs_count}).*0.08;
                            r_diff{diff_runs_count}= [JSeG{k,kchi}(:,1),JSeG{k,kchi}(:,2)];

                 
                        elseif (mean(JSeG{k,kchi}(:,4))==1) % finds processive runs
                            proc_runs_count=proc_runs_count+1;
                            pos_of_proc_run{proc_runs_count}= (pos_of_run');
                            dir_of_run_p{proc_runs_count}= sum(diff(pos_of_proc_run{proc_runs_count})); % pos2-pos1+ pos3-pos2 + ... + pos(n+1)-pos(n)
                            proc_run_time{proc_runs_count}= length(pos_of_proc_run{proc_runs_count}).*0.08;
                            r_proc{proc_runs_count}= [JSeG{k,kchi}(:,1),JSeG{k,kchi}(:,2)];

                            
                        else
                            junk_runs_count=junk_runs_count+1;
                            Junk_JSeG{junk_runs_count}=JSeG{k,kchi};
                            disp_j{junk_runs_count} = hypot(diff(Junk_JSeG{junk_runs_count}(:,1)), diff(Junk_JSeG{junk_runs_count}(:,2)));
                            disp_tot_junk{junk_runs_count} = sum(disp_j{junk_runs_count});
                          
                        end

            end
            
        end
    
        end
        
end
%% Displacement of processive runs

jplus_proc= find(cellfun(@(xi) (xi>0), dir_of_run_p)); %finds plus-end proc runs (displacement)
jminus_proc= find(cellfun(@(xi) (xi<0), dir_of_run_p)); % finds minus-end proc runs (displacement)
JPLUS_proc= abs(cell2mat(dir_of_run_p(jplus_proc))); % plus-end proc runs run length (nm) (displacement)
JMINUS_proc= abs(cell2mat(dir_of_run_p(jminus_proc))); % minus-end proc runs run length (nm) (displacement)

TPLUS_proc= sum(cell2mat(proc_run_time(jplus_proc))); % total plus-end proc run run time (frames)
TMINUS_proc= sum(cell2mat(proc_run_time(jminus_proc))); % total minus-end proc run run time (frames)
TTOT_proc= TPLUS_proc+TMINUS_proc; % total plus- and minus- end  runs run time (frames)

int_proc=50;
xpproc=25:int_proc:3000; % for histograms always start with half the value of the intervals (i.e. 25:int_prc:3000 with int_prc=50)
npproc1_1=hist(JPLUS_proc,xpproc);
npproc2_1=hist(JMINUS_proc,xpproc);

%Plus Processive
[plusp,gof_plusp,output_plusp]=fit(xpproc',npproc1_1','exp1');
proc_plus= abs(TPLUS_proc./TTOT_proc);

%Minus Processive
[minusp,gof_minusp,output_minusp]=fit(xpproc',npproc2_1','exp1');
proc_minus= abs(TMINUS_proc./TTOT_proc);


figure(k+3) % histogram of plus- and minus- end processive runs; freq. vs. displacement (nm)
bar(xpproc',npproc1_1','FaceColor','k','facealpha',1);
hold on,
bar(xpproc',-npproc2_1','FaceColor','k','facealpha',0.6);
xlim([-50 3000]);
ylabel('frequency'),xlabel('displacement (nm)');
title('Processive runs');

%% Displacement of diffusive runs

%displacement
jplus_diff= find(cellfun(@(xi) (xi>0), dir_of_run_d)); 
jminus_diff=find(cellfun(@(xi) (xi<0), dir_of_run_d));
JPLUS_diff= cell2mat(dir_of_run_d(jplus_diff)); % plus-end diff runs run length (nm)
JMINUS_diff=abs(cell2mat(dir_of_run_d(jminus_diff))); % minus-end diff runs run length (nm)


TPLUS_diff= sum(cell2mat(diff_run_time(jplus_diff)));
TMINUS_diff= sum(cell2mat(diff_run_time(jminus_diff)));
TTOT_diff= TPLUS_diff+TMINUS_diff; % total number of frames 

xpdiff=35:int_proc:3000;
npdiff1_1=hist(JPLUS_diff,xpdiff);
npdiff2_1=hist(JMINUS_diff,xpdiff);

%Plus diffusive
[plusd,gof_plusd,output_plusd]=fit(xpdiff',npdiff1_1','exp1');
diff_plus= abs(TPLUS_diff./TTOT_diff);

%Minus diffusive
[minusd,gof_minusd,output_minusd]=fit(xpdiff',npdiff2_1','exp1');
diff_minus= abs(TMINUS_diff./TTOT_diff);

figure(k+4)
bar(xpdiff',npdiff1_1','FaceColor','k','facealpha',1);
hold on,
bar(xpdiff',-npdiff2_1','FaceColor','k','facealpha',0.6);
xlim([-50 3000]);
ylabel('frequency'),xlabel('displacement (nm)');
title('Diffusive runs');

                                        
 %% Figure for fraction of time proc, diff, minus and plus ended runs
figure(k+5)

subplot(1,2,2)
hb1=bar(0.1,proc_plus,'BarWidth',bw,'FaceColor','k','facealpha',1,'edgecolor','k');
hold on,
hb2=bar(0.1,-proc_minus,'BarWidth',bw,'FaceColor','k','facealpha',0.6,'edgecolor','k');

hold on,

hb3=bar(0.3,diff_plus,'BarWidth',bw,'FaceColor','r','facealpha',1,'edgecolor','k');
hold on,
hb4=bar(0.3,-diff_minus,'BarWidth',bw,'FaceColor','r','facealpha',0.6,'edgecolor','k');
 
hold on,

legend([hb1 hb2 hb3 hb4],{'proc (+)','proc (-)','diff (+)','diff (-)'},'Location', 'northeast'); 
ylabel({'Fraction of Time of Proc/ Diff Runs','Minus-End            Plus-End'});
xlim([0 1]);

%% Calculate average velocity for proc(+), proc(-), diff(+), and diff(-)
 
 jtime_proc_plus  = cell2mat(proc_run_time(jplus_proc));
 jtime_proc_minus = cell2mat(proc_run_time(jminus_proc));
 jtime_diff_plus  = cell2mat(diff_run_time(jplus_diff));
 jtime_diff_minus = cell2mat(diff_run_time(jminus_diff));
 
%plus-end proc vels 
 vel_proc_plus     = JPLUS_proc./jtime_proc_plus;
 ave_vel_proc_plus = mean(vel_proc_plus);
 
 %minus-end proc vels 
 vel_proc_minus     = JMINUS_proc./jtime_proc_minus;
 ave_vel_proc_minus = mean(vel_proc_minus);
 

 %plus-end diff vels 
 vel_diff_plus      = JPLUS_diff./jtime_diff_plus;
 ave_vel_diff_plus  = mean(vel_diff_plus);
 
%minus-end diff vels 
 vel_diff_minus     = JMINUS_diff./jtime_diff_minus;
 ave_vel_diff_minus = mean(vel_diff_minus);
 
figure(k+6)
hb1=bar(0.1,ave_vel_proc_plus,'BarWidth',bw,'FaceColor','k','facealpha',1,'edgecolor','k');
hold on,
hb2=bar(0.1,-ave_vel_proc_minus,'BarWidth',bw,'FaceColor','k','facealpha',0.6,'edgecolor','k');
hold on,
hb3=bar(0.3,ave_vel_diff_plus,'BarWidth',bw,'FaceColor','r','facealpha',1,'edgecolor','k');
hold on,
hb4=bar(0.3,-ave_vel_diff_minus,'BarWidth',bw,'FaceColor','r','facealpha',0.6,'edgecolor','k');
hold on,

legend([hb1 hb2 hb3 hb4],{'proc (+)','proc (-)','diff (+)','diff (-)', 'stat (+)', 'stat (-)'},'Location', 'northeast'); %hb5 hb6
ylabel({'Average Velocity','Minus-End            Plus-End'});
%xlim([0 1]);


  
  % save('test.mat'); %saves all variables in the workspace.           