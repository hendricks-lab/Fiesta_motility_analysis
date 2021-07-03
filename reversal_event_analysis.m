%Code to separate trajectories into runs based on reversal events
%Use this code for bi-directional motility

%Modified by D. Beaudet, Hendricks lab, 06/22/2021
clc, clear all, close all
set(0,'DefaultFigureWindowStyle','docked')
addpath('/Users/danielbeaudet/Desktop/Methods Paper/');

% Smoothing Filter: 
span=10;    
pwr=2;      

% Some variables:
gh=0;
ki=0;
kdp=0;
kb=0;

% Initialize variables: 
ab.position=[];
ab.jkr=[];
ab.rev_rate=[];
ab.run_bw_rev=[];
ab.time_bw_rev=[];
ab.Nrev=[];

roi.position=[];
roi.jkr=[];
roi.rev_rate=[];
roi.run_bw_reversal=[];
roi.time_bw_rev=[];
roi.Nrev=[];
roi.rev_rate=[];
roi.run_bw_reversal=[];
roi.Nrev=[];

flist{1}='/Users/danielbeaudet/Desktop/Methods Paper/';
DT=0.08;
Mot_file='*.mat';   
Mt_file='*.xlsx';  

Dirname='/Users/danielbeaudet/Desktop/Methods Paper/';

%% Folders
dmot=dir(Mot_file); % A new directory of only motility files with .mat config
dmt=dir(Mt_file);   % A new directory of only MT files with .xls config

%% Load all the folders
for kd=1:numel(flist)
cd(flist{1,kd});    % Load the folders
dmot=dir(Mot_file); % Pick the trajectory files in .mat format
dmt=dir(Mt_file);   % Pick the MT files in .xlsx format

%% Load all the files
for k=1:length(dmot) % k is a vector of length dmot
g=load(dmot(k).name); %load everything in dmot
display(dmot(k).name);  % display
mt_dat=xlsread(dmt(k).name);    % load everything in dmt
display(dmt(k).name);   % display

%% Load each molecule in the file
for km=1:numel(g.Molecule)
    
traj(km).position=[];
traj(km).jrev_all=[];
traj(km).jrev=[];
traj(km).rev_rate=[];
traj(km).run_bw_reversal=[];
traj(km).off_axis_displacement=[]; 
traj(km).time=[];
traj(km).time_bw_rev=[];
traj(km).Nrev=[];
ki=ki+1;

xk1=g.Molecule(km).Results(:,3);%medfilt1(g.Molecule(km).Results(:,3));
yk1=g.Molecule(km).Results(:,4);%medfilt1(g.Molecule(km).Results(:,4));
distk=g.Molecule(km).Results(:,5);    % Distance (calculated from x and y)
fwhmk=g.Molecule(km).Results(:,6);    % Full Width Half Maximum
intensityk=g.Molecule(km).Results(:,7);   % Intensity
errork=g.Molecule(km).Results(:,8);   % Position error 

if abs(max(xk1)-min(xk1))>200||abs(max(yk1)-min(yk1)) > 200
   
framek=g.Molecule(km).Results(:,1);   %Lists the frames in a trajectory
timek1=g.Molecule(km).Results(:,2);    %Lists the instantaneous times

timek=double(timek1);

xk=xk1;
yk=yk1;
timek=timek1;

mt_coords=[mt_dat(:,1),mt_dat(:,2)];

range_xk=range(xk);
range_yk=range(yk);

%% Just for the coordinates
if range_yk>range_xk
xk_old=xk;
yk_old=yk;
yk=xk_old;
xk=yk_old;
mt_coords_new=[mt_coords(:,2),mt_coords(:,1)];
mt_coords=mt_coords_new;
range_max=max(yk)-min(yk);
elseif range_xk>range_yk
mt_coords=[mt_coords(:,1),mt_coords(:,2)];
range_max=max(xk)-min(xk);
end

%% Interpolate 
jkmt=[find(diff(mt_coords(:,1))~=0);numel(mt_coords(:,1))]; 
mt_coords=mt_coords(jkmt,:);%[mt_dat(jkmt,6),mt_dat(jkmt,7)].*pixel_size;
xkm=median(xk); %Finds the median value of x-position of the molecule
ykm=median(yk); %Finds the median value of y-position of the molecule
xmt=min(mt_coords(:,1)):100:max(mt_coords(:,1));
ymt=interp1(mt_coords(:,1),mt_coords(:,2),xkm,'pchip'); %interp1(x of MT coords,y of MT coords,x-median postion of  LBC,piecewise cubic

%% Project displacement along polynomial fit to trajectory
fit_vector=[]; 

dat_vector=[0,0;diff(xk),diff(yk)]; % Create a vector of deltaX and deltaY with initial points being (0,0)
mt_vector=[0,0;diff(mt_coords)];    % Create a vector of delta MT coordinates
fit_vector(:,1)=interp1(mt_coords(:,1),mt_vector(:,1),xk,'pchip');  %Finds corresponding values of xk
fit_vector(:,2)=interp1(mt_coords(:,1),mt_vector(:,2),xk,'pchip');  
fit_norm=repmat(sqrt(sum(fit_vector.^2,2)),1,2); % Normalize the fit
p_traj=polyfit(mt_coords(:,1),mt_coords(:,2),2);
fit_vector_norm=fit_vector./fit_norm;   %fit vector normalization factor
delp=zeros(1,numel(xk));    % Create an array of zeros containing 1 column and the row length is equal to the length of xk
delp_off=zeros(1,numel(xk));    %Same as delp

for kv=1:numel(xk)
delp(kv)=dat_vector(kv,:)*fit_vector_norm(kv,:)';   %Difference in position between consecutive frames. Multiply the fit vector by a factor of fit normalization 
delp_off(kv)=dat_vector(kv,:)*(fit_vector_norm(kv,:)*[0,1;-1,0])';
end

clear kv;

%% Reversal Calculation
position=cumsum(delp);
position_off=cumsum(delp_off);

clear delp_off

position_index=1:numel(position);
xk_index=1:numel(xk);
yk_index=1:numel(yk);

%% Position smoothing
if numel(framek)<200
position=position';
position=padarray(position,200,'post'); 
position=position';
xk=xk;
xk=padarray(xk,200,'post');
yk=yk;
yk=padarray(yk,200,'post');
elseif numel(framek)>=50
position=position;
xk=xk;
yk=yk;
end

%% Savitsky Golay smoothing

position_smooth=smooth(position,span,'sgolay',pwr);  % 39 windows size is a good option, 6 is a good power option
position_smooth=position_smooth(position_index);
position=position(position_index);

xk_smooth=smooth(xk,span,'sgolay',pwr); 
xk_smooth=xk_smooth(xk_index);
xk=xk(xk_index);

yk_smooth=smooth(yk,span,'sgolay',pwr); 
yk_smooth=yk_smooth(yk_index);
yk=yk(yk_index);

delp_smooth = diff(position_smooth);
delp1=delp_smooth(1:(end-1));  %First value to the second last value
delp2=delp_smooth(2:end);  %Second value to the last value


%% Find reversal events
jrev_all=find(sign(delp1)~=sign(delp2)); %finds all reversals
jrev_all = jrev_all + 1;
pos_rev_all=position(jrev_all);
run_bw_reversal_all=diff(pos_rev_all);
j_rev_keep = find(abs(run_bw_reversal_all)>0);
jrev = jrev_all(j_rev_keep);  

jkr=[1,jrev',numel(xk)'];

Nrev=numel(jkr);

 pos_rev=position(jkr);
 time_rev=timek(jkr);
 run_bw_reversal=diff(pos_rev);

 timek_bw_reversal=diff(time_rev);   % time between reversal
 rev_rate=Nrev/(timek(end)-timek(1));

 dt_rev=DT*diff(jkr);
end

%% Figures

figure(1),
        xlabel('Frames');ylabel('Position (nm)');
        %xlim([0 4375]); ylim([-15000 6000]);
        hold on,
        plot(framek(jkr),position_smooth(jkr),'LineWidth',2.5);
                 

                     
                            traj(km).rev_rate=rev_rate;
                            gh=gh+1;
                            num_rev{gh}=Nrev;
                            traj(km).position=position;
                            traj(km).jrev_all=jrev_all;
                            traj(km).jrev=jkr;
                            traj(km).rev_rate=rev_rate;
                            traj(km).run_bw_reversal=run_bw_reversal;
                            traj(km).time=timek;
                            traj(km).time_bw_rev=dt_rev;
                            traj(km).Nrev=Nrev;
                            kdp=kdp+1;
                            pos_ab{kdp}=traj(km).position;
% Plot a single trajectory and show stationary, processive, and diffusive runs for on and off-axis plots
r= [xk_smooth, yk_smooth];
disp = sqrt(sum((r-r(1,:)).^2,2)); 
    
figure(km+6), hold on,   
subplot(2,1,2),                      
                    for jr=2:numel(jkr)
                          diff_pos_jk=disp(jkr(jr))-disp(jkr(jr-1));
                          if abs(diff_pos_jk)<=10
                              plot(framek(jkr(jr)),disp(jkr(jr)),'k.','MarkerSize',2);
                              hold on,
                              plot(framek(jkr(jr-1):jkr(jr)),disp(jkr(jr-1):jkr(jr))','Color','k','LineWidth',2);
                          elseif abs(diff_pos_jk)>10 && abs(diff_pos_jk)<200
                              plot(framek(jkr(jr)),disp(jkr(jr)),'r.','MarkerSize',2);
                              hold on,
                              plot(framek(jkr(jr-1):jkr(jr)),disp(jkr(jr-1):jkr(jr))','Color','r','LineWidth',2);
                          elseif abs(diff_pos_jk)>=200
                              plot(framek(jkr(jr)),disp(jkr(jr)),'b.','MarkerSize',2);
                              hold on,
                              plot(framek(jkr(jr-1):jkr(jr)),disp(jkr(jr-1):jkr(jr))','Color','b','LineWidth',2);
                          end
                    end
                                ylabel('Displacement')                          
end

roi.position=[roi.position,traj(1:km).position];
roi.jkr=[roi.jkr,traj(1:km).jrev];
roi.rev_rate=[roi.rev_rate,traj(1:km).rev_rate];
roi.run_bw_reversal=[roi.run_bw_reversal,traj(1:km).run_bw_reversal];
roi.time_bw_rev=[roi.time_bw_rev,traj(1:km).time_bw_rev];
roi.Nrev=[roi.Nrev,traj(1:km).Nrev];
end
end

ab.position=[ab.position,roi.position];
ab.jkr=[ab.jkr,roi.jkr];
ab.rev_rate=[ab.rev_rate,roi.rev_rate];
ab.run_bw_rev=[ab.run_bw_rev,roi.run_bw_reversal];
ab.time_bw_rev=[ab.time_bw_rev,roi.time_bw_rev];
ab.Nrev=[ab.Nrev,roi.Nrev];

%% Average run length
bw=0.1;

jkrun=find(isnan(ab.run_bw_rev)==0);
abs_run_bw_rev=abs(ab.run_bw_rev(jkrun));
jproc=find(abs_run_bw_rev>200); % finds runs that are processive > 400nm
jdiff=find(abs_run_bw_rev<200 & abs_run_bw_rev>10); %finds runs that are diffusive between
%10nm and 400nm
proc_run_bw_rev=ab.run_bw_rev(jkrun(jproc));
diff_run_bw_rev=ab.run_bw_rev(jkrun(jdiff));
jplus_proc=find(proc_run_bw_rev>0);
jminus_proc=find(proc_run_bw_rev<0);
jplus_diff=find(diff_run_bw_rev>0); %finds runs that are (+) diffusive
jminus_diff=find(diff_run_bw_rev<0); %finds runs that are (-) diffusive
int_prc=50;

xpproc=25:int_prc:1800;
npproc1_1=hist(proc_run_bw_rev(jplus_proc),xpproc);
npproc2_1=hist(-proc_run_bw_rev(jminus_proc),xpproc);

xpdiff=0:int_prc:1800;
npdiff1_1=hist(diff_run_bw_rev(jplus_diff),xpdiff);
npdiff2_1=hist(-diff_run_bw_rev(jminus_diff),xpdiff);

%Plus Processive
[plus,gof_plus,output_plus]=fit(xpproc',npproc1_1','exp1');
pr_plus=(1./abs(plus.b));

%Minus Processive
[minus,gof_minus,output_minus]=fit(xpproc',npproc2_1','exp1');
pr_minus=(1./abs(minus.b));

%Plus Diffusive
[plus,gof_plus,output_plus]=fit(xpdiff',npdiff1_1','exp1');
diff_plus=(1./abs(plus.b));

%Minus Diffusive
[minus,gof_minus,output_minus]=fit(xpdiff',npdiff2_1','exp1');
diff_minus=(1./abs(minus.b));

figure(6)
hb1=bar(0.1,pr_plus,'BarWidth',bw,'FaceColor','k','facealpha',1,'edgecolor','k');
hold on,
hb2=bar(0.1,-pr_minus,'BarWidth',bw,'FaceColor','k','facealpha',0.6,'edgecolor','k');
hold on,
hb3=bar(0.3,diff_plus,'BarWidth',bw,'FaceColor','r','facealpha',1,'edgecolor','k');
hold on,
hb4=bar(0.3,-diff_minus,'BarWidth',bw,'FaceColor','r','facealpha',0.6,'edgecolor','k');
hold on,
ylabel({'|L_{rev}| (nm)'});


%% Run length distribution Processive Runs
figure(3)
bar(xpproc',npproc1_1','FaceColor','k','facealpha',1);
hold on,
bar(xpproc',-npproc2_1','FaceColor','k','facealpha',0.6);
xlim([-50 3000]);
ylabel('frequency'),xlabel('displacement (nm)');
title('Processive runs');



%% Run length distribution Diffusive Runs
figure(4)
bar(xpdiff',npdiff1_1','FaceColor','r','facealpha',1);
hold on,
bar(xpdiff',-npdiff2_1','FaceColor','r','facealpha',0.6);
xlim([-50 3000]);
ylabel('frequency'),xlabel('displacement (nm)');
title('Diffusive runs');

%% Fraction of time of Processice vs. Diffusive runs

%Fraction of time of processive runs
proc_time_bw_rev=ab.time_bw_rev(jkrun(jproc)); 

T_plus1_proc=proc_time_bw_rev(jplus_proc);
T_minus1_proc=proc_time_bw_rev(jminus_proc);

Ttot1_proc=sum(T_plus1_proc) + sum(T_minus1_proc);

fraction_plus_final1_proc=sum(T_plus1_proc)./Ttot1_proc;
fraction_minus_final1_proc=sum(T_minus1_proc)./Ttot1_proc;

% Fraction of time of Diffusive runs
diff_time_bw_rev=ab.time_bw_rev(jkrun(jdiff)); 

T_plus1_diff=diff_time_bw_rev(jplus_diff);
T_minus1_diff=diff_time_bw_rev(jminus_diff);

Ttot1_diff=sum(T_plus1_diff) + sum(T_minus1_diff);

fraction_plus_final1_diff=sum(T_plus1_diff)./Ttot1_diff;
fraction_minus_final1_diff=sum(T_minus1_diff)./Ttot1_diff;

 %% Figure for fraction of time proc, diff, minus and plus ended runs
figure(5)

subplot(1,2,2)
hb1=bar(0.1,fraction_plus_final1_proc,'BarWidth',bw,'FaceColor','k','facealpha',1,'edgecolor','k');
hold on,
hb1=bar(0.1,-fraction_minus_final1_proc,'BarWidth',bw,'FaceColor','k','facealpha',0.6,'edgecolor','k');

hold on,

hb3=bar(0.3,fraction_plus_final1_diff,'BarWidth',bw,'FaceColor','r','facealpha',1,'edgecolor','k');
hold on,
hb4=bar(0.3,-fraction_minus_final1_diff,'BarWidth',bw,'FaceColor','r','facealpha',0.6,'edgecolor','k');

hold on,
legend([hb1 hb2 hb3 hb4],{'proc (+)','proc (-)','diff (+)','diff (-)'},'Location', 'northeast'); 
ylabel({'Fraction of Time of Proc/ Diff Runs','Minus-End            Plus-End'});
xlim([0 1]);


%% Average velocity displays the average of all proc or diff runs as length over time. 
Vplus_proc1=proc_run_bw_rev(jplus_proc)./proc_time_bw_rev(jplus_proc);
Vminus_proc1=proc_run_bw_rev(jminus_proc)./proc_time_bw_rev(jminus_proc);  

Vplus_diff1=diff_run_bw_rev(jplus_diff)./diff_time_bw_rev(jplus_diff);
Vminus_diff1=diff_run_bw_rev(jminus_diff)./diff_time_bw_rev(jminus_diff);  

figure,
subplot(1,2,1)
avg_velplus_proc=bar(0.1,mean(Vplus_proc1),'BarWidth',bw,'FaceColor','k','facealpha',1);
hold on,
avg_velminus_proc=bar(0.1,mean(Vminus_proc1),'BarWidth',bw,'FaceColor','k','facealpha',0.6);
xlim([0 0.5]);

hold on,
avg_velplus_diff=bar(0.3,mean(Vplus_diff1),'BarWidth',bw,'FaceColor','r','facealpha',1);
hold on,
avg_velminus_diff=bar(0.3,mean(Vminus_diff1),'BarWidth',bw,'FaceColor','r','facealpha',0.6);
xlim([0 0.5]);
ylabel('Average Velocity (nm/s)');
             
   % save('test.mat'); %saves all variables in the workspace.        
