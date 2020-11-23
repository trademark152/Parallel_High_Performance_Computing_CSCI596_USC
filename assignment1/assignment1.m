%% CSCI 596 - Assignement 1
    % Theoretical Peak Performance of a Computer
    % Author: Minh Tran, September 2020
close all; clc; clear all
%% Part 1: measuring computational complexity
    % import data
    addpath(genpath(pwd))
    table = load('mdtime.mat');
    num_atoms = table2array(table.MDtime(:,1));
    runtime = table2array(table.MDtime(:,2)); % in seconds
    
    % plot results
    figure;    
    %loglog(num_atoms, runtime);
    scatter(num_atoms, runtime)
    set(gca,'xscale','log')    
    set(gca,'yscale','log') 
    xlabel('Number of atoms')
    ylabel('Run tme (seconds)')
    
    % fit:
    hold on
    Const = polyfit(log10(num_atoms),log10(runtime), 1);
    alpha = Const(1);
    beta = Const(2);
    runtime_est = num_atoms.^alpha.*10^(beta);
    loglog(num_atoms, runtime_est);
%     
%     figure2 = Plot();
%     %figure2.Title = 'Average transveral concentration at BT'; % plot title
%     %figure2.BoxDim = [7, 5]; %[width, height]
%     figure2.FontName = 'Times'; %font style
%     figure2.YLabel = 'Runtime (seconds)'; % xlabel
%     figure2.XLabel = 'Number of Atoms'; %ylabel    
%     %figure2.XLim = [0, 40]; % set x axis limit
%     %figure2.XLim = [0, 0.5]; % set y axis limit
%     %figure2.YTick = [0, 0.5, 1]; %[tick1, tick2, .. ] % major tick
%     figure2.YMinorTick = 'on'; % 'on' or 'off' % minor tick
%     %figure2.XTick = [0, 0.5, 1]; %[tick1, tick2, .. ] % major tick
%     %figure2.XMinorTick = 'on'; % 'on' or 'off' % minor tick    
%     figure2.Colors = {[0,0,1],[0,1,0]};
%     %figure2.LineWidth = [1]; % line width   
%     %figure2.LineStyle =  {'-','--'}; % three line styles 
%     %figure2.Markers = {'o'};
%     %figure2.MarkerSpacing = [15];  
%     figure2.Legend = {'actual','fit'};
%     % Save? comment the following line if you do not want to save
%     set(gcf,'color','w');
%     export_fig 'part1.png' -q101 -native
%     
%     % display results:
%     disp('fitted value of alpha: ')
%     disp(alpha);
    
    
%% Input
num_processors = 1; % number of processors in a parallel node
num_core_per_processor = 8; % number of cores per multi-core processor
clock_speed = 2.3e9; % measured in Ghz or 1e9 cycles/second
num_fma_units = 1; % number of fused multiply-add circuit in each core
vector_size = 512; % measured in bits: number of double-precision operands held in each vector register

%% Output
% theoretical peak flop/s of a computer (unit: gigaflop/s or 1e9)
peak_flop = num_processors * num_core_per_processor * clock_speed *...
    (2*num_fma_units) * (vector_size/64)/1e9;

%% Plotting
