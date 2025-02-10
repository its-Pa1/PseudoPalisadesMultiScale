clear all;
clc;
close all;
tic
main_2D;
M_total_with_g = M_total;
S_total_with_g = S_total;
main_2D_without_g;
%%
M_total  = M_total_with_g - M_total_0g;
S_total  = S_total_with_g - S_total_0g;
toc
%% videos
% These two loops save the solution video per each 3 days

% In this section Videos are saved in Videos_diff folder
% To check if Videos_diff folder exists otherwise to create
if not(isfolder('Videos_diff'))
    mkdir('Videos_diff')
end

% to check which video profile supports available in the machine
% if mp4 is not supported then avi format will be used
profiles = VideoWriter.getProfiles();
check_mp4_support = find(ismember({profiles.Name},'MPEG-4'));

if isempty(check_mp4_support)
    video_ext = '.avi';
    v_pro  = 'Motion JPEG AVI';
else
    video_ext = '.mp4';
    v_pro = 'MPEG-4';
end

videofile = VideoWriter(strcat('Videos_diff/Tumor', video_ext),v_pro);
videofile.Quality = 100;
videofile.FrameRate = 10;
open(videofile);
for i = 1: size(M_total,3)
    figure(3)
    %     set(gcf, 'Position',  [100, 600, 500, 500])
    surf(x,y,M_total(:,:,i)')
    view(0,90)
    colorbar
    shading interp
    colormap jet
    drawnow
    if (i==1)
        title('Initial glioma cells difference', 'Fontsize', 15);
    else
        title(['Glioma cells difference at t = ', num2str(3*(i-1)), ' days '], 'Fontsize', 15);
    end
    xlabel('X' , 'Fontsize', 15);
    ylabel('Y' , 'Fontsize', 15);
    F= getframe(gcf);
    writeVideo(videofile,F);
end
close(videofile);

videofile = VideoWriter(strcat('Videos_diff/Acidity', video_ext),v_pro);
videofile.Quality = 100;
videofile.FrameRate = 10;
open(videofile);
for i = 1: size(S_total,3)
    figure(4)
    %     set(gcf, 'Position',  [600, 600, 500, 500])
    surf(x,y,S_total(:,:,i)')
    view(0,90)
    colorbar
    shading interp
    colormap jet
    drawnow
    if (i==1)
        title('Initial acidity difference', 'Fontsize', 15);
    else
        title(['Acidity difference at t = ', num2str(3*(i-1)),' days '], 'Fontsize', 15);
    end
    xlabel('X' , 'Fontsize', 15);
    ylabel('Y' , 'Fontsize', 15);
    F= getframe(gcf);
    writeVideo(videofile,F)
    
end
close(videofile);

%% This loop saves the results for each 30 days in folder Plots_diff in png
% or eps format with days in file names

% to check wheather Plots folder exists otherwise it makes a folder Plots
if not(isfolder('Plots_diff'))
    mkdir('Plots_diff')
end

for i = 1:10:size(M_total,3)
    figure(5)
    %     set(gcf, 'Position',  [100, 600, 500, 500])
    surf(x,y,M_total(:,:,i)')
    view(0,90)
    colorbar
    shading interp
    colormap jet
    drawnow
    if (i==1)
        title('Initial glioma cells difference', 'Fontsize', 15);
    else
        title(['Glioma cells difference at t = ', num2str(3*(i-1)), ' days '], 'Fontsize', 15);
    end
    xlabel('X' , 'Fontsize', 15);
    ylabel('Y' , 'Fontsize', 15);
    %     saveas(gcf,sprintf('Plots_diff/TumorDiff%ddays',3*(i-1)),'epsc');%eps
    saveas(gcf,sprintf('Plots_diff/TumorDiff%ddays.png',3*(i-1)));%png
    
    figure(6)
    %     set(gcf, 'Position',  [600, 600, 500, 500])
    surf(x,y,S_total(:,:,i)')
    view(0,90)
    colorbar
    shading interp
    colormap jet
    drawnow
    if (i==1)
        title('Initial acidity difference', 'Fontsize', 15);
    else
        title(['Acidity difference at t = ', num2str(3*(i-1)),' days '], 'Fontsize', 15);
    end
    xlabel('X' , 'Fontsize', 15);
    ylabel('Y' , 'Fontsize', 15);
    %     saveas(gcf,sprintf('Plots_diff/AcidityDiff%ddays',3*(i-1)),'epsc');
    saveas(gcf,sprintf('Plots_diff/AcidityDiff%ddays.png',3*(i-1)))
end


%% uncomment to save the workspace
if not(isfolder('Mat_files'))
    mkdir('Mat_files')
end
save('Mat_files/main_diff_pH.mat');
%%
toc