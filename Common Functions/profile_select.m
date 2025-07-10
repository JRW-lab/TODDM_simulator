function [profile_sel,num_frames] = profile_select(profile_names,frame_sel)

% Parameters
flag = true;

% Set up interface to MySQL server
clc
while flag
    clc;
    fprintf("___________________________________________________\n")
    fprintf("|__________________________________________________|\n")
    fprintf("||                                                ||\n")
    fprintf("|| Saved Profiles:                                ||\n")
    fprintf("||________________________________________________||\n")
    for i = 1:length(profile_names)
        fprintf("||                                                ||\n")
        row = sprintf("|| %d:",i);
        fprintf("%s",row)
        for j = 1:(50-length(char(row)))
            fprintf(" ")
        end
        fprintf("||\n")
        row = sprintf("|| %s",profile_names{i});
        fprintf("%s",row)
        for j = 1:(50-length(char(row)))
            fprintf(" ")
        end
        fprintf("||\n")

        1;

    end
    fprintf("||________________________________________________||\n")
    fprintf("||________________________________________________||\n")

    % Enter values
    fprintf("\n")
    profile_sel = str2double(input(' > Select profile: ', 's'));
    if frame_sel
        num_frames = str2double(input(' > Enter number of frames: ', 's'));
    end

    % Error display
    if profile_sel <= length(profile_names) && profile_sel > 0 && num_frames >= 0
        flag = false;
    else
        fprintf("\n   ❌ Enter valid values (>=0 frames)!\n");
        input('', 's')
    end

end
clc