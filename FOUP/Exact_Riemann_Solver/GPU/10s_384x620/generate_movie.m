function [u] = generate_movie()

NO_RESULTS = 10;
M(NO_RESULTS) = struct('cdata', [], 'colormap', []);
f = figure;
f.Visible = 'off';
% Save the movie to a fcile
% X-dir_Velocity
% Temperature
% Mass_Fraction
% Relative_Humidity

v = VideoWriter('Mass_Fraction_Exact_Riemann_Variable_CV_10s_GPU');
v.FrameRate = 1;    % Set frames per second (optional)
v.Quality = 100;     % Set quality from 0 to 100 (optional)
open(v);
for i = 1:1:NO_RESULTS
    filename = sprintf('Results_frame_%s.txt', generate_string(i-1));

A = load(filename);
    surf(reshape(A(:,1), 620, 384), reshape(A(:,2), 620, 384), reshape(A(:,8), 620, 384)) 
    view(2)
    shading flat
    colorbar;
    axis equal; 
    axis([0 0.959 0 1.55]);
    xticks(0.1:0.1:0.959);
    yticks(0.1:0.1:1.55);
    xlabel('X (m)');
    ylabel('Y (m)');
    % xticks(0.1:0.1:0.959);
    % yticks(0.1:0.1:1.55);
    title('Mass Fraction between 0s ~ 10s (Exact_Riemann & Variable CV & GPU)');
    M(i) = getframe(f);
    pause(0.01)
    writeVideo(v, M(i));

end
f.Visible = 'on'

% Good for testing without writing the mpeg file
% movie(M);

close(v);

end

function [num_string] = generate_string(number)
    if (number < 10)
        num_string = sprintf('00%d',number);
    elseif (number < 100)
        num_string = sprintf('0%d', number);
    else
        num_string = sprintf('%d', number);
    end
end