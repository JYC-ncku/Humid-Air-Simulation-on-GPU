function [u] = generate_movie()

NO_RESULTS = 20;
M(NO_RESULTS) = struct('cdata', [], 'colormap', []);
f = figure;
f.Visible = 'off';
% Save the movie to a fcile
% X-dir_Velocity
% Temperature
% Mass_Fraction
% Relative_Humidity

v = VideoWriter('X-dir_Velocity_Constant_CV.mp4');
v.FrameRate = 5;    % Set frames per second (optional)
v.Quality = 100;     % Set quality from 0 to 100 (optional)
open(v);
for i = 1:1:NO_RESULTS
    filename = sprintf('Results_frame_%s.txt', generate_string(i-1));

A = load(filename);
    surf(reshape(A(:,1), 800, 1600), reshape(A(:,2), 800, 1600), reshape(A(:,3), 800, 1600)) 
    view(2)
    shading flat
    colorbar;
    xticks(0.1:0.1:0.9);
    yticks(0.05:0.05:0.45);
    title('X-dir Velocity between 0s ~ 1s (HLL & Constant CV)');
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
    end
end