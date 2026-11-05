function [u] = generate_movie()

    NO_RESULTS = 200;
    M(NO_RESULTS) = struct('cdata', [], 'colormap', []);
    
    f = figure;
    f.Visible = 'off';

    v = VideoWriter('Velocity_Streamlines_with_Wafers_10s.avi');
    v.FrameRate = 20;    
    v.Quality = 100;     
    open(v);

    % --- 網格參數設定 ---
    NX = 96;
    NY = 155;
    dx = 0.959 / NX;
    dy = 1.55 / NY;

    % =========================================================
    % 🦸‍♂️ 神救援區塊：直接在 MATLAB 內部把 B_mat (遮罩) 刻出來！
    % 這樣就不需要讀取 txt 的第 10 欄，保住了你跑了幾天的數據！
    % =========================================================
    B_mat = zeros(NY, NX);
    for j = 1:NY
        for k = 1:NX
            % 換算當前網格的中心點 X, Y 座標
            X = (k - 0.5) * dx;
            Y = (j - 0.5) * dy;
            
            % 1. 判斷是否為 FOUP 牆壁
            if (X > 0.0 && X < 0.344 && ((Y > 0.0 && Y < 0.662) || (Y > 0.982 && Y < 1.55)))
                B_mat(j, k) = 1.0;
            % 2. 判斷是否為 10 片晶圓
            elseif (X >= 0.039 && X <= 0.339)
                num_wafers = 10;
                start_Y = 0.68;
                pitch_Y = 0.3 / num_wafers;
                wafer_thickness = dy; % 強制為一個網格厚度
                
                for w = 0:(num_wafers-1)
                    y_bottom = start_Y + w * pitch_Y;
                    y_top = y_bottom + wafer_thickness;
                    if (Y >= y_bottom && Y <= y_top)
                        B_mat(j, k) = 1.0;
                        break;
                    end
                end
            end
        end
    end
    % =========================================================

    for i = 1:1:NO_RESULTS
        filename = sprintf('Results_frame_%s.txt', generate_string(i-1));
        
        if ~isfile(filename)
            fprintf('找不到檔案 %s，跳過...\n', filename);
            continue;
        end

        A = load(filename);
        
        X_mat = reshape(A(:,1), NY, NX);
        Y_mat = reshape(A(:,2), NY, NX);
        U_mat = reshape(A(:,4), NY, NX);  
        V_mat = reshape(A(:,5), NY, NX);  
        
        % 💡 這裡已經把讀取 A(:,10) 的那行刪掉了，直接用我們上面算好的 B_mat！
        
        cla; 
        streamslice(X_mat, Y_mat, U_mat, V_mat, 1.5);
        hold on;

        % 畫出灰色實體障礙物
        [row, col] = find(B_mat == 1.0);
        for k_idx = 1:length(row)
            rect_x = X_mat(row(k_idx), col(k_idx)) - dx/2;
            rect_y = Y_mat(row(k_idx), col(k_idx)) - dy/2;
            rectangle('Position', [rect_x, rect_y, dx, dy], ...
                      'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none');
        end

        hold off;

        axis equal; 
        axis([0 0.959 0 1.55]);
        xticks(0.1:0.1:0.959);
        yticks(0.1:0.1:1.55);
        xlabel('X (m)');
        ylabel('Y (m)');
        title(sprintf('FOUP Velocity Streamlines between 0s ~ 10s (Frame %d)', i-1));
        box on;
        
        M(i) = getframe(f);
        writeVideo(v, M(i));
        
        fprintf('已成功繪製並截圖: Frame %d\n', i-1);
    end

    f.Visible = 'on';
    close(v);
    disp('影片輸出大功告成！');

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