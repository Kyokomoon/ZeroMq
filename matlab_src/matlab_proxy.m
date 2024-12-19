
clear java;
javaaddpath('/home/kkulakov/Downloads/jeromq-0.6.0/target/jeromq-0.6.0.jar')

import org.zeromq.ZMQ.*;
import org.zeromq.*;

port_api = 2111;
context = ZMQ.context(1);
socket_api_proxy = context.socket(ZMQ.REP);
socket_api_proxy.bind(sprintf('tcp://*:%d', port_api));

fprintf("Start")
figure(1);
global pauseFlag;
pauseFlag = false;
uicontrol('Style', 'pushbutton', 'String', 'Pause/Resume', ...
              'Position', [20, 20, 100, 30], ...
              'Callback', @(src, event) togglePause());

global iter_set_dist;

iter_set_dist = 0;

global index_dl;
index_dl = 1;
global dist_list;
%dist_list = [10, 50, 100, 500, 1000, 4000];
dist_list = [10, 100, 500, 1000, 2000, 3000, 5000, 100000, 150000, 200000, 250000, 300000];

global noise_list;
%noise_list = [10, -10, -30, -50, -70, -80, -90, -100, -110, -120, 40];
%noise_list = [-110, -100, -90, -80, -70, -50, -30, -20, -10, 0];
noise_list = [-100, -90, -80, -70];


global iter_set_noise;
iter_set_noise = 0;
global index_noise_list;
index_noise_list = 1;

while true
    
    if ~pauseFlag
        msg = socket_api_proxy.recv();
        if ~isempty(msg)
            % fprintf('received message [%d]\n', length(msg));
            if(length(msg) > 1000)
                % process_data(msg);
                msg = transmission_channel_model(msg);
                % msg = simulaion(msg);
                msg = convert_to_byte(msg);
                % fprintf("send to proxe: %d\n", length(msg));

            end
            socket_api_proxy.send(msg);
        end
    else
        pause(0.1);
    end
end



function togglePause()
    global pauseFlag;
    pauseFlag = ~pauseFlag;
end


function out_data = convert_to_byte(complex_array)
    single_array = single(complex_array);
    realPart = real(single_array);
    imaginaryPart = imag(single_array);
    floatArray = zeros(1, 2 * length(single_array));
    floatArray(1:2:end) = realPart;    
    floatArray(2:2:end) = imaginaryPart; 
    out_data = typecast(single(floatArray), 'uint8');
    
end

function out_data = CostHata(data, h_enb, h_ue, d)

    fc = 2560; % Частота в МГц
    hte = h_enb; % Высота передающей антенны в метрах
    hre = h_ue; % Высота приемной антенны в метрах
    Cm = 0; % Поправочный коэффициент для средних городов и пригородов
    
    % Расчет поправочного коэффициента для высоты приемной антенны
    a_hre = (1.1 * log10(fc) - 0.7) * hre - (1.56 * log10(fc) - 0.8);
    
    % Расчет потерь сигнала
    L = 46.3 + 33.9 * log10(fc) - 13.82 * log10(hte) - a_hre + (44.9 - 6.55 * log10(hte)) * log10(d) + Cm;
    
    % res = 10^(L/10);
    fprintf("d = %d, L = %f\n", d, L);
    out_data = data / L;
end


function out_data = simulaion(data_raw)
    global iter_set_dist;
    global index_dl;
    global dist_list;
    data_slice = data_raw;
    floatArray = typecast(data_slice, 'single');
    data = complex(floatArray(1:2:end), floatArray(2:2:end));
    % fprintf("size data = %d\n", length(data));
    CON = 2;
    pos_ENB = [100, 200, 50];
    pos_UE = [200, 200, 1.5];
    iter_set_dist = iter_set_dist + 1;
    distance = 1;
    if iter_set_dist > 1000
        distance = dist_list(index_dl);
        fprintf("new distance: %d, index: %d\n", distance, index_dl);
        index_dl = index_dl + 1;
        if index_dl > length(dist_list)
            index_dl = 1;
            
        end
        iter_set_dist = 0;
    end
    %distance = 400;
    if(CON == 1)
        distance = sqrt((pos_UE(1) - pos_ENB(1))^2 + (pos_UE(2) - pos_ENB(2))^2);
    end
    %fprintf("dist = %f\n", distance);
    mu = 0; % Среднее значение
    sigma = distance; % Стандартное отклонение
    n = length(data); % Количество точек
    % Генерация нормально распределённого шума
    noise = mu + sigma * randn(n, 1);
    % Ограничение шума до ±100
    noise = max(min(noise, 100), -100);
    noise = noise + noise * 1i;
    % data_cost = CostHata(data, pos_ENB(3), pos_UE(3), distance);
    data_cost = data / (distance / 10) + noise;
    out_data = data_cost;
    % out_data = data;
    %fprintf("data_cost = %f\n", data_cost);
end

function out_data = transmission_channel_model(data_raw)
    global iter_set_dist;
    global index_dl;
    global dist_list;
    global noise_list;
    global iter_set_noise;
    global index_noise_list;

    data_slice = data_raw;
    floatArray = typecast(data_slice, 'single');
    data = complex(floatArray(1:2:end), floatArray(2:2:end));
    %{
    if iter_set_dist > 1000
        distance = dist_list(index_dl);
        fprintf("new distance: %d, index: %d\n", distance, index_dl);
        index_dl = index_dl + 1;
        if index_dl > length(dist_list)
            index_dl = 1;
            
        end
        iter_set_dist = 0;
    end
    %}


    c = 3 * 1e8;
    Nb = 20;
    f0 = 2.56 * 1e9;
    B = 23 * 1e6;
    D1 = 10;
    Dn = 200;
    N0 = -100;
    Ts = 1/B;
    
    
    PRINT_DEBUG_INFO = 0;

    L = length(data); 
    % fprintf("L = %d\n", L);
    Smpy = zeros(1, length(data));
    
    %
    iter_set_dist = iter_set_dist + 1;
    distance = dist_list(index_dl);
    if iter_set_dist > 500
        distance = dist_list(index_dl);
        fprintf("new distance: %d, index: %d\n", distance, index_dl);
        index_dl = index_dl + 1;
        if index_dl > length(dist_list)
            index_dl = length(dist_list);
            
        end
        iter_set_dist = 0;
    end
    %
    D = randi([D1, Dn], 1, Nb);
    %{
    iter_set_noise = iter_set_noise + 1;
    N0 = noise_list(index_noise_list);
    if iter_set_noise > 2000
        N0 = noise_list(index_noise_list);
        index_noise_list = index_noise_list + 1;
        fprintf("new N0: %d, index: %d\n", N0, index_noise_list);
        
        if index_noise_list > length(noise_list)
            index_noise_list = length(noise_list);
        end
        iter_set_noise = 0;

    end
    %}
    distance = distance / 1000;
    %data = CostHata(data, 50, 1, distance);

    for i = 1:Nb
        if PRINT_DEBUG_INFO
            fprintf("i = %d\n", i);
        end
        tau = round((D(i) - D1) / (c * Ts));
        %fprintf("tau = %d\n", tau);
        G = c / (4 * pi * D(i) * f0);
        %k = L + round(tau);
        Si = data;
        if PRINT_DEBUG_INFO
            %fprintf("Di = %d, tau = %d, G = %f\n", D(i), tau, G);
        end
        for k = 1:(L+tau)
           if(k <= tau)
                Si(k) = 0;
           else
               Si(k) = data(k - tau);
           end
        end
        %}
        Si = Si .* G;
        if PRINT_DEBUG_INFO
            %fprintf("len(Si) = %d\n", length(Si));
        end
        Smpy = sum_array(Smpy, Si);
        % Smpy = Smpy + Si;
    end
    
    n = transpose(wgn(length(Smpy), 1, N0));
    Smpy = Smpy + n + (n * 1i);

    out_data = Smpy;
end

function res = sum_array(A, B)
    A = A(:).';
    B = B(:).';
    min_len = min(length(A), length(B));
    res = A(1:min_len) + B(1:min_len);
    if length(A) > min_len
        res = [res, A(min_len+1:end)];
    elseif length(B) > min_len
        res = [res, B(min_len+1:end)];
    end
end
