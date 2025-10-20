clear; 
close all; 
clc;

data = readtable('Taraz_Routes.csv');

tracking_data = struct();
for i = 1:height(data)
    tracking_data(i).latitude  = data.latitude(i);
    tracking_data(i).longitude = data.longitude(i);
    tracking_data(i).speed     = data.speed(i);
    tracking_data(i).timestamp = data.timestamp(i);
end

speed_threshold = 1.0;
is_stationary   = [tracking_data.speed] < speed_threshold;

cluster_id = 0;
current_cluster = [];
all_clusters = {};

for i = 1:length(is_stationary)
    if is_stationary(i)
        current_cluster = [current_cluster, i];
    else
        if numel(current_cluster) >= 3
            cluster_id = cluster_id + 1;
            all_clusters{cluster_id} = current_cluster; %#ok<SAGROW>
        end
        current_cluster = [];
    end
end
if numel(current_cluster) >= 3
    cluster_id = cluster_id + 1;
    all_clusters{cluster_id} = current_cluster;
end

fprintf('Found %d stationary clusters with ≥ 3 points\n', numel(all_clusters));

num_clusters = numel(all_clusters);
cluster_stats = struct('cluster_id',[],'point_count',[],'mean_lat',[],'mean_lng',[], ...
                       'mean_error',[],'max_error',[],'std_dev',[]);
cluster_stats(num_clusters).cluster_id = [];
all_errors = [];

for c = 1:num_clusters
    idxs = all_clusters{c};
    if any(idxs > numel(tracking_data)) || any(idxs <= 0)
        warning('Cluster %d contains invalid indices - skipping', c);
        continue;
    end

    lats = zeros(1, numel(idxs));
    lngs = zeros(1, numel(idxs));
    for j = 1:numel(idxs)
        lats(j) = tracking_data(idxs(j)).latitude;
        lngs(j) = tracking_data(idxs(j)).longitude;
    end

    mean_lat = mean(lats);
    mean_lng = mean(lngs);

    R = 6371000;
    errs = zeros(1, numel(idxs));
    for k = 1:numel(idxs)
        lat = tracking_data(idxs(k)).latitude;
        lng = tracking_data(idxs(k)).longitude;
        dLat = deg2rad(lat - mean_lat);
        dLng = deg2rad(lng - mean_lng);
        a = sin(dLat/2).^2 + cos(deg2rad(mean_lat)) .* cos(deg2rad(lat)) .* sin(dLng/2).^2;
        c_val = 2 * atan2(sqrt(a), sqrt(1 - a));
        errs(k) = R * c_val;
    end

    cluster_stats(c).cluster_id  = c;
    cluster_stats(c).point_count = numel(idxs);
    cluster_stats(c).mean_lat    = mean_lat;
    cluster_stats(c).mean_lng    = mean_lng;
    cluster_stats(c).mean_error  = mean(errs);
    cluster_stats(c).max_error   = max(errs);
    cluster_stats(c).std_dev     = std(errs);

    all_errors = [all_errors, errs]; %#ok<AGROW>
end

if isempty(all_errors)
    error('No stationary-cluster errors were computed. Check speed threshold and data.');
end

mErr = mean(all_errors);
xErr = max(all_errors);
sErr = std(all_errors);

fprintf('\nGPS Error Statistics from stationary clusters:\n');
fprintf('  Mean Error: %.2f meters\n', mErr);
fprintf('  Maximum Error: %.2f meters\n', xErr);
fprintf('  Standard Deviation: %.2f meters\n', sErr);

sorted_errors   = sort(all_errors);
percentiles     = [50, 68, 95, 99];
percentile_vals = prctile(sorted_errors, percentiles);

fprintf('\nError Percentiles:\n');
for i = 1:numel(percentiles)
    fprintf('  %dth Percentile: %.2f meters\n', percentiles(i), percentile_vals(i));
end

w1 = mean(all_errors <= (mErr + sErr)) * 100;
w2 = mean(all_errors <= (mErr + 2*sErr)) * 100;
w3 = mean(all_errors <= (mErr + 3*sErr)) * 100;

fprintf('\nGaussian Distribution Test:\n');
fprintf('  Mean: %.2f meters\n', mErr);
fprintf('  Standard Deviation: %.2f meters\n', sErr);
fprintf('  Within mean + 1σ (%.2fm): %.2f%% (expected: ~68%%)\n', mErr + sErr, w1);
fprintf('  Within mean + 2σ (%.2fm): %.2f%% (expected: ~95%%)\n', mErr + 2*sErr, w2);
fprintf('  Within mean + 3σ (%.2fm): %.2f%% (expected: ~99.7%%)\n', mErr + 3*sErr, w3);

% figure('Name','GPS Error Histogram');
% histogram(all_errors, 10, 'Normalization', 'count');
% grid on; xlabel('Error Distance (meters)'); ylabel('Frequency');
% title('GPS Error Distribution');

figure('Name','GPS Error CDF');
x  = sorted_errors(:).';
n  = numel(x);
f  = (1:n) / n;
x0 = [0, x];
f0 = [0, f];
stairs(x0, f0*100, 'LineWidth', 2); hold on;
for i = 1:numel(percentiles)
    xl = percentile_vals(i);
    plot([xl xl], [0 100], 'r--');
    ylbl = min(98, percentiles(i) + 3);
    text(xl, ylbl, sprintf('%d%% (%.1fm)', percentiles(i), xl), ...
        'HorizontalAlignment','right','VerticalAlignment','bottom');
end
xlabel('Error Distance (meters)');
ylabel('Cumulative Probability (%)');
title('GPS Error Cumulative Distribution');
grid on; 
ylim([0 100]);
xlim([0, max(25, max(x)*1.10)]);

% 
% figure('Name','GPS Position Scatter (Stationary Clusters)'); hold on; grid on;
% colors = hsv(max(1, numel(all_clusters)));
% theta = linspace(0, 2*pi, 200);
% for r = [5 10 20 30]
%     plot(r*cos(theta), r*sin(theta), 'k--', 'LineWidth', 0.5);
%     text(r, 0, sprintf('  %dm', r), 'FontSize', 8);
% end
% R = 6371000;
% for c = 1:num_clusters
%     idxs = all_clusters{c};
%     if isempty(idxs), continue; end
%     xpts = zeros(1, numel(idxs));
%     ypts = zeros(1, numel(idxs));
%     lat0 = cluster_stats(c).mean_lat;
%     lng0 = cluster_stats(c).mean_lng;
%     for k = 1:numel(idxs)
%         lat = tracking_data(idxs(k)).latitude;
%         lng = tracking_data(idxs(k)).longitude;
%         dlat = deg2rad(lat - lat0);
%         dlng = deg2rad(lng - lng0);
%         xpts(k) = R * cos(deg2rad(lat0)) * dlng;
%         ypts(k) = R * dlat;
%     end
%     scatter(xpts, ypts, 36, colors(c,:), 'filled');
% end
% axis equal;
% xlim([-40 40]); ylim([-40 40]);
% xlabel('East–West Offset (meters)');
% ylabel('North–South Offset (meters)');
% title('GPS Position Scatter Plot (Stationary Clusters)');
% 
% models = {'Your Device','Basic Smartphone','Automotive GPS','High-precision GNSS'};
% cep_values        = [percentile_vals(1), 6, 4, 2];
% accuracy95_values = [percentile_vals(3), 12, 7.5, 3.5];
% 
% figure('Name','GPS Accuracy Comparison');
% bar([cep_values; accuracy95_values]'); grid on;
% xticklabels(models); ylabel('Error Distance (meters)');
% title('GPS Accuracy Comparison'); legend('CEP (50%)','95% Confidence');
% 
% fprintf('\nCluster Details:\n');
% fprintf('%-10s %-10s %-15s %-15s %-15s %-15s %-15s\n', ...
%     'Cluster', 'Points', 'Mean Lat', 'Mean Lng', 'Mean Error', 'Max Error', 'Std Dev');
% for c = 1:numel(cluster_stats)
%     if isempty(cluster_stats(c).cluster_id), continue; end
%     fprintf('%-10d %-10d %-15.5f %-15.5f %-15.2f %-15.2f %-15.2f\n', ...
%         cluster_stats(c).cluster_id, cluster_stats(c).point_count, ...
%         cluster_stats(c).mean_lat, cluster_stats(c).mean_lng, ...
%         cluster_stats(c).mean_error, cluster_stats(c).max_error, cluster_stats(c).std_dev);
% end
% 
% fprintf('\nTheoretical GPS Accuracy Models (reference):\n');
% fprintf('1) Basic Smartphone GPS (no SBAS): CEP 2–10 m, 95%% 4–20 m\n');
% fprintf('2) Typical Automotive GPS: CEP 3–5 m, 95%% 5–10 m\n');
% fprintf('3) High-precision GNSS (SBAS/WAAS): CEP 1–3 m, 95%% 2–5 m\n');
% 
% fprintf('\nESP GPS Device Performance:\n');
% fprintf('   - CEP (50%%): %.2f meters\n', percentile_vals(1));
% fprintf('   - 95%% accuracy: %.2f meters\n', percentile_vals(3));
% if percentile_vals(1) <= 3
%     category = 'High-precision GNSS';
% elseif percentile_vals(1) <= 5
%     category = 'Typical Automotive GPS';
% else
%     category = 'Basic Smartphone GPS';
% end
% fprintf('   - Category (approx): %s\n', category);
% 
% figure('Name','ECDF (exact)'); 
% x = sort(all_errors(:));
% f = (1:numel(x)).'/numel(x);
% plot([0; x], [0; f]*100, 'LineWidth', 2); grid on;
% xlim([0, max(25, max(x)*1.1)]); ylim([0 100]);
% xlabel('Error Distance (meters)'); ylabel('Cumulative Probability (%)');
% title('GPS Error Cumulative Distribution (connected ECDF)');
% 
% figure('Name','ECDF (kernel-smoothed)');
% x = sort(all_errors(:));
% xx = linspace(0, max(25, max(x)*1.1), 400);
% [Fc, xx] = ksdensity(x, xx, 'Function', 'cdf');
% plot(xx, Fc*100, 'LineWidth', 2); grid on;
% xlim([0, max(xx)]); ylim([0 100]);
% xlabel('Error Distance (meters)'); ylabel('Cumulative Probability (%)');
% title('GPS Error Cumulative Distribution (kernel-smoothed)');
% 
% figure('Name','ECDF (lognormal fit)');
% x = sort(all_errors(:));
% pd = fitdist(x,'lognormal');
% xx = linspace(0, max(25, max(x)*1.1), 400);
% Fc = cdf(pd, xx);
% plot(xx, Fc*100, 'LineWidth', 2); grid on;
% xlim([0, max(xx)]); ylim([0 100]);
% xlabel('Error Distance (meters)'); ylabel('Cumulative Probability (%)');
% title('GPS Error Cumulative Distribution (lognormal fit)');
