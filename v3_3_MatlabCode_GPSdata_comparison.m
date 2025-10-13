clear; close all; clc;

csvFile              = 'routes_v2.csv';
gapMinutesNewTrip    = 15;
minClusterPoints     = 3;
stationarySpeedKmh   = 1.0;
stationaryStepDist_m = 5.0;

figure(1); axis off;
text(0.5, 0.5, 'Analysis based on a SINGLE IoT Tracker for Medical Waste Container', ...
    'HorizontalAlignment','center','FontSize',14,'FontWeight','bold');
title('Figure 1: IoT Tracker Device');

opts = detectImportOptions(csvFile);
T    = readtable(csvFile, opts);

lowerNames = lower(T.Properties.VariableNames);
latVar = find(ismember(lowerNames, {'latitude','lat'}), 1);
lonVar = find(ismember(lowerNames, {'longitude','lon','lng'}), 1);
tsVar  = find(ismember(lowerNames, {'timestamp','createdat','updatedat','time'}), 1);
spdVar = find(ismember(lowerNames, {'speed','velocity'}), 1);
devVar = find(ismember(lowerNames, {'deviceid','device','device_id'}), 1);

assert(~isempty(latVar) && ~isempty(lonVar) && ~isempty(tsVar));

data = table;
data.latitude  = T{:, latVar};
data.longitude = T{:, lonVar};

rawTs = T{:, tsVar};
if isdatetime(rawTs)
    data.timestamp = rawTs;
else
    try
        data.timestamp = datetime(rawTs, 'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''','TimeZone','UTC');
    catch
        data.timestamp = datetime(rawTs, 'TimeZone','UTC');
    end
end

if ~isempty(devVar)
    devCol = T{:, devVar};
    if ~iscellstr(devCol) && ~isstring(devCol), devCol = string(devCol); end
    data.deviceId = string(devCol);
else
    data.deviceId = repmat("DEVICE-1", height(data), 1);
end

if ~isempty(spdVar)
    spd = T{:, spdVar};
    spd = double(spd);
    if median(spd(~isnan(spd))) < 20, spd = spd * 3.6; end
    data.speed_kmh = spd;
else
    data.speed_kmh = nan(height(data),1);
end

data = sortrows(data, {'deviceId','timestamp'});

data.dt_s     = nan(height(data),1);
data.trip_id  = zeros(height(data),1);

tripCounter = 0;
ids = unique(data.deviceId);
for d = 1:numel(ids)
    idx = find(data.deviceId == ids(d));
    ts  = data.timestamp(idx);
    dt  = [nan; seconds(diff(ts))];
    data.dt_s(idx) = dt;
    newTrip = isnan(dt) | (dt > gapMinutesNewTrip*60);
    grp     = cumsum(newTrip);
    data.trip_id(idx) = grp + tripCounter;
    tripCounter = tripCounter + max(grp);
end

data.stepDist_m = nan(height(data),1);
for k = 2:height(data)
    if data.trip_id(k) == data.trip_id(k-1)
        data.stepDist_m(k) = haversine_m(data.latitude(k-1), data.longitude(k-1), ...
                                         data.latitude(k),   data.longitude(k));
    end
end
if all(isnan(data.speed_kmh))
    data.speed_kmh = (data.stepDist_m ./ data.dt_s) * 3.6;
end

data.is_stationary_point = (data.speed_kmh < stationarySpeedKmh) | (data.stepDist_m <= stationaryStepDist_m);
data.stationary_run_id   = -1 * ones(height(data),1);

runId = 0;
G = findgroups(data.deviceId, data.trip_id);
keys = unique(G);
for g = 1:numel(keys)
    idx = find(G == keys(g));
    runStart = NaN; runCount = 0;
    for ii = 1:numel(idx)
        k = idx(ii);
        if data.is_stationary_point(k)
            if isnan(runStart), runStart = k; runCount = 1;
            else, runCount = runCount + 1; end
        else
            if ~isnan(runStart) && runCount >= minClusterPoints
                runId = runId + 1;
                data.stationary_run_id(runStart : runStart+runCount-1) = runId;
            end
            runStart = NaN; runCount = 0;
        end
    end
    if ~isnan(runStart) && runCount >= minClusterPoints
        runId = runId + 1;
        data.stationary_run_id(runStart : runStart+runCount-1) = runId;
    end
end

valid = data.stationary_run_id >= 0;
errAll = [];
cluster_stats = struct('cluster_id',{},'trip_id',{},'deviceId',{}, ...
                       'point_count',{},'mean_lat',{},'mean_lng',{}, ...
                       'mean_error',{},'max_error',{},'std_dev',{});
cidList = unique(data.stationary_run_id(valid));

for c = 1:numel(cidList)
    cid = cidList(c);
    idx = find(data.stationary_run_id == cid);
    lat = data.latitude(idx); lon = data.longitude(idx);
    lat0 = mean(lat); lon0 = mean(lon);
    errs = arrayfun(@(i) haversine_m(lat0, lon0, lat(i), lon(i)), 1:numel(idx));
    cs.cluster_id  = cid;
    cs.trip_id     = mode(data.trip_id(idx));
    cs.deviceId    = data.deviceId(idx(1));
    cs.point_count = numel(idx);
    cs.mean_lat    = lat0;
    cs.mean_lng    = lon0;
    cs.mean_error  = mean(errs);
    cs.max_error   = max(errs);
    cs.std_dev     = std(errs);
    cluster_stats(c) = cs;
    errAll = [errAll, errs];
end

if isempty(errAll)
    error('No stationary clusters found. Try lowering the threshold or minClusterPoints.');
end

mean_error      = mean(errAll);
max_error       = max(errAll);
std_dev_error   = std(errAll);
pctLevels       = [50 68 95 99];
pctVals         = prctile(errAll, pctLevels);

fprintf('\nGPS Error Statistics across ALL trips (multi-iteration):\n');
fprintf('  Mean: %.2f m | Max: %.2f m | Std: %.2f m\n', mean_error, max_error, std_dev_error);
fprintf('  CEP50: %.2f m | P68: %.2f m | P95: %.2f m | P99: %.2f m\n', pctVals(1), pctVals(2), pctVals(3), pctVals(4));

within_1sigma = mean(errAll <= (mean_error + std_dev_error)) * 100;
within_2sigma = mean(errAll <= (mean_error + 2*std_dev_error)) * 100;
within_3sigma = mean(errAll <= (mean_error + 3*std_dev_error)) * 100;
fprintf('\nGaussian Test (empirical within mean + kσ): 1σ %.1f%% | 2σ %.1f%% | 3σ %.1f%%\n', ...
    within_1sigma, within_2sigma, within_3sigma);

nPoints     = 100;
truePath    = [linspace(0,1000,nPoints)', zeros(nPoints,1)];
sigmaUrban  = 15;
noiseUrban  = sigmaUrban * randn(nPoints,2);
measuredUrban = truePath + noiseUrban;
errorUrban    = sqrt(sum((measuredUrban - truePath).^2,2));

figure(4);
subplot(1,2,1);
h_theo = histogram(errorUrban, 15, 'FaceColor',[0.2 0.4 0.8], 'FaceAlpha',0.5, 'EdgeColor','none');
xlabel('Error Distance (m)'); ylabel('Frequency'); title('Theoretical'); grid on; hold on; ylim([0 16]);
theo_mean = mean(errorUrban); [~,i_theo] = min(abs(h_theo.BinEdges(1:end-1) - theo_mean));
y_theo = h_theo.Values(max(i_theo,1)); text(theo_mean, y_theo, sprintf('Mean: %.1f m', theo_mean),'Color','b','FontSize',12);

subplot(1,2,2);
h_exp = histogram(errAll, 15, 'FaceColor',[0.8 0.2 0.2], 'FaceAlpha',0.5, 'EdgeColor','none');
xlabel('Error Distance (m)'); ylabel('Frequency'); title('Experimental'); grid on; hold on; ylim([0 16]);
exp_mean = mean(errAll); [~,i_exp] = min(abs(h_exp.BinEdges(1:end-1) - exp_mean));
y_exp = h_exp.Values(max(i_exp,1)); text(exp_mean, y_exp, sprintf('Mean: %.1f m', exp_mean),'Color','r','FontSize',12);

sorted_theo = sort(errorUrban);  cdf_theo = (1:numel(sorted_theo))/numel(sorted_theo)*100;
sorted_exp  = sort(errAll);      cdf_exp  = (1:numel(sorted_exp))/numel(sorted_exp)*100;

figure(5); plot(sorted_theo, cdf_theo,'b-','LineWidth',2); hold on;
plot(sorted_exp,  cdf_exp, 'r--','LineWidth',2);
xlabel('Error Distance (m)'); ylabel('Cumulative Probability (%)');
title('GPS Error CDF: Theoretical (Urban) vs. Experimental'); legend('Theoretical','Experimental','Location','best'); grid on;

theo_p50 = prctile(errorUrban,50); exp_p50 = prctile(errAll,50);
plot(theo_p50, 50, 'ks','MarkerFaceColor','w'); text(theo_p50, 50, sprintf(' %.1f m',theo_p50),'Color','b');
plot(exp_p50,  50, 'ks','MarkerFaceColor','w'); text(exp_p50,  50, sprintf(' %.1f m',exp_p50), 'Color','r');

custom_normpdf = @(x, mu, sigma) exp(-0.5*((x-mu)/sigma).^2) ./ (sigma*sqrt(2*pi));

nTheo = numel(errorUrban); hTheo = 1.06 * std(errorUrban) * nTheo^(-1/5);
xi_theo = linspace(min(errorUrban), max(errorUrban), 200);
f_theo  = arrayfun(@(x) mean(custom_normpdf(errorUrban, x, hTheo)), xi_theo);

nExp  = numel(errAll);      hExp  = 1.06 * std(errAll) * nExp^(-1/5);
xi_exp = linspace(min(errAll), max(errAll), 200);
f_exp  = arrayfun(@(x) mean(custom_normpdf(errAll, x, hExp)), xi_exp);

figure(6); plot(xi_theo, f_theo,'b-','LineWidth',2); hold on;
plot(xi_exp, f_exp, 'r--','LineWidth',2);
xlabel('Error Distance (m)'); ylabel('Density');
title('GPS Error Density: Theoretical (Urban) vs. Experimental'); legend('Theoretical','Experimental'); grid on;

[peak_theo, idxT] = max(f_theo); text(xi_theo(idxT), peak_theo, sprintf(' Theo Peak: %.1f m', xi_theo(idxT)), 'Color','b');
[peak_exp,  idxE] = max(f_exp ); text(xi_exp(idxE),  peak_exp,  sprintf(' Exp Peak: %.1f m', xi_exp(idxE)),  'Color','r');

models = {'Our IoT module','Smartphone','Automotive GPS','High-precision GNSS'};
cep_values     = [prctile(errAll,50), 6, 4, 2];
accuracy95_val = [prctile(errAll,95), 12, 7.5, 3.5];

figure(7);
bar([cep_values; accuracy95_val]', 'grouped');
xticklabels(models); ylabel('Error Distance (m)');
title('GPS Accuracy Comparison'); legend('CEP (50%)','95% Confidence'); grid on;

fprintf('\nCluster Details (multi-iteration):\n');
fprintf('%-8s %-8s %-12s %-12s %-12s %-12s %-12s\n', ...
    'Cluster','Trip','Mean Lat','Mean Lng','Points','MeanErr','MaxErr');
for c = 1:numel(cluster_stats)
    cs = cluster_stats(c);
    fprintf('%-8d %-8d %-12.6f %-12.6f %-12d %-12.2f %-12.2f\n', ...
        cs.cluster_id, cs.trip_id, cs.mean_lat, cs.mean_lng, ...
        cs.point_count, cs.mean_error, cs.max_error);
end

function d = haversine_m(lat1, lon1, lat2, lon2)
R = 6371000;
dLat = deg2rad(lat2 - lat1);
dLon = deg2rad(lon2 - lon1);
a = sin(dLat/2).^2 + cos(deg2rad(lat1)).*cos(deg2rad(lat2)).*sin(dLon/2).^2;
d = 2*R*atan2(sqrt(a), sqrt(1-a));
end
