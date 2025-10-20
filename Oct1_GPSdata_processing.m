clear; 
close all; 
clc;

csvFile = 'Taraz_Routes.csv';
gapMin = 15;
minPts = 3;
vThr = 1.0;
stepThr = 5.0;

% figure(1); axis off;
% text(0.5,0.5,'Analysis based on a SINGLE IoT Tracker for Medical Waste Container',...
%     'HorizontalAlignment','center','FontSize',14,'FontWeight','bold');
% title('Figure 1: IoT Tracker Device');

opts = detectImportOptions(csvFile);
T = readtable(csvFile, opts);

lowerNames = lower(T.Properties.VariableNames);
latVar = find(ismember(lowerNames, {'latitude','lat'}), 1);
lonVar = find(ismember(lowerNames, {'longitude','lon','lng'}), 1);
tsVar  = find(ismember(lowerNames, {'timestamp','createdat','updatedat','time'}), 1);
spdVar = find(ismember(lowerNames, {'speed','velocity'}), 1);
devVar = find(ismember(lowerNames, {'deviceid','device','device_id'}), 1);

assert(~isempty(latVar) && ~isempty(lonVar) && ~isempty(tsVar), ...
    'CSV must contain latitude, longitude, and timestamp columns.');

data = table;
data.latitude  = T{:, latVar};
data.longitude = T{:, lonVar};

rawTs = T{:, tsVar};
if isdatetime(rawTs)
    data.timestamp = rawTs;
else
    try
        data.timestamp = datetime(rawTs,'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''','TimeZone','UTC');
    catch
        data.timestamp = datetime(rawTs,'TimeZone','UTC');
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
    spd = double(T{:, spdVar});
    if median(spd(~isnan(spd))) < 20, spd = spd * 3.6; end
    data.speed_kmh = spd;
else
    data.speed_kmh = nan(height(data),1);
end

data = sortrows(data, {'deviceId','timestamp'});

data.dt_s = nan(height(data),1);
data.trip_id = zeros(height(data),1);
tripCounter = 0;
ids = unique(data.deviceId);
for d = 1:numel(ids)
    idx = find(data.deviceId == ids(d));
    ts  = data.timestamp(idx);
    dt  = [nan; seconds(diff(ts))];
    data.dt_s(idx) = dt;
    newTrip = isnan(dt) | (dt > gapMin*60);
    grp = cumsum(newTrip);
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

data.is_stationary_point = (data.speed_kmh < vThr) | (data.stepDist_m <= stepThr);
data.stationary_run_id = -1 * ones(height(data),1);

runId = 0;
G = findgroups(data.deviceId, data.trip_id);
keys = unique(G);
for g = 1:numel(keys)
    idx = find(G == keys(g));
    runStart = NaN; runCount = 0;
    for ii = 1:numel(idx)
        k = idx(ii);
        if data.is_stationary_point(k)
            if isnan(runStart), runStart = k; runCount = 1; else, runCount = runCount + 1; end
        else
            if ~isnan(runStart) && runCount >= minPts
                runId = runId + 1;
                data.stationary_run_id(runStart:runStart+runCount-1) = runId;
            end
            runStart = NaN; runCount = 0;
        end
    end
    if ~isnan(runStart) && runCount >= minPts
        runId = runId + 1;
        data.stationary_run_id(runStart:runStart+runCount-1) = runId;
    end
end

valid = data.stationary_run_id >= 0;
errAll = [];
cluster_stats = struct('cluster_id',{},'trip_id',{},'deviceId',{}, ...
    'point_count',{},'mean_lat',{},'mean_lng',{},'mean_error',{},'max_error',{},'std_dev',{});
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
    cluster_stats(c) = cs; %#ok<SAGROW>
    errAll = [errAll, errs]; %#ok<AGROW>
end

if isempty(errAll)
    error('No stationary clusters found. Adjust thresholds or minPts.');
end

mErr = mean(errAll);
xErr = max(errAll);
sErr = std(errAll);
pct = [50 68 95 99];
pv  = prctile(errAll, pct);

fprintf('\nGPS Error Statistics (all trips):\n');
fprintf('  Mean: %.2f m | Max: %.2f m | Std: %.2f m\n', mErr, xErr, sErr);
fprintf('  CEP50: %.2f m | P68: %.2f m | P95: %.2f m | P99: %.2f m\n', pv(1), pv(2), pv(3), pv(4));

w1 = mean(errAll <= (mErr + sErr)) * 100;
w2 = mean(errAll <= (mErr + 2*sErr)) * 100;
w3 = mean(errAll <= (mErr + 3*sErr)) * 100;
fprintf('\nGaussian Check (within mean + kσ):\n');
fprintf('  1σ: %.1f%% | 2σ: %.1f%% | 3σ: %.1f%%\n', w1, w2, w3);

nPoints = 100;
truePath = [linspace(0,1000,nPoints)', zeros(nPoints,1)];
sigmaUrban = 15;
noiseUrban = sigmaUrban * randn(nPoints,2);
measuredUrban = truePath + noiseUrban;
errorUrban = sqrt(sum((measuredUrban - truePath).^2,2));

figure(4);
subplot(1,2,1);
h1 = histogram(errorUrban, 15, 'FaceColor',[0.2 0.4 0.8], 'FaceAlpha',0.5, 'EdgeColor','none');
xlabel('Error (m)'); ylabel('Freq'); title('Theoretical'); 
grid on; 

hold on; 
ylim([0 16]);
mu1 = mean(errorUrban); [~,i1] = min(abs(h1.BinEdges(1:end-1) - mu1));
y1 = h1.Values(max(i1,1)); text(mu1, y1, sprintf('Mean: %.1f m', mu1),'Color','b','FontSize',12);

subplot(1,2,2);
h2 = histogram(errAll, 15, 'FaceColor',[0.8 0.2 0.2], 'FaceAlpha',0.5, 'EdgeColor','none');
xlabel('Error (m)'); ylabel('Freq'); title('Experimental'); grid on; hold on; ylim([0 16]);
mu2 = mean(errAll); [~,i2] = min(abs(h2.BinEdges(1:end-1) - mu2));
y2 = h2.Values(max(i2,1)); text(mu2, y2, sprintf('Mean: %.1f m', mu2),'Color','r','FontSize',12);

st = sort(errorUrban); ct = (1:numel(st))/numel(st)*100;
se = sort(errAll);    ce = (1:numel(se))/numel(se)*100;

figure(5);
plot(st, ct,'b-','LineWidth',2); hold on;
plot(se, ce,'r--','LineWidth',2);
xlabel('Error (m)'); ylabel('Cumulative (%)');
title('CDF: Theoretical vs Experimental'); legend('Theoretical','Experimental','Location','best'); grid on;
tp50 = prctile(errorUrban,50); ep50 = prctile(errAll,50);
plot(tp50,50,'ks','MarkerFaceColor','w'); text(tp50,50,sprintf(' %.1f m',tp50),'Color','b');
plot(ep50,50,'ks','MarkerFaceColor','w'); text(ep50,50,sprintf(' %.1f m',ep50),'Color','r');

custom_normpdf = @(x, mu, sigma) exp(-0.5*((x-mu)/sigma).^2) ./ (sigma*sqrt(2*pi));
nTheo = numel(errorUrban); hTheo = 1.06 * std(errorUrban) * nTheo^(-1/5);
xi_theo = linspace(min(errorUrban), max(errorUrban), 200);
f_theo  = arrayfun(@(x) mean(custom_normpdf(errorUrban, x, hTheo)), xi_theo);

nExp  = numel(errAll); hExp = 1.06 * std(errAll) * nExp^(-1/5);
xi_exp = linspace(min(errAll), max(errAll), 200);
f_exp  = arrayfun(@(x) mean(custom_normpdf(errAll, x, hExp)), xi_exp);

figure(6);
plot(xi_theo, f_theo,'b-','LineWidth',2); hold on;
plot(xi_exp,  f_exp, 'r--','LineWidth',2);
xlabel('Error (m)'); ylabel('Density');
title('KDE: Theoretical vs Experimental'); 


legend('Theoretical','Experimental'); grid on;
[pt, it] = max(f_theo); text(xi_theo(it), pt, sprintf(' Theo Peak: %.1f m', xi_theo(it)), 'Color','b');
[pe, ie] = max(f_exp ); text(xi_exp(ie),  pe, sprintf(' Exp Peak: %.1f m',  xi_exp(ie)),  'Color','r');

models = {'Our IoT module','Smartphone','Automotive GPS','High-precision GNSS'};
cep_values = [prctile(errAll,50), 6, 4, 2];
p95_values = [prctile(errAll,95), 12, 7.5, 3.5];

figure(7);
bar([cep_values; p95_values]','grouped');
xticklabels(models); ylabel('Error (m)');
title('Accuracy Comparison'); legend('CEP50','P95'); grid on;

fprintf('\nCluster Details:\n');
fprintf('%-8s %-8s %-12s %-12s %-12s %-12s %-12s\n', ...
    'Cluster','Trip','Mean Lat','Mean Lng','Points','MeanErr','MaxErr');
for c = 1:numel(cluster_stats)
    cs = cluster_stats(c);
    fprintf('%-8d %-8d %-12.6f %-12.6f %-12d %-12.2f %-12.2f\n', ...
        cs.cluster_id, cs.trip_id, cs.mean_lat, cs.mean_lng, ...
        cs.point_count, cs.mean_error, cs.max_error);
end

fprintf('\nTheoretical Models:\n');
fprintf('  Smartphone: CEP 2–10 m, P95 4–20 m\n');
fprintf('  Automotive: CEP 3–5 m, P95 5–10 m\n');
fprintf('  SBAS/WAAS: CEP 1–3 m, P95 2–5 m\n');

fprintf('\nESP Device (stationary clusters):\n');
fprintf('  CEP50: %.2f m | P95: %.2f m\n', prctile(errAll,50), prctile(errAll,95));
cep50_now = prctile(errAll,50);
if cep50_now <= 3
    cat = 'High-precision GNSS';
elseif cep50_now <= 5
    cat = 'Typical Automotive GPS';
else
    cat = 'Basic Smartphone GPS';
end
fprintf('  Category: %s\n', cat);

function d = haversine_m(lat1, lon1, lat2, lon2)
R = 6371000;
dLat = deg2rad(lat2 - lat1);
dLon = deg2rad(lon2 - lon1);
a = sin(dLat/2).^2 + cos(deg2rad(lat1)).*cos(deg2rad(lat2)).*sin(dLon/2).^2;
d = 2*R*atan2(sqrt(a), sqrt(1-a));
end
