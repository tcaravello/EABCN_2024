function series_stat = stat_transform(series,method)

if method == 1
    series_stat = regcyc(series);
elseif method == 2
    [~,series_stat] = hpfilter(series,1600);
elseif method == 3
    [~,series_stat] = ls_detrend(series,2);
elseif method == 4
    [~,series_stat] = ls_detrend(series,3);
elseif method == 5
    [~,series_stat] = ls_detrend(series,4);
elseif method == 6
    [~,series_stat] = ls_detrend(series,5);
end