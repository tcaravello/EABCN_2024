%% Load the informationally-robust proxy

if strcmp(dataFrequency,'M')
    load('../../instrument/appendix/OilSurprisesMLogRefined')
    
    if instRefine
        proxyRaw = [oilProxiesWTIMrefined(:,ncontract)]; 
    else
        proxyRaw = [oilProxiesWTIMcensored(:,ncontract)]; 
    end

    smplStartProxyInd = find(strcmp(sampleDatesProxy,smplStartProxy));
    smplEndProxyInd   = find(strcmp(sampleDatesProxy,smplEndProxy));

    smplStartProxyVARInd = find(strcmp(sampleDates,smplStartProxy));
    smplEndProxyVARInd   = find(strcmp(sampleDates,smplEndProxy));
end

proxy = proxyRaw(smplStartProxyInd:smplEndProxyInd,:);   

[T,np] = size(proxy);
k = 1; % index of variable(s) to be instrumented