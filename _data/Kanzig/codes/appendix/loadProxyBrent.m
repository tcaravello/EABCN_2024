%% Load the proxy based on Brent futures

if strcmp(dataFrequency,'M')
    load('../../instrument/appendix/OilSurprisesMLogBrent')

    proxyRaw = [oilProxiesBrentM(:,ncontract)]; 

    smplStartProxyInd = find(strcmp(sampleDatesProxy,smplStartProxy));
    smplEndProxyInd   = find(strcmp(sampleDatesProxy,smplEndProxy));

    smplStartProxyVARInd = find(strcmp(sampleDates,smplStartProxy));
    smplEndProxyVARInd   = find(strcmp(sampleDates,smplEndProxy));

end

proxy = proxyRaw(smplStartProxyInd:smplEndProxyInd,:);   % we loose the first p values in the VAR

[T,np] = size(proxy);
k = 1; % index of variable(s) to be instrumented