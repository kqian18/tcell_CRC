function fc = log2FC_ClusterVsNonCluster_byRegion(agg, metricName)

    T = agg(:, {'slideName','MMR_status','RegionComp',metricName});

    % Wide: one row per slide
    wide = unstack(T, metricName, 'RegionComp', ...
        'GroupingVariables', {'slideName','MMR_status'});

    epsVal = 1e-6;

    fc = wide(:, {'slideName','MMR_status'});
    fc.log2FC_Tumor = log2((wide.Tumor_Cluster + epsVal) ./ ...
                           (wide.Tumor_NonCluster + epsVal));
    fc.log2FC_NonTumor = log2((wide.Non_Tumor_Cluster + epsVal) ./ ...
                              (wide.Non_Tumor_NonCluster + epsVal));
end

