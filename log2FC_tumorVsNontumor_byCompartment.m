function fc = log2FC_tumorVsNontumor_byCompartment(agg, metricName)

    T = agg(:, {'slideName','MMR_status','RegionComp',metricName});

    wide = unstack(T, metricName, 'RegionComp', ...
        'GroupingVariables', {'slideName','MMR_status'});

    epsVal = 1e-6;

    fc = wide(:, {'slideName','MMR_status'});
    fc.log2FC_Cluster = log2((wide.Tumor_Cluster + epsVal) ./ ...
                             (wide.Non_Tumor_Cluster + epsVal));
    fc.log2FC_NonCluster = log2((wide.Tumor_NonCluster + epsVal) ./ ...
                                (wide.Non_Tumor_NonCluster + epsVal));
end

