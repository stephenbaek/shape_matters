function caesar = load_raw_data(data_path)
    % Read the excel sheets
    warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
    fprintf('Reading demographics data...'); tic
    demographics = readtable(fullfile(data_path, 'Demographics.xls'));
    fprintf('DONE (%.4f seconds)\n', toc)

    fprintf('Reading measurement data...'); tic
    measurements = readtable(fullfile(data_path, 'Measurements.xls'));
    fprintf('DONE (%.4f seconds)\n', toc)

    % Combine them into one table
    caesar = join(demographics, measurements);

    % Remove subjects that do not have 3D scan data
    todel = [150; 264; 268; 307; 550; 1303; 2178; 2236]; % subject ids with missing 3D scans
    caesar(todel,:) = [];

    % read autoencoded (AE) body shape parameters
    load(fullfile(data_path, 'bodyshape.mat'));
    caesar = [caesar array2table(double(AE))]; % append the AE parameters to the table
end