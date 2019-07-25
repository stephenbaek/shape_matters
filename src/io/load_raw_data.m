function caesar = load_raw_data(data_path)
    % Read the excel sheets
    % TODO: Warning: Variable names were modified to make them valid MATLAB identifiers. The original names are saved in the VariableDescriptions property.
    fprintf('Reading demographics data...'); tic
    demographics = readtable(fullfile(data_path, 'Demographics.xls'));
    fprintf('DONE (%.4f seconds)', toc)

    fprintf('Reading measurement data...'); tic
    measurements = readtable(fullfile(data_path, 'Measurements.xls'));
    fprintf('DONE (%.4f seconds)', toc)

    % Combine them into one table
    caesar = join(demographics, measurements);

    % Remove subjects that do not have 3D scan data
    todel = [150; 264; 268; 307; 550; 1303; 2178; 2236]; % subject ids with missing 3D scans
    caesar(todel,:) = [];

    % read autoencoded (AE) body shape parameters
    load(fullfile(data_path, 'bodyshape.mat'));
    caesar = [caesar array2table(double(AE))]; % append the AE parameters to the table

    % TODO: This should be removed later (these columns are removed from the
    % excel sheet
    % clean and rearrange columns 
    caesar = removevars(caesar, {'Country', 'Date', 'Time', 'Civilian', 'SubgroupNumber', 'Recorder', 'Measurer'});
    caesar.Properties.VariableNames{'Age_years_'} = 'Age';
end