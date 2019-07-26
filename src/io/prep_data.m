function caesar = prep_data(DATA_PATH)
    caesar = load_raw_data(DATA_PATH);

    fprintf('Reformatting table...'); tic
    % Beautify column names
    for i=1:length(caesar.Properties.VariableNames)
        caesar.Properties.VariableNames(i) = strrep(caesar.Properties.VariableNames(i), '_mm_', '');
        caesar.Properties.VariableNames(i) = strrep(caesar.Properties.VariableNames(i), '_kg_', '');    
        if caesar.Properties.VariableNames{i}(end) == '_'
            caesar.Properties.VariableNames{i}(end) = '';
        end
    end

    % Reorder columns
    caesar = movevars(caesar, 'Gender', 'After', 'SubjectNumber');
    caesar = movevars(caesar, 'Age_years', 'After', 'Gender');
    caesar.Properties.VariableNames{3} = 'Age';

    % Reformat columns
    caesar = reformat_family_income(caesar); % range string --> integer value
    caesar = reformat_birth_state(caesar);   % individual state name --> midwest, ...
    caesar = reformat_occupation(caesar);
    caesar = reformat_education(caesar);
    caesar = reformat_number_of_children(caesar);
    caesar = reformat_fitness(caesar);
    caesar = reformat_race(caesar);
    caesar = reformat_marital_status(caesar);
    caesar = reformat_car_model(caesar);
    caesar = reformat_shoe_size(caesar);
    caesar = reformat_jacket_size(caesar);
    caesar = reformat_pants_size(caesar);
    caesar = reformat_bra_size(caesar);
    
    % Reordering rows... This was necessary to set the reference group for
    % dummy variables
    caesar.BirthState = reordercats(categorical(caesar.BirthState),{'Midwest', 'Foreign',  'Northeast', 'South', 'West'}); % changing reference group for BirthState to Midwest
    caesar.CarModel = reordercats(categorical(caesar.CarModel),{'Non-sedan', 'Sedan',});
    
    % Additional variables
    BMI = caesar.Weight./(caesar.Stature*0.001).^2;
    ReportedBMI = caesar.ReportedWeight./(caesar.ReportedHeight*0.001).^2;
    AgeSquared = caesar.Age.^2;
    Exp = max(caesar.Age - caesar.Education - 6,0);
    ExpSquared = Exp.^2;
    WeightError = caesar.ReportedWeight - caesar.Weight;
    BMIError = ReportedBMI - BMI;
    HeightError = caesar.ReportedHeight - caesar.Stature;

    caesar = [caesar, array2table([BMI, ReportedBMI, WeightError, BMIError, HeightError, AgeSquared,  Exp, ExpSquared], 'VariableNames', {'BMI', 'ReportedBMI', 'WeightError', 'BMIError', 'HeightError', 'AgeSquared', 'Experience', 'ExperienceSquared'})];

    fprintf('DONE (%.4f seconds)\n', toc)
end
