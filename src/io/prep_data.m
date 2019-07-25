clear all
close all
clc

DATA_PATH = '../../Data';

caesar = load_raw_data(DATA_PATH);

%% reformat family income
for i=1:length(caesar.Properties.VariableNames)
    caesar.Properties.VariableNames(i) = strrep(caesar.Properties.VariableNames(i), '_mm_', '');
    caesar.Properties.VariableNames(i) = strrep(caesar.Properties.VariableNames(i), '_kg_', '');    
    if caesar.Properties.VariableNames{i}(end) == '_'
        caesar.Properties.VariableNames{i}(end) = '';
    end
end
caesar = movevars(caesar, 'Gender', 'After', 'SubjectNumber');
caesar = movevars(caesar, 'Age', 'After', 'Gender');
caesar.FamilyIncome(strcmp(caesar.FamilyIncome, 'Less than 10k$/10MLira')) = {'7500'};
caesar.FamilyIncome(strcmp(caesar.FamilyIncome, '10k$-14.9k$')) = {'12500'};
caesar.FamilyIncome(strcmp(caesar.FamilyIncome, '15k$-19.9k$')) = {'17500'};
caesar.FamilyIncome(strcmp(caesar.FamilyIncome, '20k$-29.9k$')) = {'25000'};
caesar.FamilyIncome(strcmp(caesar.FamilyIncome, '30k$-44.9k$')) = {'37500'};
caesar.FamilyIncome(strcmp(caesar.FamilyIncome, '45k$-59.9k$')) = {'52500'};
caesar.FamilyIncome(strcmp(caesar.FamilyIncome, '60k$-79.9k$')) = {'70000'};
caesar.FamilyIncome(strcmp(caesar.FamilyIncome, '80k$-100k$')) = {'90000'};
caesar.FamilyIncome(strcmp(caesar.FamilyIncome, 'Over 100k$')) = {'150000'};
caesar.FamilyIncome(strcmp(caesar.FamilyIncome, 'Don''t Know')) = {''};
caesar.FamilyIncome(strcmp(caesar.FamilyIncome, 'No Response')) = {''};
caesar.FamilyIncome = str2double(caesar.FamilyIncome);

%% reformat site
% caesar.Site(strcmp(caesar.Site, 'Ames IA')) = {'Ames IA'};
caesar.Site(strcmp(caesar.Site, 'Dayton Oh')) = {'Dayton OH'};
caesar.Site(strcmp(caesar.Site, 'Dayton2 Oh')) = {'Dayton OH'};

%% reformat birth state
caesar.BirthState(strcmp(caesar.BirthState, 'Not Born in the US')) = {'Foreign'};
caesar.BirthState(strcmp(caesar.BirthState, 'US Territory')) = {'Foreign'};
 caesar.BirthState(strcmp(caesar.BirthState, 'No Response')) = {''};
caesar.BirthState(strcmp(caesar.BirthState, 'Alabama')) = {'South'};
caesar.BirthState(strcmp(caesar.BirthState, 'Alaska')) = {'West'};
caesar.BirthState(strcmp(caesar.BirthState, 'Arizona')) = {'West'};
caesar.BirthState(strcmp(caesar.BirthState, 'Arkansas')) = {'South'};
caesar.BirthState(strcmp(caesar.BirthState, 'California')) = {'West'};
caesar.BirthState(strcmp(caesar.BirthState, 'Colorado')) = {'West'};
caesar.BirthState(strcmp(caesar.BirthState, 'Connecticut')) = {'Northeast'};
caesar.BirthState(strcmp(caesar.BirthState, 'Florida')) = {'South'};
caesar.BirthState(strcmp(caesar.BirthState, 'Georgia')) = {'South'};
caesar.BirthState(strcmp(caesar.BirthState, 'Hawaii')) = {'West'};
caesar.BirthState(strcmp(caesar.BirthState, 'Idaho')) = {'West'};
caesar.BirthState(strcmp(caesar.BirthState, 'Illinois')) = {'Midwest'};
caesar.BirthState(strcmp(caesar.BirthState, 'Indiana')) = {'Midwest'};
caesar.BirthState(strcmp(caesar.BirthState, 'Iowa')) = {'Midwest'};
caesar.BirthState(strcmp(caesar.BirthState, 'Kansas')) = {'Midwest'};
caesar.BirthState(strcmp(caesar.BirthState, 'Kentucky')) = {'South'};
caesar.BirthState(strcmp(caesar.BirthState, 'Louisiana')) = {'South'};
caesar.BirthState(strcmp(caesar.BirthState, 'Maine')) = {'Northeast'};
caesar.BirthState(strcmp(caesar.BirthState, 'Maryland')) = {'South'};
caesar.BirthState(strcmp(caesar.BirthState, 'Massachusetts')) = {'Northeast'};
caesar.BirthState(strcmp(caesar.BirthState, 'Michigan')) = {'Midwest'};
caesar.BirthState(strcmp(caesar.BirthState, 'Minnesota')) = {'Midwest'};
caesar.BirthState(strcmp(caesar.BirthState, 'Mississippi')) = {'South'};
caesar.BirthState(strcmp(caesar.BirthState, 'Missouri')) = {'Midwest'};
caesar.BirthState(strcmp(caesar.BirthState, 'Montana')) = {'West'};
caesar.BirthState(strcmp(caesar.BirthState, 'Nebraska')) = {'Midwest'};
caesar.BirthState(strcmp(caesar.BirthState, 'Nevada')) = {'West'};
caesar.BirthState(strcmp(caesar.BirthState, 'New Hampshire')) = {'Northeast'};
caesar.BirthState(strcmp(caesar.BirthState, 'New Jersey')) = {'Northeast'};
caesar.BirthState(strcmp(caesar.BirthState, 'New Mexico')) = {'West'};
caesar.BirthState(strcmp(caesar.BirthState, 'New York')) = {'Northeast'};
caesar.BirthState(strcmp(caesar.BirthState, 'North Carolina')) = {'South'};
caesar.BirthState(strcmp(caesar.BirthState, 'North Dakota')) = {'Midwest'};
caesar.BirthState(strcmp(caesar.BirthState, 'Ohio')) = {'Midwest'};
caesar.BirthState(strcmp(caesar.BirthState, 'Oklahoma')) = {'South'};
caesar.BirthState(strcmp(caesar.BirthState, 'Oregon')) = {'West'};
caesar.BirthState(strcmp(caesar.BirthState, 'Pennsylvania')) = {'Northeast'};
caesar.BirthState(strcmp(caesar.BirthState, 'South Carolina')) = {'South'};
caesar.BirthState(strcmp(caesar.BirthState, 'South Dakota')) = {'Midwest'};
caesar.BirthState(strcmp(caesar.BirthState, 'Tennessee')) = {'South'};
caesar.BirthState(strcmp(caesar.BirthState, 'Texas')) = {'South'};
caesar.BirthState(strcmp(caesar.BirthState, 'Utah')) = {'West'};
caesar.BirthState(strcmp(caesar.BirthState, 'Vermont')) = {'Northeast'};
caesar.BirthState(strcmp(caesar.BirthState, 'Virginia')) = {'South'};
caesar.BirthState(strcmp(caesar.BirthState, 'Washington')) = {'West'};
caesar.BirthState(strcmp(caesar.BirthState, 'Washington DC')) = {'South'};
caesar.BirthState(strcmp(caesar.BirthState, 'West Virginia')) = {'South'};
caesar.BirthState(strcmp(caesar.BirthState, 'Wisconsin')) = {'Midwest'};
caesar.BirthState(strcmp(caesar.BirthState, 'Wyoming')) = {'West'};

%% reformat occupation
caesar.Occupation(strcmp(caesar.Occupation, 'No Response')) = {''};
caesar.Occupation(strcmp(caesar.Occupation, 'Homemaker')) = {''};
caesar.Occupation(strcmp(caesar.Occupation, 'Student')) = {''};
caesar.Occupation(strcmp(caesar.Occupation, 'Training/Continuing Education')) = {''};
caesar.Occupation(strcmp(caesar.Occupation, 'Unemployed')) = {''};
caesar.Occupation(strcmp(caesar.Occupation, 'Retired')) = {''};
caesar.Occupation(strcmp(caesar.Occupation, 'Armed Services')) = {''};

caesar.Occupation(strcmp(caesar.Occupation, 'Administative Support')) = {'White Collar'};
caesar.Occupation(strcmp(caesar.Occupation, 'Administrator')) = {'White Collar'};
caesar.Occupation(strcmp(caesar.Occupation, 'Attorney or Judge')) = {'White Collar'};
caesar.Occupation(strcmp(caesar.Occupation, 'Classroom Teacher')) = {'White Collar'};
caesar.Occupation(strcmp(caesar.Occupation, 'Computer Programmer/Software Engineer')) = {'White Collar'};
caesar.Occupation(strcmp(caesar.Occupation, 'Construction')) = {'Blue Collar'};
caesar.Occupation(strcmp(caesar.Occupation, 'Degreed Engineer')) = {'White Collar'};
caesar.Occupation(strcmp(caesar.Occupation, 'Farm Occupation')) = {'Blue Collar'};
caesar.Occupation(strcmp(caesar.Occupation, 'Forestry or Fishing Occupation')) = {'Blue Collar'};
caesar.Occupation(strcmp(caesar.Occupation, 'Health Diagnosing Occupation')) = {'White Collar'};
caesar.Occupation(strcmp(caesar.Occupation, 'Health Non-Diagnosing Occupation')) = {'White Collar'};
caesar.Occupation(strcmp(caesar.Occupation, 'Machine Operator')) = {'Blue Collar'};
caesar.Occupation(strcmp(caesar.Occupation, 'Management')) = {'Management'};
caesar.Occupation(strcmp(caesar.Occupation, 'Material Handler')) = {'Blue Collar'};
caesar.Occupation(strcmp(caesar.Occupation, 'Mechanic')) = {'Blue Collar'};
caesar.Occupation(strcmp(caesar.Occupation, 'Other Legal/Judicial Occupation')) = {'White Collar'};
caesar.Occupation(strcmp(caesar.Occupation, 'Other Specialty Occupation')) = {'White Collar'};
caesar.Occupation(strcmp(caesar.Occupation, 'Sales/Marketing')) = {'Service'};
caesar.Occupation(strcmp(caesar.Occupation, 'Scientist')) = {'White Collar'};
caesar.Occupation(strcmp(caesar.Occupation, 'Service Occupation')) = {'Service'};
caesar.Occupation(strcmp(caesar.Occupation, 'Supervisor')) = {'White Collar'};
caesar.Occupation(strcmp(caesar.Occupation, 'Technician')) = {'Blue Collar'};
caesar.Occupation(strcmp(caesar.Occupation, 'Transportation Occupation')) = {'Blue Collar'};

%% reformat education
caesar.Education(strcmp(caesar.Education, 'No Response')) = {''};
caesar.Education(strcmp(caesar.Education, 'None of the above')) = {''};
caesar.Education(strcmp(caesar.Education, 'High School')) = {'12'};
caesar.Education(strcmp(caesar.Education, 'Technical Training')) = {'14'};
caesar.Education(strcmp(caesar.Education, 'Some College')) = {'14'};
caesar.Education(strcmp(caesar.Education, 'Associates')) = {'14'};
caesar.Education(strcmp(caesar.Education, 'Bachelors')) = {'16'};
caesar.Education(strcmp(caesar.Education, 'Masters')) = {'18'};
caesar.Education(strcmp(caesar.Education, 'Doctorate')) = {'21'};
caesar.Education(strcmp(caesar.Education, 'Post-Doctoral Studies')) = {'24'};
caesar.Education = str2double(caesar.Education);

%% reformat number of children
caesar.NumberOfChildren(strcmp(caesar.NumberOfChildren, 'No Response')) = {''};
caesar.NumberOfChildren(strcmp(caesar.NumberOfChildren, '7 or more')) = {'7'};
caesar.NumberOfChildren = str2double(caesar.NumberOfChildren);

%% reformat fitness
caesar.Fitness(strcmp(caesar.Fitness, 'No Response')) = {''};
caesar.Fitness(strcmp(caesar.Fitness, '0-1')) = {'0.5'};
caesar.Fitness(strcmp(caesar.Fitness, '2-3')) = {'2.5'};
caesar.Fitness(strcmp(caesar.Fitness, '4-6')) = {'5'};
caesar.Fitness(strcmp(caesar.Fitness, '6-10')) = {'8'};
caesar.Fitness(strcmp(caesar.Fitness, 'More than 10')) = {'10'};
caesar.Fitness = str2double(caesar.Fitness);

%% reformat race
caesar.Race(strcmp(caesar.Race, 'No Respose')) = {''};
caesar.Race(strcmp(caesar.Race, 'Other Not Listed Above')) = {''};
caesar.Race(strcmp(caesar.Race, 'Other Mixed Race')) = {''};

caesar.Race(strcmp(caesar.Race, 'Native American/Alaskan')) = {''};

% caesar.Race(strcmp(caesar.Race, 'White')) = {'White'};

caesar.Race(strcmp(caesar.Race, 'African American')) = {'Black'};
caesar.Race(strcmp(caesar.Race, 'Asian/Pacific Islander Asian Indian')) = {'Asian'};
caesar.Race(strcmp(caesar.Race, 'Asian/Pacific Islander Chinese')) = {'Asian'};
caesar.Race(strcmp(caesar.Race, 'Asian/Pacific Islander Filipano')) = {'Asian'};
caesar.Race(strcmp(caesar.Race, 'Asian/Pacific Islander Japanese')) = {'Asian'};
caesar.Race(strcmp(caesar.Race, 'Asian/Pacific Islander Korean')) = {'Asian'};
caesar.Race(strcmp(caesar.Race, 'Asian/Pacific Islander Vietnamese')) = {'Asian'};
caesar.Race(strcmp(caesar.Race, 'Asian/Pacific Islander Other')) = {'Asian'};
caesar.Race(strcmp(caesar.Race, 'Spanish/Hispanic Cuban')) = {'Hispanic'};
caesar.Race(strcmp(caesar.Race, 'Spanish/Hispanic Mexican American')) = {'Hispanic'};
caesar.Race(strcmp(caesar.Race, 'Spanish/Hispanic Other')) = {'Hispanic'};
caesar.Race(strcmp(caesar.Race, 'Spanish/Hispanic Puerto Rican')) = {'Hispanic'};

% caesar.Race(strcmp(caesar.Race, 'African American')) = {'NonWhite'};
% caesar.Race(strcmp(caesar.Race, 'Asian/Pacific Islander Asian Indian')) = {'NonWhite'};
% caesar.Race(strcmp(caesar.Race, 'Asian/Pacific Islander Chinese')) = {'NonWhite'};
% caesar.Race(strcmp(caesar.Race, 'Asian/Pacific Islander Filipano')) = {'NonWhite'};
% caesar.Race(strcmp(caesar.Race, 'Asian/Pacific Islander Japanese')) = {'NonWhite'};
% caesar.Race(strcmp(caesar.Race, 'Asian/Pacific Islander Korean')) = {'NonWhite'};
% caesar.Race(strcmp(caesar.Race, 'Asian/Pacific Islander Vietnamese')) = {'NonWhite'};
% caesar.Race(strcmp(caesar.Race, 'Asian/Pacific Islander Other')) = {'NonWhite'};
% caesar.Race(strcmp(caesar.Race, 'Spanish/Hispanic Cuban')) = {'NonWhite'};
% caesar.Race(strcmp(caesar.Race, 'Spanish/Hispanic Mexican American')) = {'NonWhite'};
% caesar.Race(strcmp(caesar.Race, 'Spanish/Hispanic Other')) = {'NonWhite'};
% caesar.Race(strcmp(caesar.Race, 'Spanish/Hispanic Puerto Rican')) = {'NonWhite'};

%% reformat marital status
caesar.MaritalStatus(strcmp(caesar.MaritalStatus, 'No Response')) = {''};

caesar.MaritalStatus(strcmp(caesar.MaritalStatus, 'Divorced')) = {'Divorced/Widowed'};
caesar.MaritalStatus(strcmp(caesar.MaritalStatus, 'Engaged')) = {'Single'};
% caesar.MaritalStatus(strcmp(caesar.MaritalStatus, 'Married')) = {'Married'};
caesar.MaritalStatus(strcmp(caesar.MaritalStatus, 'Single')) = {'Single'};
caesar.MaritalStatus(strcmp(caesar.MaritalStatus, 'Widowed')) = {'Divorced/Widowed'};

% caesar.MaritalStatus(strcmp(caesar.MaritalStatus, 'No Response')) = {''};
% 
% caesar.MaritalStatus(strcmp(caesar.MaritalStatus, 'Divorced')) = {'Married'};
% caesar.MaritalStatus(strcmp(caesar.MaritalStatus, 'Engaged')) = {'NeverMarried'};
% % caesar.MaritalStatus(strcmp(caesar.MaritalStatus, 'Married')) = {'Married'};
% caesar.MaritalStatus(strcmp(caesar.MaritalStatus, 'Single')) = {'NeverMarried'};
% caesar.MaritalStatus(strcmp(caesar.MaritalStatus, 'Widowed')) = {'Married'};

%% reformat car make
caesar.CarMake(strcmp(caesar.CarMake, 'No Response')) = {''};
caesar.CarMake(strcmp(caesar.CarMake, 'Other')) = {''};

caesar.CarMake(strcmp(caesar.CarMake, 'Acura')) = {'Luxury'};
caesar.CarMake(strcmp(caesar.CarMake, 'Audi')) = {'Luxury'};
caesar.CarMake(strcmp(caesar.CarMake, 'BMW')) = {'Luxury'};
caesar.CarMake(strcmp(caesar.CarMake, 'Buick')) = {'Economy'};
caesar.CarMake(strcmp(caesar.CarMake, 'Cadilac')) = {'Luxury'};
caesar.CarMake(strcmp(caesar.CarMake, 'Chevrolet')) = {'Economy'};
caesar.CarMake(strcmp(caesar.CarMake, 'Chrysler')) = {'Luxury'};
caesar.CarMake(strcmp(caesar.CarMake, 'Dodge')) = {'Economy'};
caesar.CarMake(strcmp(caesar.CarMake, 'Eagle')) = {'Economy'};
caesar.CarMake(strcmp(caesar.CarMake, 'Ford')) = {'Economy'};
caesar.CarMake(strcmp(caesar.CarMake, 'GMC')) = {'Economy'};
caesar.CarMake(strcmp(caesar.CarMake, 'Honda')) = {'Economy'};
caesar.CarMake(strcmp(caesar.CarMake, 'Hyundai')) = {'Economy'};
caesar.CarMake(strcmp(caesar.CarMake, 'Infiniti')) = {'Luxury'};
caesar.CarMake(strcmp(caesar.CarMake, 'Isuzu')) = {'Economy'};
caesar.CarMake(strcmp(caesar.CarMake, 'Jeep')) = {'Economy'};
caesar.CarMake(strcmp(caesar.CarMake, 'Lexus')) = {'Luxury'};
caesar.CarMake(strcmp(caesar.CarMake, 'Lincoln')) = {'Luxury'};
caesar.CarMake(strcmp(caesar.CarMake, 'Mazda')) = {'Economy'};
caesar.CarMake(strcmp(caesar.CarMake, 'Mercedes-Benz')) = {'Luxury'};
caesar.CarMake(strcmp(caesar.CarMake, 'Mercury')) = {'Economy'};
caesar.CarMake(strcmp(caesar.CarMake, 'Mitsubishi')) = {'Economy'};
caesar.CarMake(strcmp(caesar.CarMake, 'Nissan')) = {'Economy'};
caesar.CarMake(strcmp(caesar.CarMake, 'Oldsmobile')) = {'Economy'};
caesar.CarMake(strcmp(caesar.CarMake, 'Plymouth')) = {'Economy'};
caesar.CarMake(strcmp(caesar.CarMake, 'Pontiac')) = {'Economy'};
caesar.CarMake(strcmp(caesar.CarMake, 'Porsche')) = {'Luxury'};
caesar.CarMake(strcmp(caesar.CarMake, 'Saab')) = {'Economy'};
caesar.CarMake(strcmp(caesar.CarMake, 'Saturn')) = {'Economy'};
caesar.CarMake(strcmp(caesar.CarMake, 'Subaru')) = {'Economy'};
caesar.CarMake(strcmp(caesar.CarMake, 'Suzuki')) = {'Economy'};
caesar.CarMake(strcmp(caesar.CarMake, 'Toyota')) = {'Economy'};
caesar.CarMake(strcmp(caesar.CarMake, 'Volkswagen')) = {'Economy'};
caesar.CarMake(strcmp(caesar.CarMake, 'Volvo')) = {'Luxury'};

%% reformat car model
% caesar.CarModel(strcmp(caesar.CarModel, 'No Response')) = {''};
% caesar.CarModel(strcmp(caesar.CarModel, 'Don''t Know')) = {''};
% caesar.CarModel(strcmp(caesar.CarModel, 'Other')) = {''};
% 
% caesar.CarModel(strcmp(caesar.CarModel, 'Compact')) = {'Compact/Economy'};
% caesar.CarModel(strcmp(caesar.CarModel, 'Economy')) = {'Compact/Economy'};
% caesar.CarModel(strcmp(caesar.CarModel, 'Full Size 4-Dr')) = {'Full Size'};
% caesar.CarModel(strcmp(caesar.CarModel, 'Full size 2-Dr')) = {'Full Size'};
% % caesar.CarModel(strcmp(caesar.CarModel, 'Intermediate')) = {'Intermediate'};
% caesar.CarModel(strcmp(caesar.CarModel, 'Luxury')) = {'Full Size'};
% caesar.CarModel(strcmp(caesar.CarModel, 'Minivan')) = {'SUV/Minivan'};
% caesar.CarModel(strcmp(caesar.CarModel, 'SUV')) = {'SUV/Minivan'};
% caesar.CarModel(strcmp(caesar.CarModel, 'Sports Car')) = {'Intermediate'};
% caesar.CarModel(strcmp(caesar.CarModel, 'Station Wagon')) = {'SUV/Minivan'};
% caesar.CarModel(strcmp(caesar.CarModel, 'Truck')) = {'Truck/Van'};
% caesar.CarModel(strcmp(caesar.CarModel, 'Van')) = {'Truck/Van'};

caesar.CarModel(strcmp(caesar.CarModel, 'No Response')) = {''};
caesar.CarModel(strcmp(caesar.CarModel, 'Don''t Know')) = {''};
caesar.CarModel(strcmp(caesar.CarModel, 'Other')) = {''};

caesar.CarModel(strcmp(caesar.CarModel, 'Compact')) = {'Sedan'};
caesar.CarModel(strcmp(caesar.CarModel, 'Economy')) = {'Sedan'};
caesar.CarModel(strcmp(caesar.CarModel, 'Full Size 4-Dr')) = {'Sedan'};
caesar.CarModel(strcmp(caesar.CarModel, 'Full size 2-Dr')) = {'Sedan'};
 caesar.CarModel(strcmp(caesar.CarModel, 'Intermediate')) = {'Sedan'};
caesar.CarModel(strcmp(caesar.CarModel, 'Luxury')) = {'Sedan'};
caesar.CarModel(strcmp(caesar.CarModel, 'Minivan')) = {'Non-sedan'};
caesar.CarModel(strcmp(caesar.CarModel, 'SUV')) = {'Non-sedan'};
caesar.CarModel(strcmp(caesar.CarModel, 'Sports Car')) = {'Sedan'};
caesar.CarModel(strcmp(caesar.CarModel, 'Station Wagon')) = {'Non-sedan'};
caesar.CarModel(strcmp(caesar.CarModel, 'Truck')) = {'Non-sedan'};
caesar.CarModel(strcmp(caesar.CarModel, 'Van')) = {'Non-sedan'};


%% reformat shoe size
caesar.ShoeSize(strcmp(caesar.ShoeSize, '14 or Larger')) = {'14'};
caesar.ShoeSize(strcmp(caesar.ShoeSize, '5 or Smaller')) = {'5'};
caesar.ShoeSize(strcmp(caesar.ShoeSize, 'Don''t Know')) = {''};
caesar.ShoeSize(strcmp(caesar.ShoeSize, 'No Response')) = {''};
caesar.ShoeSize = str2double(caesar.ShoeSize);

%% reformat jacket size
caesar.JacketSize(strcmp(caesar.JacketSize, '30 or Smaller')) = {'30'};
caesar.JacketSize(strcmp(caesar.JacketSize, '48 or Larger')) = {'48'};
caesar.JacketSize(strcmp(caesar.JacketSize, 'Don''t Know')) = {''};
caesar.JacketSize(strcmp(caesar.JacketSize, 'No Response')) = {''};
caesar.JacketSize = str2double(caesar.JacketSize);

%% reformat pants size waist
caesar.PantsSizeWaist(strcmp(caesar.PantsSizeWaist, '28 or Smaller')) = {'28'};
caesar.PantsSizeWaist(strcmp(caesar.PantsSizeWaist, '46 or Larger')) = {'46'};
caesar.PantsSizeWaist(strcmp(caesar.PantsSizeWaist, 'Don''t Know')) = {''};
caesar.PantsSizeWaist(strcmp(caesar.PantsSizeWaist, 'No Response')) = {''};
caesar.PantsSizeWaist = str2double(caesar.PantsSizeWaist);

%% reformat pants size inseam
caesar.PantsSizeInseam(strcmp(caesar.PantsSizeInseam, '28 or Smaller')) = {'28'};
caesar.PantsSizeInseam(strcmp(caesar.PantsSizeInseam, 'Don''t Know')) = {''};
caesar.PantsSizeInseam(strcmp(caesar.PantsSizeInseam, 'No Response')) = {''};
caesar.PantsSizeInseam = str2double(caesar.PantsSizeInseam);

%% reformat blouse size
caesar.BlouseSize(strcmp(caesar.BlouseSize, '4 or Smaller')) = {'4'};
caesar.BlouseSize(strcmp(caesar.BlouseSize, '22 or Larger')) = {'22'};
caesar.BlouseSize(strcmp(caesar.BlouseSize, 'Don''t Know')) = {''};
caesar.BlouseSize(strcmp(caesar.BlouseSize, 'No Response')) = {''};
caesar.BlouseSize = str2double(caesar.BlouseSize);

%% reformat pants size woman
caesar.PantsSizeWoman(strcmp(caesar.PantsSizeWoman, '2 or Smaller')) = {'2'};
caesar.PantsSizeWoman(strcmp(caesar.PantsSizeWoman, '20 or Larger')) = {'20'};
caesar.PantsSizeWoman(strcmp(caesar.PantsSizeWoman, 'Don''t Know')) = {''};
caesar.PantsSizeWoman(strcmp(caesar.PantsSizeWoman, 'No Response')) = {''};
caesar.PantsSizeWoman = str2double(caesar.PantsSizeWoman);

%% reformat bra size (TODO: split the cup size?)
% tbl = [caesar, array2table(BMI, 'VariableNames', {'BMI'})];
caesar = [caesar, array2table(repmat({''}, size(caesar,1), 1), 'VariableNames', {'CupSize'})];
caesar.CupSize(strcmp(caesar.BraSize, '32a')) = {'a'};
caesar.CupSize(strcmp(caesar.BraSize, '34a')) = {'a'};
caesar.CupSize(strcmp(caesar.BraSize, '36a')) = {'a'};
caesar.CupSize(strcmp(caesar.BraSize, '38a')) = {'a'};
caesar.CupSize(strcmp(caesar.BraSize, '40a')) = {'a'};
caesar.CupSize(strcmp(caesar.BraSize, '42a')) = {'a'};
caesar.CupSize(strcmp(caesar.BraSize, '44a')) = {'a'};
caesar.CupSize(strcmp(caesar.BraSize, '46a')) = {'a'};

caesar.CupSize(strcmp(caesar.BraSize, '32b')) = {'b'};
caesar.CupSize(strcmp(caesar.BraSize, '34b')) = {'b'};
caesar.CupSize(strcmp(caesar.BraSize, '36b')) = {'b'};
caesar.CupSize(strcmp(caesar.BraSize, '38b')) = {'b'};
caesar.CupSize(strcmp(caesar.BraSize, '40b')) = {'b'};
caesar.CupSize(strcmp(caesar.BraSize, '42b')) = {'b'};
caesar.CupSize(strcmp(caesar.BraSize, '44b')) = {'b'};
caesar.CupSize(strcmp(caesar.BraSize, '46b')) = {'b'};

caesar.CupSize(strcmp(caesar.BraSize, '32c')) = {'c'};
caesar.CupSize(strcmp(caesar.BraSize, '34c')) = {'c'};
caesar.CupSize(strcmp(caesar.BraSize, '36c')) = {'c'};
caesar.CupSize(strcmp(caesar.BraSize, '38c')) = {'c'};
caesar.CupSize(strcmp(caesar.BraSize, '40c')) = {'c'};
caesar.CupSize(strcmp(caesar.BraSize, '42c')) = {'c'};
caesar.CupSize(strcmp(caesar.BraSize, '44c')) = {'c'};
caesar.CupSize(strcmp(caesar.BraSize, '46c')) = {'c'};

caesar.CupSize(strcmp(caesar.BraSize, '32d')) = {'d'};
caesar.CupSize(strcmp(caesar.BraSize, '34d')) = {'d'};
caesar.CupSize(strcmp(caesar.BraSize, '36d')) = {'d'};
caesar.CupSize(strcmp(caesar.BraSize, '38d')) = {'d'};
caesar.CupSize(strcmp(caesar.BraSize, '40d')) = {'d'};
caesar.CupSize(strcmp(caesar.BraSize, '42d')) = {'d'};
caesar.CupSize(strcmp(caesar.BraSize, '44d')) = {'d'};
caesar.CupSize(strcmp(caesar.BraSize, '46d')) = {'d'};

caesar.CupSize(strcmp(caesar.BraSize, '32dd')) = {'dd'};
caesar.CupSize(strcmp(caesar.BraSize, '34dd')) = {'dd'};
caesar.CupSize(strcmp(caesar.BraSize, '36dd')) = {'dd'};
caesar.CupSize(strcmp(caesar.BraSize, '38dd')) = {'dd'};
caesar.CupSize(strcmp(caesar.BraSize, '40dd')) = {'dd'};
caesar.CupSize(strcmp(caesar.BraSize, '42dd')) = {'dd'};
caesar.CupSize(strcmp(caesar.BraSize, '44dd')) = {'dd'};
caesar.CupSize(strcmp(caesar.BraSize, '46dd')) = {'dd'};

caesar.CupSize(strcmp(caesar.BraSize, '32ddd')) = {'ddd'};
caesar.CupSize(strcmp(caesar.BraSize, '34ddd')) = {'ddd'};
caesar.CupSize(strcmp(caesar.BraSize, '36ddd')) = {'ddd'};
caesar.CupSize(strcmp(caesar.BraSize, '38ddd')) = {'ddd'};
caesar.CupSize(strcmp(caesar.BraSize, '40ddd')) = {'ddd'};
caesar.CupSize(strcmp(caesar.BraSize, '42ddd')) = {'ddd'};
caesar.CupSize(strcmp(caesar.BraSize, '44ddd')) = {'ddd'};
caesar.CupSize(strcmp(caesar.BraSize, '46ddd')) = {'ddd'};


caesar.BraSize(strcmp(caesar.BraSize, '32a')) = {'32'};
caesar.BraSize(strcmp(caesar.BraSize, '34a')) = {'34'};
caesar.BraSize(strcmp(caesar.BraSize, '36a')) = {'36'};
caesar.BraSize(strcmp(caesar.BraSize, '38a')) = {'38'};
caesar.BraSize(strcmp(caesar.BraSize, '40a')) = {'40'};
caesar.BraSize(strcmp(caesar.BraSize, '42a')) = {'42'};
caesar.BraSize(strcmp(caesar.BraSize, '44a')) = {'44'};
caesar.BraSize(strcmp(caesar.BraSize, '46a')) = {'46'};

caesar.BraSize(strcmp(caesar.BraSize, '32b')) = {'32'};
caesar.BraSize(strcmp(caesar.BraSize, '34b')) = {'34'};
caesar.BraSize(strcmp(caesar.BraSize, '36b')) = {'36'};
caesar.BraSize(strcmp(caesar.BraSize, '38b')) = {'38'};
caesar.BraSize(strcmp(caesar.BraSize, '40b')) = {'40'};
caesar.BraSize(strcmp(caesar.BraSize, '42b')) = {'42'};
caesar.BraSize(strcmp(caesar.BraSize, '44b')) = {'44'};
caesar.BraSize(strcmp(caesar.BraSize, '46b')) = {'46'};

caesar.BraSize(strcmp(caesar.BraSize, '32c')) = {'32'};
caesar.BraSize(strcmp(caesar.BraSize, '34c')) = {'34'};
caesar.BraSize(strcmp(caesar.BraSize, '36c')) = {'36'};
caesar.BraSize(strcmp(caesar.BraSize, '38c')) = {'38'};
caesar.BraSize(strcmp(caesar.BraSize, '40c')) = {'40'};
caesar.BraSize(strcmp(caesar.BraSize, '42c')) = {'42'};
caesar.BraSize(strcmp(caesar.BraSize, '44c')) = {'44'};
caesar.BraSize(strcmp(caesar.BraSize, '46c')) = {'46'};

caesar.BraSize(strcmp(caesar.BraSize, '32d')) = {'32'};
caesar.BraSize(strcmp(caesar.BraSize, '34d')) = {'34'};
caesar.BraSize(strcmp(caesar.BraSize, '36d')) = {'36'};
caesar.BraSize(strcmp(caesar.BraSize, '38d')) = {'38'};
caesar.BraSize(strcmp(caesar.BraSize, '40d')) = {'40'};
caesar.BraSize(strcmp(caesar.BraSize, '42d')) = {'42'};
caesar.BraSize(strcmp(caesar.BraSize, '44d')) = {'44'};
caesar.BraSize(strcmp(caesar.BraSize, '46d')) = {'46'};

caesar.BraSize(strcmp(caesar.BraSize, '32dd')) = {'32'};
caesar.BraSize(strcmp(caesar.BraSize, '34dd')) = {'34'};
caesar.BraSize(strcmp(caesar.BraSize, '36dd')) = {'36'};
caesar.BraSize(strcmp(caesar.BraSize, '38dd')) = {'38'};
caesar.BraSize(strcmp(caesar.BraSize, '40dd')) = {'40'};
caesar.BraSize(strcmp(caesar.BraSize, '42dd')) = {'42'};
caesar.BraSize(strcmp(caesar.BraSize, '44dd')) = {'44'};
caesar.BraSize(strcmp(caesar.BraSize, '46dd')) = {'46'};

caesar.BraSize(strcmp(caesar.BraSize, '32ddd')) = {'32'};
caesar.BraSize(strcmp(caesar.BraSize, '34ddd')) = {'34'};
caesar.BraSize(strcmp(caesar.BraSize, '36ddd')) = {'36'};
caesar.BraSize(strcmp(caesar.BraSize, '38ddd')) = {'38'};
caesar.BraSize(strcmp(caesar.BraSize, '40ddd')) = {'40'};
caesar.BraSize(strcmp(caesar.BraSize, '42ddd')) = {'42'};
caesar.BraSize(strcmp(caesar.BraSize, '44ddd')) = {'44'};
caesar.BraSize(strcmp(caesar.BraSize, '46ddd')) = {'46'};

caesar.BraSize(strcmp(caesar.BraSize, '30 or Smaller')) = {'30'};
caesar.BraSize(strcmp(caesar.BraSize, '48 or Larger')) = {'48'};
caesar.BraSize(strcmp(caesar.BraSize, 'Don''t Know')) = {''};
caesar.BraSize(strcmp(caesar.BraSize, 'No Response')) = {''};
caesar.BraSize = str2double(caesar.BraSize);

caesar.CupSize(strcmp(caesar.CupSize, 'a')) = {'1'};
caesar.CupSize(strcmp(caesar.CupSize, 'b')) = {'2'};
caesar.CupSize(strcmp(caesar.CupSize, 'c')) = {'3'};
caesar.CupSize(strcmp(caesar.CupSize, 'd')) = {'4'};
caesar.CupSize(strcmp(caesar.CupSize, 'dd')) = {'5'};
caesar.CupSize(strcmp(caesar.CupSize, 'ddd')) = {'6'};
caesar.CupSize = str2double(caesar.CupSize);