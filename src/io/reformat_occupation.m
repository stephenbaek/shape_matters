function table = reformat_occupation(table)
    table = rename_items(table, 'Occupation', 'No Response', '');
    table = rename_items(table, 'Occupation', 'Homemaker', '');
    table = rename_items(table, 'Occupation', 'Student', '');
    table = rename_items(table, 'Occupation', 'Training/Continuing Education', '');
    table = rename_items(table, 'Occupation', 'Unemployed', '');
    table = rename_items(table, 'Occupation', 'Retired', '');
    table = rename_items(table, 'Occupation', 'Armed Services', '');

    table = rename_items(table, 'Occupation', 'Management', 'Management');
    
    table = rename_items(table, 'Occupation', 'Administative Support', 'White Collar');
    table = rename_items(table, 'Occupation', 'Administrator', 'White Collar');
    table = rename_items(table, 'Occupation', 'Attorney or Judge', 'White Collar');
    table = rename_items(table, 'Occupation', 'Classroom Teacher', 'White Collar');
    table = rename_items(table, 'Occupation', 'Computer Programmer/Software Engineer', 'White Collar');
    table = rename_items(table, 'Occupation', 'Degreed Engineer', 'White Collar');
    table = rename_items(table, 'Occupation', 'Health Diagnosing Occupation', 'White Collar');
    table = rename_items(table, 'Occupation', 'Health Non-Diagnosing Occupation', 'White Collar');
    table = rename_items(table, 'Occupation', 'Other Legal/Judicial Occupation', 'White Collar');
    table = rename_items(table, 'Occupation', 'Other Specialty Occupation', 'White Collar');
    table = rename_items(table, 'Occupation', 'Scientist', 'White Collar');
    table = rename_items(table, 'Occupation', 'Supervisor', 'White Collar');
    
    table = rename_items(table, 'Occupation', 'Construction', 'Blue Collar');
    table = rename_items(table, 'Occupation', 'Farm Occupation', 'Blue Collar');
    table = rename_items(table, 'Occupation', 'Forestry or Fishing Occupation', 'Blue Collar');
    table = rename_items(table, 'Occupation', 'Machine Operator', 'Blue Collar');
    table = rename_items(table, 'Occupation', 'Material Handler', 'Blue Collar');
    table = rename_items(table, 'Occupation', 'Mechanic', 'Blue Collar');
    table = rename_items(table, 'Occupation', 'Technician', 'Blue Collar');
    table = rename_items(table, 'Occupation', 'Transportation Occupation', 'Blue Collar');
    
    table = rename_items(table, 'Occupation', 'Sales/Marketing', 'Service');
    table = rename_items(table, 'Occupation', 'Service Occupation', 'Service');
end