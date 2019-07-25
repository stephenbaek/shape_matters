function table = reformat_race(table)
    table = rename_items(table, 'Race', 'No Respose', '');
    table = rename_items(table, 'Race', 'Other Not Listed Above', '');
    table = rename_items(table, 'Race', 'Other Mixed Race', '');

    table = rename_items(table, 'Race', 'Native American/Alaskan', '');

    table = rename_items(table, 'Race', 'African American', 'Black');

    table = rename_items(table, 'Race', 'Asian/Pacific Islander Asian Indian', 'Asian');
    table = rename_items(table, 'Race', 'Asian/Pacific Islander Chinese', 'Asian');
    table = rename_items(table, 'Race', 'Asian/Pacific Islander Filipano', 'Asian');
    table = rename_items(table, 'Race', 'Asian/Pacific Islander Japanese', 'Asian');
    table = rename_items(table, 'Race', 'Asian/Pacific Islander Korean', 'Asian');
    table = rename_items(table, 'Race', 'Asian/Pacific Islander Vietnamese', 'Asian');
    table = rename_items(table, 'Race', 'Asian/Pacific Islander Other', 'Asian');
    
    table = rename_items(table, 'Race', 'Spanish/Hispanic Cuban', 'Hispanic');
    table = rename_items(table, 'Race', 'Spanish/Hispanic Mexican American', 'Hispanic');
    table = rename_items(table, 'Race', 'Spanish/Hispanic Other', 'Hispanic');
    table = rename_items(table, 'Race', 'Spanish/Hispanic Puerto Rican', 'Hispanic');
end