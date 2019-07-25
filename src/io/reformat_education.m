function table = reformat_education(table)
    table = rename_items(table, 'Education', 'No Response', '');
    table = rename_items(table, 'Education', 'None of the above', '');
    table = rename_items(table, 'Education', 'High School', '12');
    table = rename_items(table, 'Education', 'Technical Training', '14');
    table = rename_items(table, 'Education', 'Some College', '14');
    table = rename_items(table, 'Education', 'Associates', '14');
    table = rename_items(table, 'Education', 'Bachelors', '16');
    table = rename_items(table, 'Education', 'Masters', '18');
    table = rename_items(table, 'Education', 'Doctorate', '21');
    table = rename_items(table, 'Education', 'Post-Doctoral Studies', '24');
    table.Education = str2double(table.Education);
end