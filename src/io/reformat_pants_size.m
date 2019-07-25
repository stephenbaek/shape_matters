function table = reformat_pants_size(table)
    table = rename_items(table, 'PantsSizeWaist', '28 or Smaller', '28');
    table = rename_items(table, 'PantsSizeWaist', '46 or Larger', '46');
    table = rename_items(table, 'PantsSizeWaist', 'Don''t Know', '');
    table = rename_items(table, 'PantsSizeWaist', 'No Response', '');
    table.PantsSizeWaist = str2double(table.PantsSizeWaist);
    
    table = rename_items(table, 'PantsSizeInseam', '28 or Smaller', '28');
    table = rename_items(table, 'PantsSizeInseam', 'Don''t Know', '');
    table = rename_items(table, 'PantsSizeInseam', 'No Response', '');
    table.PantsSizeInseam = str2double(table.PantsSizeInseam);
    
    table = rename_items(table, 'PantsSizeWoman', '2 or Smaller', '2');
    table = rename_items(table, 'PantsSizeWoman', '20 or Larger', '20');
    table = rename_items(table, 'PantsSizeWoman', 'Don''t Know', '');
    table = rename_items(table, 'PantsSizeWoman', 'No Response', '');
    table.PantsSizeWoman = str2double(table.PantsSizeWoman);
end