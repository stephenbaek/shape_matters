function table = reformat_jacket_size(table)
    table = rename_items(table, 'JacketSize', '30 or Smaller', '30');
    table = rename_items(table, 'JacketSize', '48 or Larger', '48');
    table = rename_items(table, 'JacketSize', 'Don''t Know', '');
    table = rename_items(table, 'JacketSize', 'No Response', '');
    table.JacketSize = str2double(table.JacketSize);
    
    table = rename_items(table, 'BlouseSize', '4 or Smaller', '4');
    table = rename_items(table, 'BlouseSize', '22 or Larger', '22');
    table = rename_items(table, 'BlouseSize', 'Don''t Know', '');
    table = rename_items(table, 'BlouseSize', 'No Response', '');
    table.BlouseSize = str2double(table.BlouseSize);
end