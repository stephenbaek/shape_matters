function table = reformat_shoe_size(table)
table = rename_items(table, 'ShoeSize', '14 or Larger', '14');
table = rename_items(table, 'ShoeSize', '5 or Smaller', '5');
table = rename_items(table, 'ShoeSize', 'Don''t Know', '');
table = rename_items(table, 'ShoeSize', 'No Response', '');
table.ShoeSize = str2double(table.ShoeSize);
end