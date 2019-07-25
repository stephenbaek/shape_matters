function table = reformat_number_of_children(table)
    table = rename_items(table, 'NumberOfChildren', 'No Response', '');
    table = rename_items(table, 'NumberOfChildren', '7 or more', '7');
    table.NumberOfChildren = str2double(table.NumberOfChildren);
end