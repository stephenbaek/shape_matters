function table = reformat_fitness(table)
    table = rename_items(table, 'Fitness', 'No Response', '');
    table = rename_items(table, 'Fitness', '0-1', '0.5');
    table = rename_items(table, 'Fitness', '2-3', '2.5');
    table = rename_items(table, 'Fitness', '4-6', '5');
    table = rename_items(table, 'Fitness', '6-10', '8');
    table = rename_items(table, 'Fitness', 'More than 10', '10');
    table.Fitness = str2double(table.Fitness);
end