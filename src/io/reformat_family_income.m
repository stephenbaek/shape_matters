function table = reformat_family_income(table)
    table = rename_items(table, 'FamilyIncome', 'Less than 10k$/10MLira', '7500');
    table = rename_items(table, 'FamilyIncome', '10k$-14.9k$', '12500');
    table = rename_items(table, 'FamilyIncome', '15k$-19.9k$', '17500');
    table = rename_items(table, 'FamilyIncome', '20k$-29.9k$', '25000');
    table = rename_items(table, 'FamilyIncome', '30k$-44.9k$', '37500');
    table = rename_items(table, 'FamilyIncome', '45k$-59.9k$', '52500');
    table = rename_items(table, 'FamilyIncome', '60k$-79.9k$', '70000');
    table = rename_items(table, 'FamilyIncome', '80k$-100k$', '90000');
    table = rename_items(table, 'FamilyIncome', 'Over 100k$', '150000');
    table = rename_items(table, 'FamilyIncome', 'Don''t Know', '');
    table = rename_items(table, 'FamilyIncome', 'No Response', '');
    table.FamilyIncome = str2double(table.FamilyIncome);
end