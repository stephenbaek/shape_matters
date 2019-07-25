function table = reformat_marital_status(table)
    table = rename_items(table, 'MaritalStatus', 'No Response', '');

    table = rename_items(table, 'MaritalStatus', 'Divorced', 'Divorced/Widowed');
    table = rename_items(table, 'MaritalStatus', 'Engaged', 'Single');
    table = rename_items(table, 'MaritalStatus', 'Single', 'Single');
    table = rename_items(table, 'MaritalStatus', 'Widowed', 'Divorced/Widowed');
end