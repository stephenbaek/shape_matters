function table = rename_items(table, FieldName, from, to)
% Finds field items that matches with from string and replace them with to
    table.(FieldName)(strcmp(table.(FieldName), from)) = {to};
end