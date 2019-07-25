function table = reformat_bra_size(table)
    table = [table, array2table(repmat({''}, size(table,1), 1), 'VariableNames', {'CupSize'})];
    table.CupSize(strcmp(table.BraSize, '32a')) = {'a'};
    table.CupSize(strcmp(table.BraSize, '34a')) = {'a'};
    table.CupSize(strcmp(table.BraSize, '36a')) = {'a'};
    table.CupSize(strcmp(table.BraSize, '38a')) = {'a'};
    table.CupSize(strcmp(table.BraSize, '40a')) = {'a'};
    table.CupSize(strcmp(table.BraSize, '42a')) = {'a'};
    table.CupSize(strcmp(table.BraSize, '44a')) = {'a'};
    table.CupSize(strcmp(table.BraSize, '46a')) = {'a'};

    table.CupSize(strcmp(table.BraSize, '32b')) = {'b'};
    table.CupSize(strcmp(table.BraSize, '34b')) = {'b'};
    table.CupSize(strcmp(table.BraSize, '36b')) = {'b'};
    table.CupSize(strcmp(table.BraSize, '38b')) = {'b'};
    table.CupSize(strcmp(table.BraSize, '40b')) = {'b'};
    table.CupSize(strcmp(table.BraSize, '42b')) = {'b'};
    table.CupSize(strcmp(table.BraSize, '44b')) = {'b'};
    table.CupSize(strcmp(table.BraSize, '46b')) = {'b'};

    table.CupSize(strcmp(table.BraSize, '32c')) = {'c'};
    table.CupSize(strcmp(table.BraSize, '34c')) = {'c'};
    table.CupSize(strcmp(table.BraSize, '36c')) = {'c'};
    table.CupSize(strcmp(table.BraSize, '38c')) = {'c'};
    table.CupSize(strcmp(table.BraSize, '40c')) = {'c'};
    table.CupSize(strcmp(table.BraSize, '42c')) = {'c'};
    table.CupSize(strcmp(table.BraSize, '44c')) = {'c'};
    table.CupSize(strcmp(table.BraSize, '46c')) = {'c'};

    table.CupSize(strcmp(table.BraSize, '32d')) = {'d'};
    table.CupSize(strcmp(table.BraSize, '34d')) = {'d'};
    table.CupSize(strcmp(table.BraSize, '36d')) = {'d'};
    table.CupSize(strcmp(table.BraSize, '38d')) = {'d'};
    table.CupSize(strcmp(table.BraSize, '40d')) = {'d'};
    table.CupSize(strcmp(table.BraSize, '42d')) = {'d'};
    table.CupSize(strcmp(table.BraSize, '44d')) = {'d'};
    table.CupSize(strcmp(table.BraSize, '46d')) = {'d'};

    table.CupSize(strcmp(table.BraSize, '32dd')) = {'dd'};
    table.CupSize(strcmp(table.BraSize, '34dd')) = {'dd'};
    table.CupSize(strcmp(table.BraSize, '36dd')) = {'dd'};
    table.CupSize(strcmp(table.BraSize, '38dd')) = {'dd'};
    table.CupSize(strcmp(table.BraSize, '40dd')) = {'dd'};
    table.CupSize(strcmp(table.BraSize, '42dd')) = {'dd'};
    table.CupSize(strcmp(table.BraSize, '44dd')) = {'dd'};
    table.CupSize(strcmp(table.BraSize, '46dd')) = {'dd'};

    table.CupSize(strcmp(table.BraSize, '32ddd')) = {'ddd'};
    table.CupSize(strcmp(table.BraSize, '34ddd')) = {'ddd'};
    table.CupSize(strcmp(table.BraSize, '36ddd')) = {'ddd'};
    table.CupSize(strcmp(table.BraSize, '38ddd')) = {'ddd'};
    table.CupSize(strcmp(table.BraSize, '40ddd')) = {'ddd'};
    table.CupSize(strcmp(table.BraSize, '42ddd')) = {'ddd'};
    table.CupSize(strcmp(table.BraSize, '44ddd')) = {'ddd'};
    table.CupSize(strcmp(table.BraSize, '46ddd')) = {'ddd'};


    table.BraSize(strcmp(table.BraSize, '32a')) = {'32'};
    table.BraSize(strcmp(table.BraSize, '34a')) = {'34'};
    table.BraSize(strcmp(table.BraSize, '36a')) = {'36'};
    table.BraSize(strcmp(table.BraSize, '38a')) = {'38'};
    table.BraSize(strcmp(table.BraSize, '40a')) = {'40'};
    table.BraSize(strcmp(table.BraSize, '42a')) = {'42'};
    table.BraSize(strcmp(table.BraSize, '44a')) = {'44'};
    table.BraSize(strcmp(table.BraSize, '46a')) = {'46'};

    table.BraSize(strcmp(table.BraSize, '32b')) = {'32'};
    table.BraSize(strcmp(table.BraSize, '34b')) = {'34'};
    table.BraSize(strcmp(table.BraSize, '36b')) = {'36'};
    table.BraSize(strcmp(table.BraSize, '38b')) = {'38'};
    table.BraSize(strcmp(table.BraSize, '40b')) = {'40'};
    table.BraSize(strcmp(table.BraSize, '42b')) = {'42'};
    table.BraSize(strcmp(table.BraSize, '44b')) = {'44'};
    table.BraSize(strcmp(table.BraSize, '46b')) = {'46'};

    table.BraSize(strcmp(table.BraSize, '32c')) = {'32'};
    table.BraSize(strcmp(table.BraSize, '34c')) = {'34'};
    table.BraSize(strcmp(table.BraSize, '36c')) = {'36'};
    table.BraSize(strcmp(table.BraSize, '38c')) = {'38'};
    table.BraSize(strcmp(table.BraSize, '40c')) = {'40'};
    table.BraSize(strcmp(table.BraSize, '42c')) = {'42'};
    table.BraSize(strcmp(table.BraSize, '44c')) = {'44'};
    table.BraSize(strcmp(table.BraSize, '46c')) = {'46'};

    table.BraSize(strcmp(table.BraSize, '32d')) = {'32'};
    table.BraSize(strcmp(table.BraSize, '34d')) = {'34'};
    table.BraSize(strcmp(table.BraSize, '36d')) = {'36'};
    table.BraSize(strcmp(table.BraSize, '38d')) = {'38'};
    table.BraSize(strcmp(table.BraSize, '40d')) = {'40'};
    table.BraSize(strcmp(table.BraSize, '42d')) = {'42'};
    table.BraSize(strcmp(table.BraSize, '44d')) = {'44'};
    table.BraSize(strcmp(table.BraSize, '46d')) = {'46'};

    table.BraSize(strcmp(table.BraSize, '32dd')) = {'32'};
    table.BraSize(strcmp(table.BraSize, '34dd')) = {'34'};
    table.BraSize(strcmp(table.BraSize, '36dd')) = {'36'};
    table.BraSize(strcmp(table.BraSize, '38dd')) = {'38'};
    table.BraSize(strcmp(table.BraSize, '40dd')) = {'40'};
    table.BraSize(strcmp(table.BraSize, '42dd')) = {'42'};
    table.BraSize(strcmp(table.BraSize, '44dd')) = {'44'};
    table.BraSize(strcmp(table.BraSize, '46dd')) = {'46'};

    table.BraSize(strcmp(table.BraSize, '32ddd')) = {'32'};
    table.BraSize(strcmp(table.BraSize, '34ddd')) = {'34'};
    table.BraSize(strcmp(table.BraSize, '36ddd')) = {'36'};
    table.BraSize(strcmp(table.BraSize, '38ddd')) = {'38'};
    table.BraSize(strcmp(table.BraSize, '40ddd')) = {'40'};
    table.BraSize(strcmp(table.BraSize, '42ddd')) = {'42'};
    table.BraSize(strcmp(table.BraSize, '44ddd')) = {'44'};
    table.BraSize(strcmp(table.BraSize, '46ddd')) = {'46'};

    table.BraSize(strcmp(table.BraSize, '30 or Smaller')) = {'30'};
    table.BraSize(strcmp(table.BraSize, '48 or Larger')) = {'48'};
    table.BraSize(strcmp(table.BraSize, 'Don''t Know')) = {''};
    table.BraSize(strcmp(table.BraSize, 'No Response')) = {''};
    table.BraSize = str2double(table.BraSize);

    table.CupSize(strcmp(table.CupSize, 'a')) = {'1'};
    table.CupSize(strcmp(table.CupSize, 'b')) = {'2'};
    table.CupSize(strcmp(table.CupSize, 'c')) = {'3'};
    table.CupSize(strcmp(table.CupSize, 'd')) = {'4'};
    table.CupSize(strcmp(table.CupSize, 'dd')) = {'5'};
    table.CupSize(strcmp(table.CupSize, 'ddd')) = {'6'};
    table.CupSize = str2double(table.CupSize);
end