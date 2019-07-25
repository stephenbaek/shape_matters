function table = reformat_birth_state(table)
    table = rename_items(table, 'BirthState', 'No Response', '');
    
    table = rename_items(table, 'BirthState', 'Not Born in the US', 'Foreign');
    table = rename_items(table, 'BirthState', 'US Territory', 'Foreign');
    
    table = rename_items(table, 'BirthState', 'Illinois', 'Midwest');
    table = rename_items(table, 'BirthState', 'Indiana', 'Midwest');
    table = rename_items(table, 'BirthState', 'Iowa', 'Midwest');
    table = rename_items(table, 'BirthState', 'Kansas', 'Midwest');
    table = rename_items(table, 'BirthState', 'Michigan', 'Midwest');
    table = rename_items(table, 'BirthState', 'Minnesota', 'Midwest');
    table = rename_items(table, 'BirthState', 'Missouri', 'Midwest');
    table = rename_items(table, 'BirthState', 'Nebraska', 'Midwest');
    table = rename_items(table, 'BirthState', 'North Dakota', 'Midwest');
    table = rename_items(table, 'BirthState', 'Ohio', 'Midwest');
    table = rename_items(table, 'BirthState', 'South Dakota', 'Midwest');
    table = rename_items(table, 'BirthState', 'Wisconsin', 'Midwest');
        
    table = rename_items(table, 'BirthState', 'Alaska', 'West');
    table = rename_items(table, 'BirthState', 'Arizona', 'West');
    table = rename_items(table, 'BirthState', 'California', 'West');
    table = rename_items(table, 'BirthState', 'Colorado', 'West');
    table = rename_items(table, 'BirthState', 'Hawaii', 'West');
    table = rename_items(table, 'BirthState', 'Idaho', 'West');
    table = rename_items(table, 'BirthState', 'Montana', 'West');
    table = rename_items(table, 'BirthState', 'Nevada', 'West');
    table = rename_items(table, 'BirthState', 'New Mexico', 'West');
    table = rename_items(table, 'BirthState', 'Oregon', 'West');
    table = rename_items(table, 'BirthState', 'Utah', 'West');
    table = rename_items(table, 'BirthState', 'Washington', 'West');
    table = rename_items(table, 'BirthState', 'Wyoming', 'West');
    
    table = rename_items(table, 'BirthState', 'Connecticut', 'Northeast');
    table = rename_items(table, 'BirthState', 'Maine', 'Northeast');
    table = rename_items(table, 'BirthState', 'Massachusetts', 'Northeast');
    table = rename_items(table, 'BirthState', 'New Hampshire', 'Northeast');
    table = rename_items(table, 'BirthState', 'New Jersey', 'Northeast');
    table = rename_items(table, 'BirthState', 'New York', 'Northeast');
    table = rename_items(table, 'BirthState', 'Pennsylvania', 'Northeast');
    table = rename_items(table, 'BirthState', 'Vermont', 'Northeast');
    
    table = rename_items(table, 'BirthState', 'Alabama', 'South');
    table = rename_items(table, 'BirthState', 'Arkansas', 'South');
    table = rename_items(table, 'BirthState', 'Florida', 'South');
    table = rename_items(table, 'BirthState', 'Georgia', 'South');
    table = rename_items(table, 'BirthState', 'Kentucky', 'South');
    table = rename_items(table, 'BirthState', 'Louisiana', 'South');
    table = rename_items(table, 'BirthState', 'Maryland', 'South');
    table = rename_items(table, 'BirthState', 'Mississippi', 'South');    
    table = rename_items(table, 'BirthState', 'North Carolina', 'South');
    table = rename_items(table, 'BirthState', 'Oklahoma', 'South');    
    table = rename_items(table, 'BirthState', 'South Carolina', 'South');
    table = rename_items(table, 'BirthState', 'Tennessee', 'South');
    table = rename_items(table, 'BirthState', 'Texas', 'South');
    table = rename_items(table, 'BirthState', 'Virginia', 'South');
    table = rename_items(table, 'BirthState', 'Washington DC', 'South');
    table = rename_items(table, 'BirthState', 'West Virginia', 'South');
end
