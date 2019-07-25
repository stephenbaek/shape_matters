function table = reformat_car_model(table)
    table = rename_items(table, 'CarModel', 'No Response', '');
    table = rename_items(table, 'CarModel', 'Don''t Know', '');
    table = rename_items(table, 'CarModel', 'Other', '');

    table = rename_items(table, 'CarModel', 'Compact', 'Sedan');
    table = rename_items(table, 'CarModel', 'Economy', 'Sedan');
    table = rename_items(table, 'CarModel', 'Full Size 4-Dr', 'Sedan');
    table = rename_items(table, 'CarModel', 'Full size 2-Dr', 'Sedan');
     table = rename_items(table, 'CarModel', 'Intermediate', 'Sedan');
    table = rename_items(table, 'CarModel', 'Luxury', 'Sedan');
    table = rename_items(table, 'CarModel', 'Sports Car', 'Sedan');
    
    table = rename_items(table, 'CarModel', 'Minivan', 'Non-sedan');
    table = rename_items(table, 'CarModel', 'SUV', 'Non-sedan');
    table = rename_items(table, 'CarModel', 'Station Wagon', 'Non-sedan');
    table = rename_items(table, 'CarModel', 'Truck', 'Non-sedan');
    table = rename_items(table, 'CarModel', 'Van', 'Non-sedan');
end