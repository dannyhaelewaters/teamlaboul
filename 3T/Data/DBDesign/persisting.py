from models import Base, Batfly, Bats, BatSpecies
import pandas as pd
import numpy as np
import csv
from connect import engine
from sqlalchemy.orm import Session

# need to iterate through the csv and insert the data into a database

# create the engine
print('Inserting rows...')

# create the session
session = Base.metadata.create_all(bind = engine)

# create the session
session = Session(bind = engine)

# read the csv
df = pd.read_csv('bats.csv')

# iterate through the csv and insert the data into the database

for index, row in df.iterrows():
    bat = Bats(name = row['name'], description = row['description'], bat_species_id = row['bat_species_id'])
    session.add(bat)
    
    bat_species = BatSpecies(name = row['name'], description = row['description'])
    session.add(bat_species)
    
    batfly = Batfly(name = row['name'], description = row['description'])
    session.add(batfly)
    
    session.commit()