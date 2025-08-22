from models import Base, Batfly, Bats, BatSpecies
from connect import engine

# create the engine
print('Creating tables...')
Base.metadata.create_all(bind = engine)