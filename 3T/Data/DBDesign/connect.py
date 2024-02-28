from sqlalchemy import create_engine, text

#engine = create_engine('mysql+pymysql://username:password@localhost/db_name', echo=True)

engine = create_engine('sqlite:///sample.db', echo=True)

with engine.connect() as connection:
    result = connection.execute(text('select "Hello, World!"'))

    print(result.fetchall())

    