from sqlalchemy import create_engine, text

engine = create_engine('mysql+pymysql://username:password@localhost/db_name', echo=True)

with engine.connect() as connection:
    result = connection.execute(text('SELECT * FROM sqlite_master WHERE type="table"'))

    print(result.fetchall())

    