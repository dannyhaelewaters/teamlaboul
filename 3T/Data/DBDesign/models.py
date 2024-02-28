from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column
from sqlalchemy import ForeignKey, Text, Integer, DateTime, String, Column, create_engine
from sqlalchemy.sql import func
import datetime

class Base(DeclarativeBase):
    pass # gives us a base class to inherit from

class Bats(Base):
    __tablename__ = 'Bats'

    bat_id:Mapped[int] = mapped_column(primary_key=True)
    name:Mapped[str] = mapped_column(nullable = False)
    description:Mapped[str] = mapped_column()
    created_at:Mapped[datetime.datetime] = mapped_column(DateTime(timezone=True), server_default=func.now())
    updated_at:Mapped[datetime.datetime] = mapped_column(DateTime(timezone=True), onupdate=func.now())
    bat_species_id:Mapped[int] = mapped_column(ForeignKey('BatSpecies.bat_species_id'))

class BatSpecies(Base):
    __tablename__ = 'BatSpecies'

    bat_species_id:Mapped[int] = mapped_column(primary_key=True, nullable = False)
    name:Mapped[str] = mapped_column(nullable = False)

class Batfly(Base):
    __tablename__ = 'Batfly'

    batfly_id:Mapped[int] = mapped_column(primary_key=True)
    name:Mapped[str] = mapped_column(nullable = False)
    description:Mapped[str] = mapped_column()
    created_at:Mapped[datetime.datetime] = mapped_column(DateTime(timezone=True), server_default=func.now())
    updated_at:Mapped[datetime.datetime] = mapped_column(DateTime(timezone=True), onupdate=func.now())
