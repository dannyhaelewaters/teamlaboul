from Bio import Entrez
import pycountry
import pandas as pd
import numpy as np
import pycountry_convert as pc


class EpoUtils:
    """
    A utility class for working with eponyms and species data.
    """

    def __init__(self, email, db):
        """
        Initializes an instance of the EpoUtils class.

        :param name: The name of the user.
        :param email: The email of the user.
        :param db: The database to search for articles.
        """
        self.email = email
        self.db = db

    def search_species(self, species_name: str):
        if self.email:
            Entrez.email = self.email
        else:
            print("Email environment variable not set.")

        handle = Entrez.esearch(db="taxonomy", term=species_name)
        record = Entrez.read(handle)
        handle.close()

        return record

    def fetch_species_data(self, taxon_id: int):

        if self.email:
            Entrez.email = self.email
        else:
            print("Email environment variable not set.")

        handle = Entrez.efetch(db="taxonomy", id=taxon_id, retmode="xml")
        record = Entrez.read(handle)
        handle.close()
        return record

    def verify_species(self, species_name: str):
        """
        Verifies if a species name is valid.

        :param species_name: The name of the species to verify.
        :return: None
        """
        # let user know that we are searching
        print("searching....")
        search_result = self.search_species(species_name)

        # return the search result
        print("search result")
        print(search_result)

        if int(search_result["Count"]) > 0:
            taxon_id = search_result["IdList"][
                0
            ]  # Get the taxonomic ID of the first result
            species_data = self.fetch_species_data(taxon_id)
            print(
                species_data, width=5, compact=True
            )  # Set the width parameter to control the text wrapping
            print(f"Verified species: {species_name}", width=5, compact=True)
        else:
            print("Species not found.", width=5, compact=True)

    def search_articles(self, species_name: str, max_results: int = 5):

        # email passed in
        email = self.email

        # Check if email is set, if not, print a message
        if email:
            Entrez.email = email
        else:
            print(
                "A valid email associated with NCBI was not provided. Please try again"
            )

        # Search for articles
        print(f"Searching for articles on {species_name}...")
        handle = Entrez.esearch(db=self.db, term=species_name, retmax=max_results)

        # Read the search results
        record = Entrez.read(handle)

        # Close the handle
        handle.close()

        return record

    def return_articles(self, species_name: str, results: int = 5):
        """
        Returns the article IDs for a given species name

        :param species_name: The name of the species to search for
        :param max: The maximum number of articles to return
        :return: None

        """

        # Search for articles
        search_result = self.search_articles(
            species_name, max_results=results, db=self.db
        )

        # Check if any articles were found, count greater than 0
        if int(search_result["Count"]) > 0:
            print(f"Found {search_result['Count']} articles on {species_name}.")
            print("printing...")

            article_ids = search_result["IdList"]  # Get the list of article IDs

            # Print the article IDs
            for article_id in article_ids:
                print("Article ID:", article_id)

            for identifier in search_result["IdList"]:
                pubmed_entry = Entrez.efetch(db=self.db, id=identifier, retmode="xml")
                result = Entrez.read(pubmed_entry)
                article = result["PubmedArticle"][0]["MedlineCitation"]["Article"]

                print(
                    identifier,
                    '"{}" in "{}"'.format(
                        article["ArticleTitle"], article["Journal"]["Title"]
                    ),
                )

            # close the entry
            pubmed_entry.close()

        # If no articles were found, return no articles found
        else:
            print("No articles found.")

    def find_geography(self, location: str):
        """
        Finds the geographic location based on the given location name.

        Args:
            location (str): The name of the location to search for.

        Returns:
            None

        Raises:
            AttributeError: If the given location is not a valid geographic location.

        """

        @staticmethod
        def is_continent(name):
            """checks if a continent"""
            continents = [
                "Africa",
                "Antarctica",
                "Asia",
                "Europe",
                "North America",
                "Oceania",
                "South America",
            ]

            if name in continents:
                return True
            else:
                return False

        try:
            if is_continent(
                location
            ):  # if this lcoation is a continent, we can determine that it's name off of it

                # it is valid
                print(f"{location} is a valid geographic location.")

                return location

            else:
                # attempts to get a country name
                location = pycountry.countries.get(name=location)

                # it is valid
                print(f"{location} is a valid geographic location.")

                return location

        except AttributeError:
            print(f"{location} is not a valid geographic location.")
            return None
