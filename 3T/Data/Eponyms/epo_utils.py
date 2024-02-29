from Bio import Entrez
import pycountry
import pandas as pd
import numpy as np
import pycountry_convert as pc


class EpoUtils:
    """
    A utility class for working with eponyms and species data.
    """

    def __init__(self, name, email, db):
        """
        Initializes an instance of the EpoUtils class.

        :param name: The name of the user.
        :param email: The email of the user.
        :param db: The database to search for articles.
        """
        self.name = name
        self.email = email
        self.db = db

    def verify_species(self, species_name: str):
        """
        Verifies if a species name is valid.

        :param species_name: The name of the species to verify.
        :return: None
        """

        # Add your code here to verify the species name

        pass

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
        try:
            # attempts to get a country name
            location = pycountry.countries.get(name=location)
            print(f"{location} is a valid geographic location.")
            # if not, does not return
        except AttributeError:
            print(f"{location} is not a valid geographic location.")
