from epo_utils import EpoUtils
import os
import dotenv

# Load the environment variables
dotenv.load_dotenv()

email = os.getenv("email")

# Create an instance of the Eponyms class
epo = EpoUtils("email")

# Verify a species
epo.verify_species("Harmonia")

# finiid article
epo.return_articles("homo sampiens")
