import pandas as pd
import pycountry
import re
from tqdm.auto import tqdm

# Activer tqdm pour pandas
tqdm.pandas()

# Charger la liste des pays (noms officiels en anglais)
countries = pycountry.countries
country_names = [country.name.lower() for country in countries]
country_names.extend([country.official_name.lower() for country in countries if hasattr(country, 'official_name')])

def extract_first_country(text):
    if pd.isna(text):
        return 'NA'
    text_lower = text.lower()
    found_countries = []
    for name in country_names:
        if re.search(r'\b' + re.escape(name) + r'\b', text_lower):
            found_countries.append(name.title())
    return found_countries[0] if found_countries else 'NA'

# Lire ton CSV (adapte 'ton_csv.csv')
df = pd.read_csv('abstracts_pmid.csv')  # Colonnes attendues: 'pmid', 'Abstract'

# Appliquer l'extraction avec une barre de progression
df['pays'] = df['Abstract'].progress_apply(extract_first_country)

# Sauvegarder avec PMID et pays
df[['PMID', 'pays']].to_csv('pmids_pays.csv', index=False)

print(df[['PMID', 'pays']].head())
print("\nExemple de stats pays trouvés:")
print(df['pays'].value_counts().head())
