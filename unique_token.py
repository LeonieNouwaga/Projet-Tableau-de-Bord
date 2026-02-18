import pandas as pd

# Lire ton CSV (adapte 'ton_csv.csv')
df = pd.read_csv('tokens_final.csv')  # Colonnes : 'PMID' (ou 'pmid'), 'tokens' (ou nom similaire)

# Identifier les colonnes (adapte si différent)
pmid_col = 'PMID'  # Ou 'pmid'
tokens_col = 'lemmas'  # Ou nom de ta colonne tokens

# Split la chaîne tokens sur espaces → liste
df[tokens_col] = df[tokens_col].astype(str).str.split()

# Explode : 1 ligne par token, PMID répété
df_long = df.explode(tokens_col).reset_index(drop=True)

# Nettoyer : strip, drop vides/NA
df_long[tokens_col] = df_long[tokens_col].str.strip()
df_long = df_long[df_long[tokens_col] != ''].dropna(subset=[tokens_col]).reset_index(drop=True)

# Sauvegarder
df_long[[pmid_col, tokens_col]].to_csv('pmid_tokens_long.csv', index=False)

print("Aperçu :")
print(df_long[[pmid_col, tokens_col]].head(10))
print(f"\nLignes total : {len(df_long)}")
