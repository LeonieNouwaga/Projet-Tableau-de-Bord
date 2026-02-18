import pandas as pd

# Lire ton CSV
df = pd.read_csv('pmid_tokens_long.csv')  # 'PMID', 'tokens'

pmid_col = 'PMID'  # Adapte si besoin
tokens_col = 'lemmas'

# Split sur espaces
df[tokens_col] = df[tokens_col].astype(str).str.split()

# Explode
df_long = df.explode(tokens_col).reset_index(drop=True)

# Nettoyage : strip + remove [, ], ,
df_long[tokens_col] = df_long[tokens_col].astype(str).str.strip()
df_long[tokens_col] = df_long[tokens_col].str.replace(r'[\[\],]', '', regex=True)

# Drop vides/NA
df_long = df_long[df_long[tokens_col] != ''].dropna(subset=[tokens_col]).reset_index(drop=True)

# Sauvegarde
df_long[[pmid_col, tokens_col]].to_csv('pmid_tokens_long.csv', index=False)

print("Aperçu nettoyé :")
print(df_long[[pmid_col, tokens_col]].head(10))
print(f"Lignes : {len(df_long)}")

