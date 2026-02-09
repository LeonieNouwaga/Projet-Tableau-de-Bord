import pandas as pd
import glob
import os

dossier = r"c:\Users\hagaa\OneDrive\projet TB"
os.chdir(dossier)
all_csv = glob.glob("*.csv")


if len(all_csv) == 0:
    print("AUCUN CSV → vérifie tes exports PubMed")
else:
    print("Concaténation...")
    try:
        df_final = pd.concat([pd.read_csv(f) for f in all_csv], ignore_index=True)
        df_final.to_csv("pubmed_complet.csv", index=False)
        print("Sauvegardé : pubmed_complet.csv")
    except Exception as e:
        print(f"Erreur concat : {e}")
print(df_final.head())
