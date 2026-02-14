import pandas as pd
from tqdm import tqdm

# TES LISTES
FACTEURS = ["income", "family history", "stress", "urban", "socioeconomic", 
            "poverty", "country", "smoking", "alcohol", "education", "diet"]

CANCERS = ["lung cancer", "breast cancer", "colon cancer", "colorectal cancer", 
           "prostate cancer", "pancreatic cancer"]

def tout_en_un(csv_title, csv_abstract, csv_final=None):
    """Facteurs + Cancers 0/1 → 1 CSV FINAL propre."""
    df_title = pd.read_csv(csv_title)
    df_abs = pd.read_csv(csv_abstract)
    
    lignes_avant = len(df_title)
    print(f"📊 AVANT: {lignes_avant} lignes")
    
    # Fusion ANTI-Doublons
    df_final = df_title.merge(df_abs[['PMID', 'Abstract']], on='PMID', how='left')
    df_final = df_final.drop_duplicates(subset=['PMID'])
    print(f"✅ Fusion: {len(df_final)} lignes (PMID unique)")
    
    # 1. FACTEURS
    print("\n🔄 === FACTEURS RISQUE ===")
    for sujet in tqdm(FACTEURS, desc="Facteurs"):
        col_out = f"{sujet.replace(' ', '_').title()}_present"
        pattern = r'\b' + sujet.replace(' ', r'\s+') + r'\b'
        df_final[col_out] = df_final['Abstract'].str.contains(
            pattern, case=False, na=False, regex=True
        ).astype(int)
    
    # 2. CANCERS
    print("\n🔄 === TYPES CANCER ===")
    for cancer in tqdm(CANCERS, desc="Cancers"):
        col_out = f"{cancer.replace(' ', '_').replace('-', '_').title()}_present"
        pattern = r'\b' + cancer.replace(' ', r'\s+') + r'\b'
        df_final[col_out] = df_final['Abstract'].str.contains(
            pattern, case=False, na=False, regex=True
        ).astype(int)
    
    # SUPPRIME Abstract
    df_final = df_final.drop(columns=['Abstract'], errors='ignore')
    
    csv_final = csv_final or csv_title.replace(".csv", "_complet.csv")
    df_final.to_csv(csv_final, index=False, encoding='utf-8')
    print(f"\n💾 FINAL: {csv_final} ({len(df_final)} lignes)")
    print("✅ Facteurs + Cancers = 22 cols 0/1 !")
    
    return df_final

# 🎯 Lance !
if __name__ == "__main__":
    csv_precedent = "new_pub.csv"      # PMID|Title|Classe
    csv_abstracts = "abstracts_pmid.csv"  # PMID|Abstract
    
    df_final = tout_en_un(csv_precedent, csv_abstracts)
    
    # Aperçu
    cols_all = ['PMID', 'Title', 'Classe'] + \
               [f"{s.replace(' ', '_').title()}_present" for s in FACTEURS] + \
               [f"{c.replace(' ', '_').replace('-', '_').title()}_present" for c in CANCERS]
    print("\n✅ TOUT (extrait) :")
    print(df_final[cols_all[:10]].head())  # 1er 10 cols
    
    print("\n📈 Stats rapides :")
    stats = df_final[[f"{s.replace(' ', '_').title()}_present" for s in FACTEURS + CANCERS]].sum().sort_values(ascending=False)

