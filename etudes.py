import pandas as pd
import re
import pycountry
import time
import os
from xml.etree import ElementTree as ET
from Bio import Entrez

#config
Entrez.email = "ayawa-nouwaga@utoulouse.fr"
resultat_max = 80000 
eq_recherche = '(("lung cancer"[Title/Abstract]) OR ("breast cancer"[Title/Abstract]) OR ("colon cancer"[Title/Abstract]) OR ("colorectal cancer"[Title/Abstract]) OR ("prostate cancer"[Title/Abstract])) AND (("income"[Title/Abstract]) OR ("family history"[Title/Abstract]) OR ("stress"[Title/Abstract]) OR ("urban"[Title/Abstract]) OR ("socioeconomic"[Title/Abstract]) OR ("poverty"[Title/Abstract]) OR ("country"[Title/Abstract]) OR ("smoking"[Title/Abstract]) OR ("alcohol"[Title/Abstract]) OR ("education"[Title/Abstract]) OR ("diet"[Title/Abstract])) AND ("2000"[Date - Publication] : "2025"[Date - Publication])'

dict_iso = {c.name.upper(): (c.alpha_3, c.name) for c in pycountry.countries}#cree un dictionnaire avec les noms des pays en majuscules comme clés et (code iso, nom complet) comme valeurs
dict_iso.update({"USA": ("USA", "United States"), "UK": ("GBR", "United Kingdom"), "SOUTH KOREA": ("KOR", "Korea, Republic of")})

def extraire_pays(texte):
    if not texte: return []
    trouves = set()
    texte_clean = texte.upper()
    noms_tries = sorted(dict_iso.keys(), key=len, reverse=True)
    for nom in noms_tries:
        if re.search(r'\b' + re.escape(nom) + r'\b', texte_clean):#expression reguliere pour trouver le nom du pays dans le texte
            trouves.add(dict_iso[nom])
    return list(trouves)

def main():
    print("verif des fichiers et fusion")
    annees = list(range(2000, 2026))
    annees.reverse() 
    fichiers_sauves = []

    for annee in annees:
        nom_fichier = f"temp_annee_{annee}.csv"
        fichiers_sauves.append(nom_fichier)

        #pour eviter de refaire la recherche si le fichier existe deja et a une taille de 500 bytes pour eviter les fichiers vides
        if os.path.exists(nom_fichier) and os.path.getsize(nom_fichier) > 500:
            print(f"annee {annee} deja presente")
            continue

        print(f"extraction pour l'annee {annee}")
        requete_specifique = f"({eq_recherche}) AND ({annee}[Date - Publication])"
        try:
            connexion = Entrez.esearch(db='pubmed', term=requete_specifique, retmax=9999)#demande les pmids pour l'année en cours et 9999 pour eviter de depasser la limite de l'api
            res = Entrez.read(connexion)
            connexion.close()
            ids = res.get('IdList', [])#recup les pmids pour l'année en cours
            if not ids: continue
            donnees_annee = []
            for j in range(0, len(ids), 1000):#decoupe les pmids en paquets de 1000 pour eviter de depasser la limite de l'api
                paquet = ids[j:j+1000]#prend un paquet de pmids
                h_fetch = Entrez.efetch(db='pubmed', id=','.join(paquet), retmode='xml')#demande les details pour les pmids du paquet en cours
                racine = ET.fromstring(h_fetch.read())#lit la reponse de l'api et parse en xml pour pouvoir extraire les infos dont jai besoin
                h_fetch.close()
                for article in racine.findall('.//PubmedArticle'):
                    pmid = article.findtext('.//PMID')
                    abs_text = " ".join([n.text for n in article.findall('.//AbstractText') if n.text])#recup le texte de l abstract pour essayer d y trouver le pays 
                    pays = extraire_pays(abs_text)
                    if pays:
                        for iso, nom in pays:
                            donnees_annee.append({"pmid": pmid, "code_iso": iso, "nom_pays": nom})
                    else:
                        donnees_annee.append({"pmid": pmid, "code_iso": "N/A", "nom_pays": "Unknown"})
                time.sleep(1)
            pd.DataFrame(donnees_annee).to_csv(nom_fichier, index=False)
        except Exception as e:
            print(f" erreur sur l'annee {annee}: {e}")

    print("fusion des fichiers par annee")
    all_data = []
    for f in fichiers_sauves:
        if os.path.exists(f):
            df_annee = pd.read_csv(f, keep_default_na=False)# keep_default_na=False pour garder le N/A en texte et ne pas le transformer en NaN
            all_data.append(df_annee)
    
    if all_data:
        df_final = pd.concat(all_data, ignore_index=True)
        
        #remplace lesvides par N/A
        df_final['code_iso'] = df_final['code_iso'].replace('', 'N/A')
        
        # sauvegarde pays_etude.csv
        df_pays = df_final[['code_iso', 'nom_pays']].drop_duplicates()
        df_pays.to_csv("pays_etude.csv", index=False)
        
        # sauvegarde pays_articles.csv
        df_articles = df_final[['pmid', 'code_iso']].drop_duplicates()
        df_articles.to_csv("pays_articles.csv", index=False)
        
        print("les fichiers pays_etude.csv et pays_articles.csv ont été créés")

if __name__ == '__main__':
    main()