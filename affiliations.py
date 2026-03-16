import os
import time
import pandas as pd
from xml.etree import ElementTree as ET
from Bio import Entrez

#config
Entrez.email = "ayawa-nouwaga@utoulouse.fr"
eq_recherche = '(("lung cancer"[Title/Abstract]) OR ("breast cancer"[Title/Abstract]) OR ("colon cancer"[Title/Abstract]) OR ("colorectal cancer"[Title/Abstract]) OR ("prostate cancer"[Title/Abstract])) AND (("income"[Title/Abstract]) OR ("family history"[Title/Abstract]) OR ("stress"[Title/Abstract]) OR ("urban"[Title/Abstract]) OR ("socioeconomic"[Title/Abstract]) OR ("poverty"[Title/Abstract]) OR ("country"[Title/Abstract]) OR ("smoking"[Title/Abstract]) OR ("alcohol"[Title/Abstract]) OR ("education"[Title/Abstract]) OR ("diet"[Title/Abstract])) AND ("2000"[Date - Publication] : "2025"[Date - Publication])'
resultat_max = 80000
#xml ici car c'est comme ça que PubMed envoie les donnees et c'st plus simple que du texte brut
def extraire_donnees_article(article_xml):
    import re
    auteurs = article_xml.findall('.//Author')
    if not auteurs: 
        return None 
    premier = auteurs[0]
    nom = premier.findtext('LastName')
    prenom = premier.findtext('ForeName')
    initiales = premier.findtext('Initials')
    initiales_final = ""
    if initiales:
        initiales_final = "".join([c for c in initiales if c.isupper()])#ne garde que les majuscules
        
    #si les initiales sont vides ou bizarres, on les fais depuis le prnom
    if not initiales_final and prenom:
        parties = re.split(r'[\s\-]+', prenom)#dcoupe le prenom par espace ou tiret
        initiales_final = "".join([p[0].upper() for p in parties if p])
    nom_final = f"{nom} {initiales_final}" if nom and initiales_final else (nom or "N/A")#met tout ensemble

    affil_info = premier.find('.//AffiliationInfo/Affiliation')
    if affil_info is not None and affil_info.text:
        info = affil_info.text.strip()
        return {
            "nom_auteur": nom_final,
            "institution": info
        }
    return None

def main():
    print("extraction par années")
    #boucle par annee pour garantir de depasser les 9999 articles de la limite pubmed
    annees = list(range(2000, 2026))
    annees.reverse() 
    fichiers_tranches = []

    for annee in annees:
        nom_fichier = f"temp_affil_annee_{annee}.csv"
        fichiers_tranches.append(nom_fichier)

        #si le fichier existe deja, on passe pour ne pas tout re-telecharger
        if os.path.exists(nom_fichier) and os.path.getsize(nom_fichier) > 500:
            print(f"annee {annee} deja presente")
            continue

        print(f"recup annee{annee}")
        requete = f"({eq_recherche}) AND ({annee}[Date - Publication])"
        
        try:
            #cherche les id pour l'annee
            connexion = Entrez.esearch(db='pubmed', term=requete, retmax=9999)#demande les pmids pour l'année en cours et 9999 pour eviter de depasser la limite de l'api
            res = Entrez.read(connexion)
            connexion.close()
            id_annee = res.get('IdList', [])#recup les pmids pour l'année en cours
            print(f"articles trouves pour {annee} : {len(id_annee)}")
            
            if not id_annee:
                continue

            donnees_annee = []
            # telechargement par paquets de 1000 pour la stabilite
            for j in range(0, len(id_annee), 1000):
                paquet_ids = id_annee[j:j+1000]
                connexion = Entrez.efetch(db='pubmed', id=','.join(paquet_ids), retmode='xml')#demande les details pour les pmids du paquet en cours
                racine = ET.fromstring(connexion.read())#lit la reponse de l'api et parse en xml pour pouvoir extraire les infos dont jai besoin
                connexion.close()

                for article in racine.findall('.//PubmedArticle'):#parcourt les articles du paquet
                    info = extraire_donnees_article(article)
                    if info:
                        donnees_annee.append(info)
                time.sleep(1)

            # sauvegarde temporaire de l'annee
            pd.DataFrame(donnees_annee).to_csv(nom_fichier, index=False)

        except Exception as e:
            print(f"   ! Erreur annee {annee} : {e}")
            time.sleep(5)
            continue

    # fusion des fichiers par annee pour creer les fichiers finaux
    print("fusion des donnees et generation des fichiers finaux")
    tous_dfs = []
    for f in fichiers_tranches:
        if os.path.exists(f):
            tous_dfs.append(pd.read_csv(f))
    
    if not tous_dfs:
        print("aucune donnée extraite")
        return

    #regroupe tout et limite 80000
    df_global = pd.concat(tous_dfs, ignore_index=True).head(resultat_max)
    df_global[['institution']].to_csv("all_affiliations.csv", index=False)
    #affiliations.csv
    df_affil = df_global[['institution']].drop_duplicates().reset_index(drop=True)
    df_affil['id_affiliation'] = df_affil.index + 1
    df_affil[['id_affiliation', 'institution']].to_csv("affiliations.csv", index=False)

    # jointure pour recup l'id_affiliation pour chaque auteur
    df_final = pd.merge(df_global, df_affil, on='institution')

    # auteurs_affiliations.csv
    df_final[['nom_auteur', 'id_affiliation']].to_csv("auteurs_affiliations.csv", index=False)

    print("extraction terminee")
    print(f"Articles traites : {len(df_global)}")
    print(f"Affiliations uniques : {len(df_affil)}")
    print(f"all_affiliations.csv ({len(df_global)} lignes)")
    print("Fichiers générés : affiliations.csv, auteurs_affiliations.csv")

if __name__ == '__main__':
    main()