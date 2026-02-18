import os
import time
import pandas as pd
from xml.etree import ElementTree as ET #pour lire le format XML
from Bio import Entrez #librairie de python pour se connecter à l'api de pubmed

eq_recherche = '(("lung cancer"[Title/Abstract]) OR ("breast cancer"[Title/Abstract]) OR ("colon cancer"[Title/Abstract]) OR ("colorectal cancer"[Title/Abstract]) OR ("prostate cancer"[Title/Abstract])) AND (("income"[Title/Abstract]) OR ("family history"[Title/Abstract]) OR ("stress"[Title/Abstract]) OR ("urban"[Title/Abstract]) OR ("socioeconomic"[Title/Abstract]) OR ("poverty"[Title/Abstract]) OR ("country"[Title/Abstract]) OR ("smoking"[Title/Abstract]) OR ("alcohol"[Title/Abstract]) OR ("education"[Title/Abstract]) OR ("diet"[Title/Abstract])) AND ("2000"[Date - Publication] : "2025"[Date - Publication])'

max_res = 10000 #nbre max d'articles qu'on recup pr l instant

#fichiers csv
file_mots = "keyword_unique.csv"
file_liens = "pmid_keyword.csv"


# recup les pmids
def recuperer_pmids(requete, max_res=500):
    Entrez.email = "ayawa-nouwaga@utoulouse.fr" #pour s'identifier dans l'api pubmed

    if not Entrez.email:
        print("Entrez une adresse mail")

    try:
        print(f"Lancement de la recherche Entrez pour la requête (max {max_res})...")
        # esearch pour la recherche sur la base de données pubmed
        # term=requete utilise l'equation recherche.
        # retmax=max_res limite le nombre de resats
        connexion = Entrez.esearch(db='pubmed', term=requete, retmax=max_res)
        
        #lit la reponse envoyee par pubmed
        res = Entrez.read(connexion)
        
        #ferme la connexion
        connexion.close()
        
        # extrait les pmids trouvés sinon liste vide
        pmids = res.get('IdList', [])
        print(f"Recherche terminée — {len(pmids)} PMIDs trouvés.")
        return pmids
        
    except Exception as e:
        print(f"Erreur lors de la recherche Entrez: {e}")
        return []

# pour couper une liste de 10000 pmids en petits paquets de 200 car l'api ne peut pas tout envoyer d'un seul coup
def decouper_liste(liste, n):
    for i in range(0, len(liste), n):
        yield liste[i:i+n]


# extrait les mots cles d'un article
def extraire_mots_cles(article_xml):
    mots = []
    # On récupère les mots bruts
    kws = article_xml.findall('.//Keyword')
    for k in kws:
        if k is not None and k.text:
            mots.append(k.text.strip())

    donnees_nettoyees = []
    
    symboles_a_supprimer = "()[].,;:\'\""

    for terme in mots:
        #remplace chaque symbole par du vide
        t = terme
        for symbole in symboles_a_supprimer:
            t = t.replace(symbole, "")
        
        #refait un strip pour les espaces restants
        t = t.strip() #pour nettoyer

        if not t: continue
        # si le mot contient un chiffre on le jette
        if any(char.isdigit() for char in t):
            continue

        #le mot doit faire au moins 3 lettres ça élimine les mots cles bizarres
        if len(t) < 3:
            continue

        #le mot doit commencer par une lettre
        if not t[0].isalpha():
            continue

        # met en majuscule
        t_propre = t.capitalize()

        # on garde si ce n'est pas un doublon
        if t_propre not in donnees_nettoyees:
            donnees_nettoyees.append(t_propre)

    return donnees_nettoyees


def main():
    pmids = recuperer_pmids(eq_recherche, max_res=max_res)
    
    if not pmids:
        print("Aucun pmid trouvé pour la requête")
        return
    print(f"Début du traitement de {len(pmids)} PMIDs...")

    resultats = []

    #traite les pmid par paquets de 200
    for paquet in decouper_liste(pmids, 200):
        #transforme la liste ['123', '456'] en une chaine "123,456" pour l'envoyer à l'API.
        ids = ','.join(paquet)
        print(f"Traitement d'un paquet de {len(paquet)} PMIDs...")
        
        try:
            # efetch demande les détails complets  pour ces pmids
            connexion = Entrez.efetch(db='pubmed', id=ids, retmode='xml')
            data = connexion.read() # On lit les données reçues.
            connexion.close()
            print("efetch réussi pour le paquet (données reçues)")
            
        except Exception as e:
            print(f"Erreur efetch pour {ids[:50]}...: {e}")
            time.sleep(1)
            continue

        try:
            #convertit le texte reçu en xml 
            racine = ET.fromstring(data)
            nb_articles = len(racine.findall('.//PubmedArticle'))
            print(f"Parsing XML OK — {nb_articles} articles dans le paquet")
        except Exception as e:
            print(f"Erreur parsing xml: {e}")
            continue

        #parcourt chaque article reçu dans le paquet xml
        for article in racine.findall('.//PubmedArticle'):
            #cherche le pmid de l'article courant
            pmid_node = article.find('.//PMID')
            # si on le trouve, on prend le texte sinon none
            pmid = pmid_node.text if pmid_node is not None else None
            
            if not pmid:
                continue

            mots = extraire_mots_cles(article)
            if mots:
                print(f"PMID {pmid}: {len(mots)} mot(s)-clé(s) extraits")
            
            #ajoute chaque couple (pmid, mot cle) à la liste
            for mot in mots:
                resultats.append({"PMID": pmid, "Mot_clés": mot})

        
        # faire une pause pour respecter les règles d'utilisation de pubmed pour ne pas surcharger le serveur
        time.sleep(0.4)

    if resultats:
        # la liste de résultats est transformee en dataframe
        df = pd.DataFrame(resultats)
        
        # sauvegarde le fichier pmid_keyword.csv
        df.to_csv(file_liens, index=False, encoding='utf-8')# index=False evite d'ajouter une colonne de numérotation 0,1,2 etc
        print(f"fichier {file_liens} crée")

        # crée le fichier keyword.csv
        df_mots_uniques = df[['Mot_clés']].drop_duplicates().sort_values(by='Mot_clés')
        
        # sauvegarde du fichier keyword.csv
        df_mots_uniques.to_csv(file_mots, index=False, encoding='utf-8')
        print(f"fichier {file_mots} crée")
       
    else:
        print("aucune donnee recuperee.")

if __name__ == '__main__':
    main()