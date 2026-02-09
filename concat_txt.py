import os
import glob

dossier = r"c:\Users\hagaa\OneDrive\projet TB"
os.chdir(dossier)

# Trouve tous les .txt
txt_files = glob.glob("*.txt")
print(f"Trouvé {len(txt_files)} fichiers TXT : {txt_files}")

# Concatène dans un seul fichier
with open("tous_abstracts.txt", "w", encoding="utf-8") as fichier_sortie:
    for fichier in txt_files:
        with open(fichier, "r", encoding="utf-8") as f:
            contenu = f.read()
            fichier_sortie.write(contenu + "\n\n")  # Séparateur entre fichiers

print("✅ tous_abstracts.txt créé !")
