import re
import pandas as pd


file  = open("tous_abstracts.txt","r",encoding="utf8")
Q_list = file.readlines()
new_file = []
#print(Q_list[1:50])

def count_space(term, count):
    if term == "\n" :
        count += 1
        return count
    else:
        return count

def keep_abstract(new_file, n):
    count = 0
    k = 0
    if n == 1 :
        for terms in Q_list:
            if count != 4 :
                count = count_space(terms, count)
            else:
                count = count_space(terms, count)
                while count != 5:
                    new_file.append(terms)
                    break
        return new_file
    else:
        for terms in Q_list[(n-1):]:
            if count == 4 and word(terms, "Comment in") is True :
                k = 1
            if count == 5 + k :
                break
            if count != 4 + k :
                count = count_space(terms, count)
            else:
                count = count_space(terms, count)
                while count != 5 + k:
                    new_file.append(terms)
                    break
        return new_file

def word(phrase, mot):
    return mot.lower() in phrase.lower()

def delete_words(phrase, mot):
    return phrase.replace(mot, "")

def keep_pmid(new_file, n):
    for terms in Q_list[(n-1):]:
        if word(terms,"pmid") is True :
            phrase = delete_words(terms," [Indexed for MEDLINE]\n")
            new_file.append(phrase)
            return(None)

def each_doc_methods(new_file, n):
    keep_pmid(new_file, n)
    keep_abstract(new_file, n)


def begin_by(phrase):
    if not phrase:  
        return False  
    # Regex : ^début + \d+un ou + chiffres + \.point + \s+un ou + espaces
    pattern = r'^\d+\.\s'
    return bool(re.match(pattern, phrase.strip()))

count = 0
n = 0
for terms in Q_list:
    n += 1
    if begin_by(terms) is True :
        count += 1
        print("loading : ", (count*100)/80068," %")
        each_doc_methods(new_file, n)
    else:
        continue

def parser_liste_melangee(liste_phrases):
    abstracts = []
    current_abstract = []
    current_pmid = None
    
    pmid_pattern = re.compile(r'(?i)pmid[:\-]?\s*(\d+)', re.IGNORECASE)
    
    for phrase in liste_phrases:
        # C'est un PMID ?
        pmid_match = pmid_pattern.search(phrase)
        if pmid_match:
            # Sauvegarde abstract précédent
            if current_abstract and current_pmid:
                abstracts.append({
                    'PMID': current_pmid,
                    'Abstract': ' '.join(current_abstract).strip()
                })
            
            # Nouveau PMID
            current_pmid = pmid_match.group(1)
            current_abstract = []
        else:
            # Phrase normale → ajoute à l'abstract courant
            current_abstract.append(phrase.strip())
    
    # Dernier abstract
    if current_abstract and current_pmid:
        abstracts.append({
            'PMID': current_pmid,
            'Abstract': ' '.join(current_abstract).strip()
        })
    
    return pd.DataFrame(abstracts)

# TON USAGE
ta_liste = new_file

df = parser_liste_melangee(ta_liste)
df.to_csv("abstracts_pmid.csv", index=False)
print(df)
print("✅ CSV créé !")
























