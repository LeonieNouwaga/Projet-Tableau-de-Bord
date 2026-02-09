import nltk
nltk.download('punkt_tab')
nltk.download('stopwords')
nltk.download('wordnet')
from nltk.tokenize import word_tokenize
from nltk.corpus import stopwords
from nltk.stem import PorterStemmer, WordNetLemmatizer
import pandas as pd
 
df = pd.read_csv('abstracts_pmid.csv')  
data = []
count = 0
for identity in df["PMID"]:
    count += 1
    print((count * 100) / (80000))
    for cell in df['Abstract'].head(count):
        # Tokenisation
        texte = cell
        tokens = word_tokenize(texte.lower())

        # Stopwords 
        stop_words = set(stopwords.words('english'))
        tokens_filtres = [w for w in tokens if w not in stop_words and w.isalpha()]

        # Stemming 
        stemmer = PorterStemmer()
        stems = [stemmer.stem(w) for w in tokens_filtres]


        # Lemmatisation 
        lemmatizer = WordNetLemmatizer()
        lemmas = [lemmatizer.lemmatize(w) for w in tokens_filtres]

        article = {
            'PMID': identity,
            'tokens': tokens,
            'tokens_filtres': tokens_filtres,
            'stems': stems,
            'lemmas': lemmas
            }
        data.append(article)

df_exit = pd.DataFrame(data)
df_exit.to_csv("test.csv", index=False)  
print(df_exit.head())    












