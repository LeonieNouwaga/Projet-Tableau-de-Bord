import pandas as pd

df = pd.read_csv('pubmed_complet.csv')  
df = df.drop('DOI', axis=1) 
df = df.drop('Authors', axis=1)   
df = df.drop('Citation', axis=1)
df = df.drop('Journal/Book', axis=1)
df = df.drop('PMCID', axis=1)
df = df.drop('NIHMS ID', axis=1)
df.to_csv('new_pub.csv', index=False)  
print(df.head())  # Vérifie 

