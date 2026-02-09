import pandas as pd
df = pd.read_csv('test_powerbi.csv')  
df = df.drop('tokens_str', axis=1)
df.to_csv('tokens.csv', index=False)  
print(df.head())  # Vérifie 