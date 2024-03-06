import pandas as pd
from IPython.display import display, HTML

pd.set_option('display.max_columns', 100)  # Default is 20
pd.set_option('display.max_rows', 100)  # Default is 20

def make_clickable(url):
    if "nan" not in url.lower():
        return f'<a href="{url}" target="_blank">{url}</a>'
    else:
        return ""

df = pd.read_csv("tableS1_big_table_final.csv")

prefix="https://doi.org/"
selcolumns = ['doi (DAT)', 'doi (SERT)', 'doi (NET)']

for icol in selcolumns:
    df[icol] = df[icol].apply(lambda x: prefix +str(x))

for icol in selcolumns:
    df[icol] = df[icol].apply(make_clickable)

df = df.fillna("-")

# Display the DataFrame with clickable links
display(HTML(df.to_html(escape=False)))