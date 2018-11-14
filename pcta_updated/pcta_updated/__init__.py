import pandas as pd
from .celery import app as celery_app

__all__ = ['celery_app']

all_expr_df = pd.read_csv('user_data/pcta_expression_data.csv',index_col=0)
all_expr_df.index = all_expr_df.index.astype(str)

mra_set = pd.read_csv('user_data/PDI_171116.csv',index_col=0)
mra_set.index = mra_set.index.astype(str)
mra_cad = set(list(all_expr_df.index.tolist())).intersection(set(list(mra_set.index.tolist())))
mra_set = mra_set.loc[mra_cad]

pcta_id = pd.read_csv('user_data/pcta_id_mapping.csv',index_col=0)
pcta_id.index = pcta_id.index.astype(str)

