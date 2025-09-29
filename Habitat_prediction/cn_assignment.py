
import os
from os.path import join
import pandas as pd
import numpy as np
from pathlib import Path
import tqdm


data_dir=""
#get all the csv files
csv_files = list(Path(data_dir).glob('*.csv'))
#get all patient ids
patient_ids = [csv_file.stem.split('_')[0] for csv_file in csv_files]
#randomly split the patient ids into training and validation sets
np.random.seed(5)
np.random.shuffle(patient_ids)
train_ids = patient_ids[:int(len(patient_ids)*0.8)]
val_ids = patient_ids[int(len(patient_ids)*0.8):]
train_csvs=[]
for train_id in train_ids:
    train_csv=pd.read_csv(join(data_dir,f'{train_id}_labels.csv'))
    train_csv['image_path'] = train_csv['image_name'].apply(lambda x: join(data_dir,train_id,x))
    train_csvs.append(train_csv)
train_csvs=pd.concat(train_csvs)
val_csvs=[]
for val_id in val_ids:
    val_csv=pd.read_csv(join(data_dir,f'{val_id}_labels.csv'))
    val_csv['image_path'] = val_csv['image_name'].apply(lambda x: join(data_dir,val_id,x))
    val_csvs.append(val_csv)
val_csvs=pd.concat(val_csvs)
label_columns = [f'label_{i}' for i in range(1, 11)]
#get cell_count <=5
train_csvs_others = train_csvs[train_csvs['cell_count'] <= 5].copy()
train_csvs_others['label'] = 0
#get cell_count >15
train_csvs = train_csvs[train_csvs['cell_count'] > 15].copy()
train_csvs = train_csvs[train_csvs[label_columns].max(axis=1) >= 0.6]
#get labels of each row using the maximum value of label columns
train_csvs['label'] = train_csvs[label_columns].idxmax(axis=1).apply(lambda x: int(x.split('_')[1]))
train_csvs['label'] = train_csvs['label'].astype(int)
train_csvs = pd.concat([train_csvs, train_csvs_others])
train_csvs.reset_index(drop=True, inplace=True)

#get validation csvs
val_csvs_others = val_csvs[val_csvs['cell_count'] <= 5].copy()
val_csvs_others['label'] = 0
#get cell_count >15
val_csvs = val_csvs[val_csvs['cell_count'] > 15].copy()
val_csvs = val_csvs[val_csvs[label_columns].max(axis=1) >= 0.6]
val_csvs['label'] = val_csvs[label_columns].idxmax(axis=1).apply(lambda x: int(x.split('_')[1]))
val_csvs['label'] = val_csvs['label'].astype(int)
val_csvs = pd.concat([val_csvs, val_csvs_others])
val_csvs.reset_index(drop=True, inplace=True)

#save train and validation csvs
train_csvs.to_csv(join(data_dir, 'train_labels.csv'), index=False)
val_csvs.to_csv(join(data_dir, 'val_labels.csv'), index=False)