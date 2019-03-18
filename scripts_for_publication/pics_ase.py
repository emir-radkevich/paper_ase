# Импортируем все необходиме модули
# pandas -- для создания датафреймов, как в r, в питоне
# os -- для того, чтобы обращаться к папкам и вообще работать с файловой системой
# seaborn -- для того, чтобы делать красивые картинки
# matplotlib -- этот модуль нужен для seaborn

import numpy as np
import pandas as pd
import os

import seaborn as sns
import matplotlib.pyplot as plt

# Напишем функцию, которая будет загружать наш датасет и добавлять необходимые
# колонки для построения графика

def make_proper_df(df_name):
    # Смотрим, где мы
    where_am_i = os.getcwd()
    # Грузим в df наш датасет
    df = pd.read_csv(os.path.join(where_am_i, df_name), sep='\t')
    # Добавляем в наш df три столбика
    # ase_normal -- процент АСЭ генов в норме
    # ase_tumor -- процент АСЭ генов в опухоли
    # total -- по сути везде пишется 100 % (это нужно для построения графика, 
    # никакой смысловой нагрузки этот столик не несет)
    df['ase_normal'] = df.count_normal/df.total_normal * 100
    df['ase_tumor'] = df.count_tumor/df.total_tumor * 100
    df['total'] = df.total_normal/df.total_normal * 100
    return df


df = make_proper_df('statistics_ext.tsv')

### Сделаем картинку, где будут отображены гены, где частота АСЭ в опухоли выше, чем в норме

# Это то, что будет ИСКЛЮЧЕНО из датасета
# Что исключаем?

# (df.total_normal <= 7) | (df.total_tumor <= 7) -- количество образцов в норме и опухоли,
# которые подавали mbased, больше 7 (число из головы)

# (df.ase_tumor/df.ase_normal <= 2) -- частота АСЭ в опухоли в 2 раза и больше выше,
# чем в норме (просто хотел показать, что таких вариантов очень много)

# (df.ase_normal >= 100) | (df.ase_tumor >= 100) -- это проверка на то, чтобы количество
# АСЭ образцов было НЕ больше, чем тестируемых образцов (такое действительно есть!)

statement = (df.total_normal <= 7) | (df.total_tumor <= 7) |     (df.ase_tumor/df.ase_normal <= 2) | (df.ase_normal >= 100) | (df.ase_tumor >= 100)

df = df.drop(df[statement].index)
df = df.reset_index(drop=True)
pd.set_option('display.max_rows', 10)

# Создаем картинку вида stacked barplot в seaborn

sns.set(style="whitegrid")

f, ax = plt.subplots(figsize=(5, 15))

sns.set_color_codes("muted")
sns.barplot(x='total', y="gene", data=df,
            label="total", color="r")

sns.set_color_codes("muted")
sns.barplot(x='ase_tumor', y="gene", data=df,
            label="ase_tumor", color='pink')

sns.set_color_codes("pastel")
sns.barplot(x="ase_normal", y="gene", data = df,
            label="ase_normal", color="b")

ax.legend(ncol=1, loc="upper right", frameon=True)
ax.set(ylabel="gene",
       xlabel="fraction, %")
sns.despine(left=True, bottom=True)

plt.savefig('high_tumor_ase_frequency.png', dpi=400, bbox_inches='tight')


### Теперь сделаем картинку, где отобразим все гены, у которых хотя бы в половине образцов наблюдается АСЭ (и в норме, и в опухоли) + частота АСЭ в опухоли выше, чем в норме

df = make_proper_df('statistics_ext.tsv')

# Это то, что будет ИСКЛЮЧЕНО из датасета
# Что исключаем?

# Большинство условий такие же, как в предыдущем условии
# (df.ase_tumor/df.ase_normal <= 1.0) -- частота АСЭ в опухоли выше, чем в норме

# (df.count_normal/df.total_normal <= 0.5) | (df.count_tumor/df.total_tumor <= 0.5) -- 
# хотя бы в половине образцов наблюдается АСЭ (и в норме, и в опухоли)

statement2 = (df.total_normal <= 7) | (df.total_tumor <= 7) | (df.count_tumor/df.total_tumor <= 0.5) |     (df.count_normal/df.total_normal <= 0.5) | (df.ase_normal >= 100) | (df.ase_tumor >= 100) |     (df.ase_tumor/df.ase_normal <= 1.0)
     
df2 = df.drop(df[statement2].index)
df2 = df2.reset_index(drop=True)
pd.set_option('display.max_rows', 10)

# Делаем картинку

sns.set(style="whitegrid")

f, ax = plt.subplots(figsize=(5, 13))

sns.set_color_codes("muted")
sns.barplot(x='total', y="gene", data=df2,
            label="total", color="r")

sns.set_color_codes("muted")
sns.barplot(x='ase_tumor', y="gene", data=df2,
            label="ase_tumor", color='pink')

sns.set_color_codes("pastel")
sns.barplot(x="ase_normal", y="gene", data = df2,
            label="ase_normal", color="b")

ax.legend(ncol=1, loc="upper right", frameon=True)
ax.set(ylabel="gene",
       xlabel="fraction, %")
sns.despine(left=True, bottom=True)

plt.savefig('high_ase_frequency.png', dpi=400, bbox_inches='tight')

