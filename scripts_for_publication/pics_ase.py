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
    
    # Добавляем в наш df четыре столбика
    # ase_normal -- процент АСЭ генов в норме
    # ase_tumor -- процент АСЭ генов в опухоли
    # total_percent -- по сути везде пишется 100 % (это нужно для построения графика, 
    # никакой смысловой нагрузки этот столик не несет)
    # total_count -- количество образцов (и нормы, и опухоли), которые тестировали
    df['ase_normal'] = df.count_normal/df.total_normal * 100
    df['ase_tumor'] = df.count_tumor/df.total_tumor * 100
    df['total_percent'] = df.total_normal/df.total_normal * 100
    df['total_count'] = df.total_normal + df.total_tumor
    
    return df



pd.set_option('display.max_rows', 10)

# Грузим список наших генов
df = make_proper_df('statistics_ext.tsv')

# Грузим сегментрано дуплицированные гены
where_am_i = os.getcwd()
df_segm_dupl = pd.read_csv(os.path.join(where_am_i, 'segm_dupl.tsv'), sep='\t')

# Тут перечислены все гены, которые не являются супер дуплицированными
df_nsd = df.query("gene not in @df_segm_dupl.gene")
df_nsd = df_nsd.reset_index(drop=True)

# Тут перечислены все гены, для которых показан геномный импринтинг
# взял с geneimprint.com
df_imprinted = pd.read_csv(os.path.join(where_am_i, 'imprinted_genes.tsv'))

# Потенциально импринтированные гены
df_nsd_imp = df_nsd.query("gene in @df_imprinted.Gene")
df_nsd_imp = df_nsd_imp.reset_index(drop=True)

### Сделаем картинку, где будут отображены гены, где частота АСЭ в опухоли выше, чем в норме

# Это то, что будет оставлено в датасете

# (total_normal > 7) & (total_tumor > 7) -- количество образцов в норме и опухоли,
# которые подавали mbased, больше 7 (число из головы)

# (ase_tumor/ase_normal <= 2) -- частота АСЭ в опухоли в 2 раза и больше выше,
# чем в норме (просто хотел показать, что таких вариантов очень много)

# (ase_normal < 100) & (ase_tumor < 100) -- это проверка на то, чтобы количество
# АСЭ образцов было НЕ больше, чем тестируемых образцов (такое действительно есть!)

df_htase = df_nsd.query('(total_normal > 7) & (total_tumor > 7) & (ase_tumor/ase_normal > 2) & (ase_normal < 100) & (ase_tumor < 100)')
df_htase = df_htase.reset_index(drop=True)

# Создаем картинку вида stacked barplot в seaborn

sns.set(style="whitegrid")

f, ax = plt.subplots(figsize=(5, 15))

sns.set_color_codes("muted")
sns.barplot(x='total_percent', y="gene", data=df_htase,
            label="total", color="r")

sns.set_color_codes("muted")
sns.barplot(x='ase_tumor', y="gene", data=df_htase,
            label="ase_tumor", color='pink')

sns.set_color_codes("pastel")
sns.barplot(x="ase_normal", y="gene", data = df_htase,
            label="ase_normal", color="b")

ax.legend(ncol=1, loc="upper right", frameon=True)
ax.set(ylabel="gene",
       xlabel="fraction, %")
sns.despine(left=True, bottom=True)

plt.savefig('high_tumor_ase_frequency.png', dpi=400, bbox_inches='tight')


### Гены, где частота АСЭ в норме выше, чем в опухоли

# Все гены, у которых (ase_normal/ase_tumor > 2)
df_hnase = df_nsd.query('(total_normal > 7) & (total_tumor > 7) & (ase_normal/ase_tumor > 2) & (ase_normal < 100) & (ase_tumor < 100)')
df_hnase = df_hnase.reset_index(drop=True)

# Создаем картинку вида stacked barplot в seaborn

sns.set(style="whitegrid")

f, ax = plt.subplots(figsize=(3, 5))

sns.set_color_codes("muted")
sns.barplot(x='total_percent', y="gene", data=df_hnase,
            label="total", color="r")

sns.set_color_codes("pastel")
sns.barplot(x="ase_normal", y="gene", data = df_hnase,
            label="ase_normal", color="b")

sns.set_color_codes("muted")
sns.barplot(x='ase_tumor', y="gene", data=df_hnase,
            label="ase_tumor", color='pink')

ax.legend(ncol=1, loc="upper right", frameon=True)
ax.set(ylabel="gene",
       xlabel="fraction, %")
sns.despine(left=True, bottom=True)

plt.savefig('high_normal_ase_frequency.png', dpi=400, bbox_inches='tight')


### Теперь сделаем картинку, где отобразим все гены, у которых хотя бы в половине образцов наблюдается АСЭ (и в норме, и в опухоли) + частота АСЭ в опухоли выше, чем в норме

# Это то, что будет оставлено в датасете

# Большинство условий такие же, как в предыдущем условии
# (ase_tumor/ase_normal > 1.0) -- частота АСЭ в опухоли выше, чем в норме

# (count_normal/total_normal > 0.5) & (count_tumor/total_tumor > 0.5) -- 
# хотя бы в половине образцов наблюдается АСЭ (и в норме, и в опухоли)
 
df_hase = df_nsd.query('(total_normal > 7) & (total_tumor > 7) & (count_tumor/total_tumor > 0.5) &     (count_normal/total_normal > 0.5) & (ase_normal < 100) & (ase_tumor < 100) &     (ase_tumor/ase_normal > 1.0)')
df_hase = df_hase.reset_index(drop=True)

# Делаем картинку

sns.set(style="whitegrid")

f, ax = plt.subplots(figsize=(3, 7))

sns.set_color_codes("muted")
sns.barplot(x='total_percent', y="gene", data=df_hase,
            label="total", color="r")

sns.set_color_codes("muted")
sns.barplot(x='ase_tumor', y="gene", data=df_hase,
            label="ase_tumor", color='pink')

sns.set_color_codes("pastel")
sns.barplot(x="ase_normal", y="gene", data = df_hase,
            label="ase_normal", color="b")

ax.legend(ncol=1, loc="upper right", frameon=True)
ax.set(ylabel="gene",
       xlabel="fraction, %")
sns.despine(left=True, bottom=True)

plt.savefig('high_ase_frequency.png', dpi=400, bbox_inches='tight')

