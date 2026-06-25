# Extra omics tasks

Код для самостоятельных заданий из конца `day5_omics_practice/README.md` и вопросов по WGBS/IGV.

Результаты не хранятся в этой папке. По умолчанию скрипты пишут таблицы и картинки в `extra_credit_results/`.

Запуск:

```bash
conda activate genome
python extra_credit_code/download_light_data.py
python extra_credit_code/run_extra_analysis.py
python extra_credit_code/make_presentation.py
```

Что делает:

- скачивает только легкие готовые RNAseq/CAGE bigWig и GTF-аннотацию;
- сравнивает параметры поиска границ доменов;
- считает связь границ с CTCF;
- извлекает промоторы из GTF;
- считает RNAseq, H3K27Ac, H3K9me3, метилирование и E1 в промоторах;
- сравнивает активные и неактивные промоторы;
- считает расстояние до ближайшего CTCF пика;
- выбирает несколько регионов и строит терминальную замену IGV через `pyGenomeTracks`.
- собирает простую редактируемую презентацию.
