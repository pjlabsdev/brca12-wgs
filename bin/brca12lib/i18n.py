"""Report language. Result JSON values stay German (stable schema); HTML/Markdown are rendered in `de` or `en`."""
import re

UI = {
    "title": ("BRCA1/BRCA2-Bericht", "BRCA1/BRCA2 report"),
    "important": ("Wichtig", "Important"),
    "reference": ("Referenz", "Reference"),
    "input": ("Eingabe", "Input"),
    "snv_caller": ("SNV/Indel-Caller", "SNV/indel caller"),
    "created": ("Erstellt", "Created"),
    "strand_minus": ("Minus-Strang", "minus strand"),
    "strand_plus": ("Plus-Strang", "plus strand"),
    "small_title": ("{gene} SNVs/Indels ({chrom}, {strand}, {transcript})",
                    "{gene} SNVs/indels ({chrom}, {strand}, {transcript})"),
    "small_headers": (["POS", "REF", "ALT", "Genotyp", "Konsequenz", "Transkript", "HGVSc", "HGVSp",
                       "ClinVar", "ACMG/AMP-Hinweise", "Qualität"],
                      ["POS", "REF", "ALT", "Genotype", "Consequence", "Transcript", "HGVSc", "HGVSp",
                       "ClinVar", "ACMG/AMP hints", "Quality"]),
    "no_small": ("Keine Varianten gegenüber der Referenz gefunden.", "No variants relative to the reference."),
    "cnv_title": ("BRCA1/BRCA2 CNVs (Deletionen/Duplikationen)", "BRCA1/BRCA2 CNVs (deletions/duplications)"),
    "cnv_headers": (["Gen", "Typ", "Bereich", "CN", "log2/Ratio", "Exons", "Hinweise"],
                    ["Gene", "Type", "Region", "CN", "log2/ratio", "Exons", "Notes"]),
    "sv_title": ("BRCA1/BRCA2 strukturelle Varianten", "BRCA1/BRCA2 structural variants"),
    "sv_headers": (["Gen", "Typ", "Position", "Partner", "GT", "PR", "SR", "FILTER", "Hinweise"],
                   ["Gene", "Type", "Position", "Mate", "GT", "PR", "SR", "FILTER", "Notes"]),
    "no_sv": ("Manta hat keine SV mit Bruchpunkt in BRCA1/BRCA2 gemeldet.",
              "Manta reported no SV with a breakpoint in BRCA1/BRCA2."),
    "exon_title": ("Coverage je Exon (MAPQ >= 20)", "Coverage per exon (MAPQ >= 20)"),
    "exon_headers": (["Gen", "Exon", "Bereich", "Mittel", "% >=10x", "% >=20x", "% >=30x", "Anteil MAPQ<20",
                      "Status", "Grund"],
                     ["Gene", "Exon", "Region", "Mean", "% >=10x", "% >=20x", "% >=30x", "Fraction MAPQ<20",
                      "Status", "Reason"]),
    "low_title": ("Unsichere Bereiche (Tiefe < {depth:g}x)", "Uncertain regions (depth < {depth:g}x)"),
    "low_headers": (["Gen", "Bereich", "Länge (bp)", "min. Tiefe"], ["Gene", "Region", "Length (bp)", "Min. depth"]),
    "no_low": ("Keine Bereiche unter der Schwelle.", "No regions below the threshold."),
    "no_entries": ("Keine Einträge.", "No entries."),
    "versions": ("Tool-Versionen", "Tool versions"),
    "version": ("Version", "Version"),
    "unknown": ("unbekannt", "unknown"),
    "commands": ("Reproduzierbare Befehle", "Reproducible commands"),
    "commands_note": ("Alle ausgeführten Befehle stehen vollständig in", "Every executed command is listed verbatim in"),
    "limitations": ("Grenzen", "Limitations"),
    "na": ("n. v.", "n/a"),
    "mosaic_title": ("Mosaik-Kandidaten (Allelanteil {lo:g} bis {hi:g} %, GATK Mutect2)",
                     "Mosaic candidates (allele fraction {lo:g} to {hi:g}%, GATK Mutect2)"),
    "mosaic_headers": (["Gen", "POS", "REF", "ALT", "Allelanteil", "Alt-Reads", "Tiefe", "Konsequenz", "HGVSc",
                        "HGVSp", "ClinVar", "Hinweise"],
                       ["Gene", "POS", "REF", "ALT", "Allele fraction", "Alt reads", "Depth", "Consequence", "HGVSc",
                        "HGVSp", "ClinVar", "Notes"]),
    "no_mosaic": ("Keine Variante mit niedrigem Allelanteil gefunden.", "No variant with a low allele fraction found."),
    "qc_title": ("Lese- und Ziel-Qualität (Assay: {assay})", "Read and target quality (assay: {assay})"),
    "qc_headers": (["Kennzahl", "Wert"], ["Metric", "Value"]),
    "qc_reads": ("Reads in den Genbereichen", "Reads in the gene regions"),
    "qc_mapped": ("kartiert (%)", "mapped (%)"),
    "qc_paired": ("korrekt gepaart (% der kartierten)", "properly paired (% of mapped)"),
    "qc_dup": ("Duplikate (% der kartierten)", "duplicates (% of mapped)"),
    "qc_mapq0": ("MAPQ 0 (% der kartierten)", "MAPQ 0 (% of mapped)"),
    "qc_insert": ("Insertgröße, Median (bp)", "insert size, median (bp)"),
    "qc_error": ("Fehlerrate (%)", "error rate (%)"),
    "qc_readlen": ("mittlere Leselänge (bp)", "mean read length (bp)"),
    "qc_basequal": ("mittlere Basenqualität", "mean base quality"),
    "qc_targets": ("Ziele in BRCA1/BRCA2", "targets in BRCA1/BRCA2"),
    "qc_target_mean": ("mittlere Zieltiefe (x)", "mean target depth (x)"),
    "qc_target_pct20": ("Zielbasen >= 20x (%)", "target bases >= 20x (%)"),
    "qc_uniformity": ("Gleichmäßigkeit: Zielbasen >= 0,2 x Mittel (%)", "uniformity: target bases >= 0.2 x mean (%)"),
    "qc_fold80": ("Fold-80-Basenstrafe", "fold-80 base penalty"),
    "qc_ontarget": ("Reads auf Zielen (% aller kartierten)", "on-target reads (% of all mapped)"),
    "qc_weak": ("schwach abgedeckte Ziele", "weakly covered targets"),
    "qc_clipped": ("Reads mit abgeschnittenem Primer (%)", "reads with a clipped primer (%)"),
    "qc_short": ("nach Primer-Abschnitt zu kurze Reads", "reads too short after primer clipping"),
    "qc_flag": ("Warnung", "Warning"),
    "weak_title": ("Schwach abgedeckte Ziele (mittlere Tiefe < {depth:g}x)", "Weakly covered targets (mean depth < {depth:g}x)"),
    "weak_headers": (["Ziel", "Bereich", "Mittel", "Minimum", "% >=20x"], ["Target", "Region", "Mean", "Minimum", "% >=20x"]),
    "no_weak": ("Alle Ziele ausreichend abgedeckt.", "All targets sufficiently covered."),
}

EXACT = {
    "nicht bestimmt": "not determined", "homozygot Referenz": "homozygous reference",
    "homozygot": "homozygous", "heterozygot": "heterozygous",
    "compound/multiallelisch": "compound/multiallelic",
    "unsicher": "uncertain", "ok": "ok",
    "Mosaik-Kandidat": "mosaic candidate", "schwach abgedeckt": "weakly covered",
    "nicht im Zielbereich des Panels": "not covered by a panel target",
    "Ion Torrent: Insertionen/Deletionen in Homopolymeren unzuverlässig":
        "Ion Torrent: insertions/deletions in homopolymers unreliable",
    "niedriger Allelanteil: Bestätigung mit orthogonaler Methode nötig":
        "low allele fraction: confirmation with an orthogonal method required",
    "Deletion": "Deletion", "Duplikation": "Duplication",
    "möglicher Verlust": "possible loss", "möglicher Zugewinn": "possible gain",
    "geringe Tiefe (DP < 15)": "low depth (DP < 15)",
    "geringe Genotypqualität (GQ < 20)": "low genotype quality (GQ < 20)",
    "Allelverhältnis auffällig (mögliches Mosaik oder Artefakt)": "unusual allele balance (possible mosaic or artefact)",
    "geringe Mapping-Qualität (MQ < 40)": "low mapping quality (MQ < 40)",
    "unpräzise Bruchpunkte": "imprecise breakpoints",
    "Bestätigung mit orthogonaler Methode nötig": "confirmation with an orthogonal method required",
    "Kandidat, Bestätigung per MLPA oder validierter CNV-Analyse nötig":
        "candidate, confirmation by MLPA or validated CNV analysis required",
    "eingefügte Sequenz mit Alu-Signatur: mögliche Alu-Insertion":
        "inserted sequence with an Alu signature: possible Alu insertion",
    "abgeschnittene Leseenden mit Alu-Sequenz: mögliche Alu-Insertion":
        "clipped read ends with Alu sequence: possible Alu insertion",
    "abgeschnittene Leseenden von beiden Seiten: mögliche Insertion":
        "clipped read ends from both sides: possible insertion",
    "Hinweis aus Leseenden, kein Variantenaufruf; Bestätigung mit orthogonaler Methode nötig":
        "hint from read ends, not a variant call; confirmation with an orthogonal method required",
    "Insertionssignatur": "insertion signature",
    "Tandemkopie-Signatur": "tandem copy signature",
    "zusätzliche 126-bp-Kopie vor BRCA1 Exon 1 (häufiges Allel)": "extra 126 bp copy upstream of BRCA1 exon 1 (common allele)",
    "bekannte zusätzliche Tandemkopie": "known extra tandem copy",
    "Markerbase der zusätzlichen Kopie in den Lesungen vorhanden": "marker base of the extra copy present in the reads",
    "gekoppelter Haplotyp der zusätzlichen Kopie vorhanden": "linked haplotype of the extra copy present",
    "bekannte intronische Alu-Insertion in BRCA2 (häufiges Allel)": "known intronic Alu insertion in BRCA2 (common allele)",
    "bekanntes Insertionsallel": "known insertion allele",
    "keine Lesung für das Referenzallel: vermutlich homozygot": "no read for the reference allele: probably homozygous",
    "Lesungen für beide Allele: vermutlich heterozygot": "reads for both alleles: probably heterozygous",
    "Lesungen passen besser zur bekannten Allelsequenz als zur Referenz (Neuausrichtung)":
        "reads fit the known allele sequence better than the reference (realignment)",
    "Haplotyp homozygot: vermutlich zwei Allele mit der zusätzlichen Kopie":
        "haplotype homozygous: probably two alleles with the extra copy",
    "Haplotyp heterozygot: vermutlich ein Allel mit der zusätzlichen Kopie":
        "haplotype heterozygous: probably one allele with the extra copy",
    "abgeschnittene Leseenden an beiden Grenzen einer Wiederholungseinheit: mögliche zusätzliche Tandemkopie":
        "clipped read ends at both borders of a repeat unit: possible extra tandem copy",
    "Bruchpunktpaar an derselben Stelle: mögliche Insertion (z. B. mobiles Element wie Alu)":
        "pair of breakends at the same site: possible insertion (e.g. a mobile element such as Alu)",
    "Exonweise Tiefe relativ zum Chromosomen-Median; nur ein Hinweis, kein CNV-Call. Exon ist zugleich als unsicher markiert: Verlust oder Abdeckungsproblem, nicht unterscheidbar.":
        "Exon depth relative to the chromosome median; a hint only, not a CNV call. The exon is also flagged as uncertain: loss or coverage problem, cannot be told apart.",
    "Exonweise Tiefe relativ zum Chromosomen-Median; nur ein Hinweis, kein CNV-Call.":
        "Exon depth relative to the chromosome median; a hint only, not a CNV call.",
}

PATTERNS = [
    (re.compile(r"^hoher Duplikatanteil \((.+) %\)$"), r"high duplicate fraction (\1%)"),
    (re.compile(r"^geringer Anteil kartierter Reads \((.+) %\)$"), r"low fraction of mapped reads (\1%)"),
    (re.compile(r"^hohe Fehlerrate \((.+) %\)$"), r"high error rate (\1%)"),
    (re.compile(r"^(\d+) Ziel\(e\) mit mittlerer Tiefe < (.+)x$"), r"\1 target(s) with mean depth < \2x"),
    (re.compile(r"^ungleichmäßige Abdeckung \(nur (.+) % der Zielbasen >= 0,2 x Mittel\)$"),
     r"uneven coverage (only \1% of target bases >= 0.2 x mean)"),
    (re.compile(r"^nur (.+) % der Reads mit abgeschnittenem Primer$"), r"only \1% of reads with a clipped primer"),
    (re.compile(r"^mittlere Tiefe (.+)x < (.+)x$"), r"mean depth \1x < \2x"),
    (re.compile(r"^nur (.+) % der Basen >= 20x$"), r"only \1% of bases >= 20x"),
    (re.compile(r"^(.+) % der Reads mit MAPQ < 20 \(Pseudogen/Repeat\)$"), r"\1% of reads with MAPQ < 20 (pseudogene/repeat)"),
    (re.compile(r"^liegt in unsicher abgedecktem Bereich (.+)$"), r"located in uncertain region \1"),
    (re.compile(r"^Exon (\d+) \(Legacy (\d+)\)$"), r"Exon \1 (legacy \2)"),
    (re.compile(r"^ENIGMA Tabelle 4 \(falls bestätigt\): (.+)$"), r"ENIGMA Table 4 (if confirmed): \1"),
]


def ui(key, lang, **kwargs):
    value = UI[key][0 if lang == "de" else 1]
    return value.format(**kwargs) if kwargs and isinstance(value, str) else value


def tr(text, lang):
    """Translate a German result value for display; unknown values are returned unchanged."""
    if lang == "de" or not isinstance(text, str):
        return text
    if text in EXACT:
        return EXACT[text]
    for pattern, replacement in PATTERNS:
        if pattern.match(text):
            return pattern.sub(replacement, text)
    return text


# --- Russian (ru) ----------------------------------------------------------------------------------
# Machine-assisted translation, not reviewed by a native speaker (see README). Values without a
# Russian entry fall back to English.
UI_RU = {
    "title": "Отчёт BRCA1/BRCA2",
    "important": "Важно",
    "reference": "Референс",
    "input": "Входные данные",
    "snv_caller": "Вызов SNV/инделов",
    "created": "Создан",
    "strand_minus": "минус-цепь",
    "strand_plus": "плюс-цепь",
    "small_title": "{gene}: SNV/инделы ({chrom}, {strand}, {transcript})",
    "small_headers": ["POS", "REF", "ALT", "Генотип", "Последствие", "Транскрипт", "HGVSc", "HGVSp",
                      "ClinVar", "Указания ACMG/AMP", "Качество"],
    "no_small": "Вариантов относительно референса не найдено.",
    "cnv_title": "BRCA1/BRCA2: CNV (делеции/дупликации)",
    "cnv_headers": ["Ген", "Тип", "Область", "CN", "log2/отношение", "Экзоны", "Примечания"],
    "sv_title": "BRCA1/BRCA2: структурные варианты",
    "sv_headers": ["Ген", "Тип", "Позиция", "Партнёр", "GT", "PR", "SR", "FILTER", "Примечания"],
    "no_sv": "Manta не сообщила о SV с точкой разрыва в BRCA1/BRCA2.",
    "exon_title": "Покрытие по экзонам (MAPQ >= 20)",
    "exon_headers": ["Ген", "Экзон", "Область", "Среднее", "% >=10x", "% >=20x", "% >=30x", "Доля MAPQ<20",
                     "Статус", "Причина"],
    "low_title": "Неопределённые области (глубина < {depth:g}x)",
    "low_headers": ["Ген", "Область", "Длина (п. н.)", "Мин. глубина"],
    "no_low": "Областей ниже порога нет.",
    "mosaic_title": "Кандидаты в мозаицизм (доля аллеля {lo:g}–{hi:g} %, GATK Mutect2)",
    "mosaic_headers": ["Ген", "POS", "REF", "ALT", "Доля аллеля", "Риды ALT", "Глубина", "Последствие", "HGVSc",
                       "HGVSp", "ClinVar", "Примечания"],
    "no_mosaic": "Вариантов с низкой долей аллеля не найдено.",
    "qc_title": "Качество ридов и целевых областей (тип анализа: {assay})",
    "qc_headers": ["Показатель", "Значение"],
    "qc_reads": "Риды в областях генов",
    "qc_mapped": "картировано (%)",
    "qc_paired": "правильно спаренные (% картированных)",
    "qc_dup": "дубликаты (% картированных)",
    "qc_mapq0": "MAPQ 0 (% картированных)",
    "qc_insert": "размер вставки, медиана (п. н.)",
    "qc_error": "частота ошибок (%)",
    "qc_readlen": "средняя длина рида (п. н.)",
    "qc_basequal": "среднее качество оснований",
    "qc_targets": "целевые области в BRCA1/BRCA2",
    "qc_target_mean": "средняя глубина целевых областей (x)",
    "qc_target_pct20": "основания целевых областей >= 20x (%)",
    "qc_uniformity": "равномерность: основания >= 0,2 x среднего (%)",
    "qc_fold80": "fold-80 base penalty",
    "qc_ontarget": "риды на целевых областях (% всех картированных)",
    "qc_weak": "слабо покрытые целевые области",
    "qc_clipped": "риды с обрезанным праймером (%)",
    "qc_short": "риды, слишком короткие после обрезки праймеров",
    "qc_flag": "Предупреждение",
    "weak_title": "Слабо покрытые целевые области (средняя глубина < {depth:g}x)",
    "weak_headers": ["Цель", "Область", "Среднее", "Минимум", "% >=20x"],
    "no_weak": "Все целевые области покрыты достаточно.",
    "no_entries": "Записей нет.",
    "versions": "Версии инструментов",
    "version": "Версия",
    "unknown": "неизвестно",
    "commands": "Воспроизводимые команды",
    "commands_note": "Все выполненные команды полностью приведены в",
    "limitations": "Ограничения",
    "na": "н/д",
}

EXACT_RU = {
    "Mosaik-Kandidat": "кандидат в мозаицизм", "schwach abgedeckt": "слабо покрыто",
    "nicht im Zielbereich des Panels": "не входит в целевые области панели",
    "niedriger Allelanteil: Bestätigung mit orthogonaler Methode nötig":
        "низкая доля аллеля: требуется подтверждение ортогональным методом",
    "Ion Torrent: Insertionen/Deletionen in Homopolymeren unzuverlässig":
        "Ion Torrent: инсерции/делеции в гомополимерах ненадёжны",
    "nicht bestimmt": "не определён", "homozygot Referenz": "гомозигота по референсу",
    "homozygot": "гомозигота", "heterozygot": "гетерозигота",
    "compound/multiallelisch": "компаунд/мультиаллельный",
    "unsicher": "неопределённый", "ok": "ok",
    "Deletion": "Делеция", "Duplikation": "Дупликация",
    "möglicher Verlust": "возможная потеря", "möglicher Zugewinn": "возможный прирост",
    "geringe Tiefe (DP < 15)": "низкая глубина (DP < 15)",
    "geringe Genotypqualität (GQ < 20)": "низкое качество генотипа (GQ < 20)",
    "Allelverhältnis auffällig (mögliches Mosaik oder Artefakt)": "необычное соотношение аллелей (возможен мозаицизм или артефакт)",
    "geringe Mapping-Qualität (MQ < 40)": "низкое качество картирования (MQ < 40)",
    "unpräzise Bruchpunkte": "неточные точки разрыва",
    "Bestätigung mit orthogonaler Methode nötig": "требуется подтверждение ортогональным методом",
    "Kandidat, Bestätigung per MLPA oder validierter CNV-Analyse nötig":
        "кандидат, требуется подтверждение методом MLPA или валидированным анализом CNV",
    "eingefügte Sequenz mit Alu-Signatur: mögliche Alu-Insertion":
        "вставленная последовательность с сигнатурой Alu: возможная вставка Alu",
    "abgeschnittene Leseenden mit Alu-Sequenz: mögliche Alu-Insertion":
        "обрезанные концы прочтений с последовательностью Alu: возможная вставка Alu",
    "abgeschnittene Leseenden von beiden Seiten: mögliche Insertion":
        "обрезанные концы прочтений с обеих сторон: возможная вставка",
    "Hinweis aus Leseenden, kein Variantenaufruf; Bestätigung mit orthogonaler Methode nötig":
        "указание по концам прочтений, не вызов варианта; требуется подтверждение ортогональным методом",
    "Insertionssignatur": "сигнатура вставки",
    "Tandemkopie-Signatur": "сигнатура тандемной копии",
    "zusätzliche 126-bp-Kopie vor BRCA1 Exon 1 (häufiges Allel)": "дополнительная копия 126 п. н. перед экзоном 1 BRCA1 (частый аллель)",
    "bekannte zusätzliche Tandemkopie": "известная дополнительная тандемная копия",
    "Markerbase der zusätzlichen Kopie in den Lesungen vorhanden": "маркерное основание дополнительной копии присутствует в прочтениях",
    "gekoppelter Haplotyp der zusätzlichen Kopie vorhanden": "сцепленный гаплотип дополнительной копии присутствует",
    "bekannte intronische Alu-Insertion in BRCA2 (häufiges Allel)": "известная интронная вставка Alu в BRCA2 (частый аллель)",
    "bekanntes Insertionsallel": "известный аллель вставки",
    "keine Lesung für das Referenzallel: vermutlich homozygot": "нет прочтений референсного аллеля: вероятно, гомозигота",
    "Lesungen für beide Allele: vermutlich heterozygot": "прочтения обоих аллелей: вероятно, гетерозигота",
    "Lesungen passen besser zur bekannten Allelsequenz als zur Referenz (Neuausrichtung)":
        "прочтения лучше соответствуют известной последовательности аллеля, чем референсу (повторное выравнивание)",
    "Haplotyp homozygot: vermutlich zwei Allele mit der zusätzlichen Kopie":
        "гаплотип гомозиготный: вероятно, два аллеля с дополнительной копией",
    "Haplotyp heterozygot: vermutlich ein Allel mit der zusätzlichen Kopie":
        "гаплотип гетерозиготный: вероятно, один аллель с дополнительной копией",
    "abgeschnittene Leseenden an beiden Grenzen einer Wiederholungseinheit: mögliche zusätzliche Tandemkopie":
        "обрезанные концы прочтений на обеих границах повторяющейся единицы: возможная дополнительная тандемная копия",
    "Bruchpunktpaar an derselben Stelle: mögliche Insertion (z. B. mobiles Element wie Alu)":
        "пара точек разрыва в одном месте: возможная вставка (например, мобильного элемента, такого как Alu)",
    "Exonweise Tiefe relativ zum Chromosomen-Median; nur ein Hinweis, kein CNV-Call. Exon ist zugleich als unsicher markiert: Verlust oder Abdeckungsproblem, nicht unterscheidbar.":
        "Глубина экзона относительно медианы хромосомы; только указание, не вызов CNV. Экзон одновременно отмечен как неопределённый: потерю и проблему покрытия различить нельзя.",
    "Exonweise Tiefe relativ zum Chromosomen-Median; nur ein Hinweis, kein CNV-Call.":
        "Глубина экзона относительно медианы хромосомы; только указание, не вызов CNV.",
}

PATTERNS_RU = [
    (re.compile(r"^hoher Duplikatanteil \((.+) %\)$"), r"высокая доля дубликатов (\1 %)"),
    (re.compile(r"^geringer Anteil kartierter Reads \((.+) %\)$"), r"низкая доля картированных ридов (\1 %)"),
    (re.compile(r"^hohe Fehlerrate \((.+) %\)$"), r"высокая частота ошибок (\1 %)"),
    (re.compile(r"^(\d+) Ziel\(e\) mit mittlerer Tiefe < (.+)x$"), r"целевых областей со средней глубиной < \2x: \1"),
    (re.compile(r"^ungleichmäßige Abdeckung \(nur (.+) % der Zielbasen >= 0,2 x Mittel\)$"),
     r"неравномерное покрытие (только \1 % оснований >= 0,2 x среднего)"),
    (re.compile(r"^nur (.+) % der Reads mit abgeschnittenem Primer$"), r"только \1 % ридов с обрезанным праймером"),
    (re.compile(r"^mittlere Tiefe (.+)x < (.+)x$"), r"средняя глубина \1x < \2x"),
    (re.compile(r"^nur (.+) % der Basen >= 20x$"), r"только \1 % оснований >= 20x"),
    (re.compile(r"^(.+) % der Reads mit MAPQ < 20 \(Pseudogen/Repeat\)$"), r"\1 % прочтений с MAPQ < 20 (псевдоген/повтор)"),
    (re.compile(r"^liegt in unsicher abgedecktem Bereich (.+)$"), r"находится в неопределённо покрытой области \1"),
    (re.compile(r"^Exon (\d+) \(Legacy (\d+)\)$"), r"Экзон \1 (устар. нумерация \2)"),
    (re.compile(r"^ENIGMA Tabelle 4 \(falls bestätigt\): (.+)$"), r"ENIGMA, таблица 4 (при подтверждении): \1"),
]

LANGUAGES = ("de", "en", "ru")
_ui_de_en = ui
_tr_de_en = tr


def ui(key, lang, **kwargs):  # noqa: F811  (extends the de/en version above)
    if lang == "ru" and key in UI_RU:
        value = UI_RU[key]
        return value.format(**kwargs) if kwargs and isinstance(value, str) else value
    return _ui_de_en(key, "de" if lang == "de" else "en", **kwargs)


def tr(text, lang):  # noqa: F811
    if lang == "ru" and isinstance(text, str):
        if text in EXACT_RU:
            return EXACT_RU[text]
        for pattern, replacement in PATTERNS_RU:
            if pattern.match(text):
                return pattern.sub(replacement, text)
        return _tr_de_en(text, "en")
    return _tr_de_en(text, lang)
