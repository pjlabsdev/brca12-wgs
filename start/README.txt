brca12-wgs Start
================

BRCA1/BRCA2 analysis from whole-genome, exome or panel data (GRCh38).
Research use only (RUO). Not a medical report and not a medical device. Results must not be
interpreted medically without clinically validated confirmation in an accredited laboratory.

This package contains the start file and the pipeline (folder "pipeline"). It does not contain
your data, the reference genome, the Ensembl VEP cache or ClinVar.

Start
-----
  macOS     double-click "brca12-wgs-start.command" (opens in Terminal).
            If macOS blocks it the first time: right-click, Open.
  Linux     open a terminal in this folder and run:  bash brca12-wgs-start.sh
  Windows   double-click "brca12-wgs-start.bat". It runs the start file in WSL 2.
            Install WSL 2 once in an administrator PowerShell:  wsl --install

Updates: at every start the start file checks GitHub for a newer release (nothing about you
or your data is sent) and asks "Update now? [y/N]". On "y" it replaces its own program files and
starts again; your results and settings stay. BRCA12_NO_UPDATE=1 switches the check off.

The start file asks for a BAM or CRAM (drag the file into the window or type its path), shows
input, reference and output folder and starts after Enter. At the end it opens the report
(brca12_report.html, English and German in one file). The output goes into a new folder
"brca12-<sample>" in your Downloads folder. Nothing is uploaded; everything runs on your computer.

Needed once
-----------
  1. Java 17 or newer and Nextflow 24.10 or newer (https://www.nextflow.io).
     On Windows install both inside WSL.
  2. The tools, one of:
       conda or mamba (default; environments are created at the first run), or
       Docker (set PROFILE=docker in the settings file).
  3. The reference genome the reads were aligned to (GRCh38, FASTA with .fai index).
  4. The Ensembl VEP cache (homo_sapiens, GRCh38) and the ClinVar VCF (clinvar.vcf.gz with index).
The start file asks for 3 and 4 once and remembers them in ~/.config/brca12-wgs/start.conf.
Further Nextflow options go into EXTRA_ARGS in that file.

Language
--------
English. German is used when the system language is German. To set it yourself start with
BRCA12_LANG=en or BRCA12_LANG=de, for example:  BRCA12_LANG=de bash brca12-wgs-start.sh

Documentation, source and citation: https://github.com/pjlabsdev/brca12-wgs
DOI: https://doi.org/10.5281/zenodo.24173638


Deutsch
=======

BRCA1/BRCA2-Auswertung aus Genom-, Exom- oder Paneldaten (GRCh38).
Nur für Forschung (RUO). Kein medizinischer Befund und kein Medizinprodukt. Ergebnisse dürfen
nicht ohne klinisch validierte Bestätigung in einem akkreditierten Labor medizinisch gedeutet werden.

Dieses Paket enthält die Startdatei und die Pipeline (Ordner "pipeline"). Es enthält weder deine
Daten noch das Referenzgenom, den Ensembl-VEP-Cache oder ClinVar.

Start
-----
  macOS     Doppelklick auf "brca12-wgs-start.command" (öffnet im Terminal).
            Falls macOS beim ersten Mal blockiert: Rechtsklick, Öffnen.
  Linux     Terminal in diesem Ordner öffnen und starten:  bash brca12-wgs-start.sh
  Windows   Doppelklick auf "brca12-wgs-start.bat". Die Startdatei läuft in WSL 2.
            WSL 2 einmalig in einer PowerShell als Administrator installieren:  wsl --install

Updates: Bei jedem Start fragt die Startdatei GitHub nach einer neueren Version (über dich und
deine Daten wird nichts gesendet) und fragt "Jetzt aktualisieren? [j/N]". Bei "j" ersetzt sie ihre
Programmdateien und startet neu; Ergebnisse und Einstellungen bleiben. BRCA12_NO_UPDATE=1 schaltet
die Prüfung ab.

Die Startdatei fragt nach einer BAM oder CRAM (Datei ins Fenster ziehen oder Pfad eingeben), zeigt
Eingabe, Referenz und Ausgabeordner und startet nach Enter. Am Ende öffnet sie den Bericht
(brca12_report.html, Englisch und Deutsch in einer Datei). Die Ausgabe liegt in einem neuen Ordner
"brca12-<Probe>" im Downloads-Ordner. Es wird nichts hochgeladen; alles läuft auf deinem Rechner.

Einmalig nötig
--------------
  1. Java 17 oder neuer und Nextflow 24.10 oder neuer (https://www.nextflow.io).
     Unter Windows beides in WSL installieren.
  2. Die Werkzeuge, eines von beiden:
       conda oder mamba (Standard; die Umgebungen entstehen beim ersten Lauf), oder
       Docker (in der Einstellungsdatei PROFILE=docker setzen).
  3. Das Referenzgenom, mit dem die Reads ausgerichtet wurden (GRCh38, FASTA mit .fai-Index).
  4. Den Ensembl-VEP-Cache (homo_sapiens, GRCh38) und die ClinVar-VCF (clinvar.vcf.gz mit Index).
Nach 3 und 4 fragt die Startdatei einmal und merkt sie sich in ~/.config/brca12-wgs/start.conf.
Weitere Nextflow-Optionen stehen dort in EXTRA_ARGS.

Sprache
-------
Englisch. Deutsch erscheint, wenn die Systemsprache Deutsch ist. Zum Festlegen mit
BRCA12_LANG=en oder BRCA12_LANG=de starten, zum Beispiel:  BRCA12_LANG=de bash brca12-wgs-start.sh

Dokumentation, Quelltext und Zitat: https://github.com/pjlabsdev/brca12-wgs
DOI: https://doi.org/10.5281/zenodo.24173638
