[![DOI](https://zenodo.org/badge/623420849.svg)](https://zenodo.org/badge/latestdoi/623420849) [![CC BY-SA 4.0][cc-by-sa-shield]][cc-by-sa]  
This work is licensed under a
[Creative Commons Attribution-ShareAlike 4.0 International License][cc-by-sa].


[cc-by-sa]: http://creativecommons.org/licenses/by-sa/4.0/
[cc-by-sa-image]: https://licensebuttons.net/l/by-sa/4.0/88x31.png
[cc-by-sa-shield]: https://img.shields.io/badge/License-CC%20BY--SA%204.0-lightgrey.svg

# Vorwort

Dies ist das Vorlesungsskript meiner Vorlesung 'Festkörperphysik II'. Sie ist eine Kursvorlesung für  Studierende im dritten Jahr des Bachelorstudiums. Bei der Auswahl und Gewichtung der Themen folgt sie sehr stark dem in Bayreuth Üblichen. 


Dieses Skript ist 'work in progress', und wahrscheinlich nie wirklich fertig.   Wenn Sie Fehler finden, sagen Sie es mir bitte. 
Die aktuellste Version des Vorlesungsskripts finden Sie auf [github](https://github.com/MarkusLippitz/Festkoerper_II). Ich habe alles unter eine CC-BY-SA-Lizenz gestellt. In meinen Worten: Sie können damit machen, was Sie wollen. Wenn Sie Ihre Arbeit der Öffentlichkeit zur Verfügung stellen, erwähnen Sie mich und verwenden Sie eine ähnliche Lizenz. 


Der Text wurde mit der LaTeX-Klasse ['tufte-book'](https://tufte-latex.github.io/tufte-latex/) von Bil Kleb, Bill Wood und Kevin Godby  gesetzt, die sich der Arbeit von [Edward Tufte](https://www.edwardtufte.com/) annähert. Ich habe viele der Modifikationen angewandt, die von Dirk Eddelbüttel im R-Paket ['tint'](https://dirk.eddelbuettel.com/code/tint.html) eingeführt wurden. Die Quelle ist vorerst LaTeX, nicht Markdown.




## Struktur des Repositorys

Die Kapitel sind einzelne TeX-Dateien im Verzeichnisbaum, die per 'include' in die Haupt-TeX-Datei eingefügt werden. Die Abbildungen werden meist mit tikz erstellt. Tikz-external wird verwendet, um Kompilierungszeit zu sparen. Es generiert die pdfs im Verzeichnis 'tikz_external'. 

