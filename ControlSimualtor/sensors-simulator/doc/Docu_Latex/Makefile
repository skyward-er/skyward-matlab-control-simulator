#per poter compilare da riga di comando, sostituite “Rinominare_questo_file.tex” con il nome del vostro file principale, ad esempio, “mioreport.tex”

all:
	latexmk -interaction=batchmode -pdf Rinominare_questo_file.tex

realtime:
	latexmk -interaction=batchmode -pvc -pdf Rinominare_questo_file.tex

clean:
	rm -f *.aux *.fdb_latexmk *.fls *.lof *.log *.lot *.nlo *.out *.pdf *.tdo *.toc
tidy:
	rm -f *.aux *.fdb_latexmk *.fls *.lof *.log *.lot *.nlo *.out *.tdo *.toc

.PHONY: clean realtime

